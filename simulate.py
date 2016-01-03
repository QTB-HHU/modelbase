# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 15:28:33 2015

@author: oliver
"""

import numpy as np
import scipy.integrate as sci
import math

import numdifftools as nd

class Simulate:

    def __init__(self, model):

        self.model = model

        def dydt(t, y, m):
            return m.model(y)

        self.dydt = dydt
        self._successful = True
        self._monitor = True
        self._warnings = False
        self.clearResults()

    def successful(self):
        return self._successful

    def doesMonitor(self, setMonitor=None):
        if setMonitor != None:
            self._monitor = setMonitor
        return self._monitor

    def clearResults(self):
        self.results = []

    def integrate(self, t, y0, integrator='lsoda', minstep=1e-8, maxstep=0.1, nsteps=500, t0=0):
        """ integration, returns variables at time t """
        step = maxstep
        numsteps = max(nsteps, 10*math.floor((t-t0)/step))

        while step >= minstep:
            r = sci.ode(self.dydt).set_integrator(integrator, max_step=step, nsteps=numsteps)
            r.set_initial_value(y0, t0)
            r.set_f_params(self.model)

            # suppress FORTRAN warnings
            if not self._warnings:
                r._integrator.iwork[2] = -1
            try:
                r.integrate(t)
                if r.successful():
                    break
            except ModelError:
                print('caught error at ',step,'. Reducing step size')

            step = step/10
            numsteps = numsteps*10

            if self._warnings:
                print('numsteps=', numsteps, ', step=', step)
                print(r.t, r.y)
                print(self.model.rates(r.y))

        self._successful = r.successful()
        return r.y

    def timeCourse(self, T, y0, integrator='lsoda', minstep=1e-8, maxstep=0.1, nsteps=500):
        """ integration over time, different integrators possible, lsoda default
            returns: array of state variables
        """

        self._successful = True

        Y = [y0]
        #print Y, type(Y)
        cnt = 1
        while cnt < len(T) and self.successful():
            if self._warnings:
                print cnt, Y
                print(T[cnt])
            Y.append(self.integrate(T[cnt], Y[cnt-1], t0=T[cnt-1],
                                    minstep=minstep,
                                    maxstep=maxstep,
                                    nsteps=nsteps,
                                    integrator=integrator))
            cnt += 1

        if self.doesMonitor() and self.successful():
            self.results.append({'t': T, 'y': np.vstack(Y)})

        return np.vstack(Y)


    def sim2SteadyState(self, y0, AbsTol = 1e-7, t0 = 0, step = 0.1, maxstep = 1000):
        '''
        Simulation until steady-state is reached (difference between two simulation steps < AbsTol) 
        or maxstep steps have been performed.
        Returns the last value of simulation.
        If unsuccessful (> maxstep simulation steps), self.successful() will return False, else True
        '''
        self._successful = True

        T = t0
        cnt = 0
        Y0 = y0
        err = np.linalg.norm(Y0,ord=2)

        while self.successful() and cnt < maxstep and err > AbsTol:
            Y = self.integrate(T+step, Y0, t0=T)
            T += step
            cnt += 1
            err = np.linalg.norm(Y-Y0, ord=2)
            if self._warnings:
                print('T=', T, ' err=', err)
            Y0 = Y

        if cnt >= maxstep:
            self._successful = False

        elif self.doesMonitor() and self.successful():
            self.results.append({'t':np.array([T]),'y':np.vstack([Y])})

        return Y



    # these two do not belong here, should be part of model.py
    # they have been introduced in model.py but kept here for compatilibity reasons

    def numericElasticities(self, y0, rate):
        '''
        y0: state vector
        rate: name of rate for which elasticities shall be determined
        '''
 
        v0 = self.model.rates(y0)

        def vi(y):
            v = self.model.rates(y)
            return v[rate]

        jac = nd.Jacobian(vi,step=y0.min()/100)

        epsilon = jac(y0)

        return epsilon

    def numericJacobian(self, y0):

        J = np.zeros([len(y0),len(y0)])

        for i in range(len(y0)):

            def fi(y):
                dydt = self.model.model(y)
                return dydt[i]

            jac = nd.Jacobian(fi,step=y0.min()/100)

            J[i,:] = jac(y0)

        return np.matrix(J)

        
