# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 15:28:33 2015

@author: oliver
"""

import numpy as np
import scipy.integrate as sci
import math

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
        print Y, type(Y)
        cnt = 1
        while cnt < len(T) and self.successful():
            print cnt, Y
            Y.append(self.integrate(T[cnt], Y[cnt-1], t0=T[cnt-1],
                                    minstep=minstep,
                                    maxstep=maxstep,
                                    nsteps=nsteps,
                                    integrator=integrator))
            print(T[cnt])
            cnt += 1

        if self.doesMonitor() and self.successful():
            self.results.append({'t': T, 'y': np.vstack(Y)})

        return np.vstack(Y)
