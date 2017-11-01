# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 15:28:33 2015

@author: oliver
"""

import numpy as np
import scipy.integrate as sci
import math

import itertools

#import numdifftools as nd

import pickle

class Simulate(object):

    def __init__(self, model, **kwargs):

        self.model = model

        def dydt(t, y, m):
            return m.model(y, t, **kwargs)

        self.dydt = dydt
        self._successful = True
        self._monitor = True
        self._warnings = False
        self.clearResults()
        self.generate_integrator()

    def successful(self):
        return self._successful

    def doesMonitor(self, setMonitor=None):
        if setMonitor != None:
            self._monitor = setMonitor
        return self._monitor

    def clearResults(self):
        self.results = []

    def storeResults(self, filename):
        '''
        stores the parameters to file FILENAME
        :input filename: FILENAME
        '''
        f = open(filename,'wb')
        pickle.dump(self.results, f)
        f.close

    def loadResults(self, filename):
        '''
        loads results from file and stores in results attribute
        NOTE: overrides old results
        '''
        res = pickle.load(open(filename,'rb'))
        self.results = res


    def generate_integrator(self, integrator='lsoda', max_step=0.1, nsteps=500):
        '''
        generates a sci.ode object used for integration
        '''
        self.integrator = sci.ode(self.dydt).set_integrator(integrator, max_step=max_step, nsteps=nsteps)
        self._integrator = integrator
        self._max_step = max_step
        self._nsteps = nsteps

    def set_initial_value(self, y0, t0=0):
        '''
        initialises the sci.ode integrator
        '''
        self.integrator.set_initial_value(y0, t0)
        self.integrator.set_f_params(self.model)

    def set_initial_value_to_last(self):
        '''
        initialises the sci.ode integrator to last values stored in results
        '''
        tlast = self.getT()[-1]
        ylast = self.getVarsByName(self.model.cpdNames)[-1,:]
        self.set_initial_value(ylast, tlast)


    def integrate(self, t, minstep=1e-8, maxstep=0.1, nsteps=500):
        """ integration, returns variables at time t """

        #if self.integrator is None:
        #    raise(ModelError, "no integrator defined. Call generate_integrator")

        r = self.integrator

        y0 = r.y
        t0 = r.t

        step = maxstep
        numsteps = max(nsteps, 10*math.floor((t-t0)/step))

        while step >= minstep:

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

            r.set_integrator(self._integrator, max_step=step, nsteps=numsteps)
            r.set_initial_value(y0, t0)
            r.set_f_params(self.model)


        if step < maxstep: # set back to standard steps
            r.set_integrator(self._integrator, max_step=self._max_step, nsteps=self._nsteps)


        self._successful = r.successful()
        return r.y


    def timeCourse(self, Torig, y0, integrator='lsoda', minstep=1e-8, maxstep=0.1, nsteps=500):
        """ integration over time, different integrators possible, lsoda default
            returns: array of state variables
        """

        self._successful = True

        T = Torig.copy()

        if y0 is not None:
            Y = [y0]
            #print Y, type(Y)

            self.set_initial_value(y0,t0=T[0])
        else:
            Y = [np.array(self.integrator.y)]
            if T[0] == 0:
                T += self.integrator.t

        cnt = 1
        while cnt < len(T) and self.successful():
            if self._warnings:
                print(cnt, Y)
                print(T[cnt])
            Y.append(self.integrate(T[cnt],
                                    minstep=minstep,
                                    maxstep=maxstep,
                                    nsteps=nsteps))
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

        self.set_initial_value(Y0, t0=T)

        while self.successful() and cnt < maxstep and err > AbsTol:
            Y = self.integrate(T+step)
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


    def estimatePeriod(self,y0,t0=0.,twait=1000.,tend=3000.,dt=0.1,osctol=1.,varno=0):
        '''
        attempts to estimate a period from a simulation running to a stable limit cycle.
        This does only work for 'smooth' oscillation. 
        TODO: a more robust method based on autocorrelations (see np.convolve)
        '''
        T = np.arange(t0,tend,dt)
        
        Y = self.timeCourse(T,y0)

        iwait = np.where(T>twait)[0].min()

        if not self.successful():
            return False, False, False

        A = Y[iwait:,varno]
        Amax = A.max()
        Amin = A.min()

        if Amax - Amin < osctol:
            return False, False, False

        m = (Amax+Amin)/2.

        if A[-1] < m:
            A = Amax + Amin - A

        i0 = np.where(A<m)[0].max()
        w1 = np.where(A[:i0]>m)[0]
        if len(w1) == 0:
            return False, False, False
        else:
            i1 = w1.max()
        w2 = np.where(A[:i1]<m)[0]
        if len(w2) == 0:
            return False, False, False
        else:
            i2 = w2.max()

        P = T[iwait+i0] - T[iwait+i2]

        ymax = Y[iwait+i2:iwait+i0,:].max(0)
        ymin = Y[iwait+i2:iwait+i0,:].min(0)

        return P, ymax, ymin

        


    # What follows are the old Results methods.
    # It does not make sense to separate Simulate and Results

    def getT(self, r=None):
        """
        :param r: range of steps of simulation for which results we are interested in
        :return: time of all results as one vector
        """
        if r is None:
            r = range(len(self.results))

        T = np.hstack([self.results[i]['t'] for i in r])

        return T

    def getY(self, r=None):
        """
        :param r: range of steps of simulation for which results we are interested in
        :return: values of all variables as one array
        """
        if r is None:
            r = range(len(self.results))

        Y = np.vstack([self.results[i]['y'] for i in r])
        return Y


    def getVar(self, j, r=None):
        """
        :param j: int of the state variable [0:PQred, ..., 5:ATPsynth]
        :param r: range of steps of simulation for which results we are interested in
        :return: concentrations of variable j as one vector
        """

        if type(j) == int:
            j = [j]

        if r is None:
            r = range(len(self.results))

        Y = self.getY(r)
        X = Y[:,j]

        #X = np.vstack([self.results[i]['y'][:,j] for i in r])
        #X = np.vstack([np.reshape(self.results[i]['y'][:, j],np.size(self.results[i]['y'],0),len(j)) for i in r])

        if np.size(X,1) == 1:
            X = np.reshape(X,[np.size(X,0)])

        return X

    def getVarByName(self, cpdName, r=None):

        if r is None:
            r = range(len(self.results))

        ids = self.model.get_argids(cpdName)

        return self.getVar(ids, r)

    def getVarsByName(self, cpdNames, r=None):

        if r is None:
            r = range(len(self.results))

        ids = self.model.get_argids(*cpdNames)

        return self.getVar(ids, r)

    def getVarsByRegexp(self, regexp, r=None):

        if r is None:
            r = range(len(self.results))

        ids = self.model.find_re_argids(regexp)

        return self.getVar(ids, r)



    def getV(self, r=None):

        if r is None:
            r = range(len(self.results))

        V = np.array([])

        for i in r:
            if not self.results[i].has_key('v'):
                
                t = self.results[i]['t']
                y = self.results[i]['y']
                rlist = self.model.rateNames()

                Vnew = []
                for j in range(len(t)):
                    vd = self.model.rates(y[j])
                    vt = np.array([vd[k] for k in rlist])
                    Vnew.append(vt)

                self.results[i]['v'] = np.vstack(Vnew)

        V = np.vstack([self.results[i]['v'] for i in r])
        return V

                       
    def getRate(self, rate, r=None):
        """
        :param rate: name of the rate
        :param r: range of steps of simulation for which results we are interested in
        :return: rate with name 'rate' of all results as one vector
        """
        if r is None:
            r = range(len(self.results))

        V = np.array([])

        for i in r:
            t = self.results[i]['t']
            y = self.results[i]['y']

            V = np.hstack([V,
                           np.array(
                           [self.model.rates(y[j])[rate] for j in range(len(t))]
                           )])

        return V



         

class AlgmSimulate(Simulate):
    '''
    subclass of Simulate, including access to derived Variables
    '''

    def getY(self, r=None):
        """
        :param r: range of steps of simulation for which results we are interested in
        :return: values of all variables (including algModel derived) as one array
        """
        Y = super(AlgmSimulate,self).getY(r)
        Yfull = self.model.fullConcVec(Y)
        return Yfull






class LabelSimulate(AlgmSimulate):
    '''
    subclass of Simulate, including several access methods for labels
    '''

    def getTotal(self, cpdBaseName, r=None):
        '''
        retrieves total concentration for compound cpdBaseName
        :param cpdBaseName: base name of compound
        :return: vector with concentrations
        '''

        if r is None:
            r = range(len(self.results))
        
        c = self.model.cpdBaseNames[cpdBaseName]

        regexp = "\A" + cpdBaseName + '.' * c + "\Z"

        Y = self.getVarsByRegexp(regexp, r)

        return Y.sum(1)

    def getLabelAtPos(self, cpdBaseName, lpos, r=None):
        '''
        retrieves total of cpdBaseName where label is at lpos
        :param cpdBaseName: base name of compound
        :param lpos: position of label
        :return: vector with concentrations
        '''

        if r is None:
            r = range(len(self.results))
        
        if type(lpos) == int:
            lpos = [lpos]

        c = self.model.cpdBaseNames[cpdBaseName]

        l = ['.'] * c
        for p in lpos:
            l[p] = '1'

        regexp = "\A" + cpdBaseName + ''.join(l) + "\Z"

        Y = self.getVarsByRegexp(regexp, r)

        return Y.sum(1)

    def getNumLabel(self, cpdBaseName, nlab, r=None):
        '''
        retrieves total of cpdBaseName with exactly nlab labels
        :param cpdBaseName: base name of compound
        :param nlab: numbers of labels
        :return: vector with concentrations
        '''

        if r is None:
            r = range(len(self.results))

        c = self.model.cpdBaseNames[cpdBaseName]
        
        lcom = itertools.combinations(range(c),nlab)

        cpdNames = []
        for i in lcom:
            l = ['0'] * c
            for p in i:
                l[p] = '1'
            cpdNames.append(cpdBaseName + ''.join(l))

        Y = self.getVarsByName(cpdNames, r)

        return Y.sum(1)


    def getTotalLabel(self, cpdBaseName, r=None):
        '''
        retrieves total labels of cpdBaseName
        :param cpdBaseName: base name of compound
        :return: vector with concentrations
        '''

        if r is None:
            r = range(len(self.results))

        c = self.model.cpdBaseNames[cpdBaseName]
        
        Ylab = []
        for lnum in range(1,c):
            Ylab.append(self.getNumLabel(cpdBaseName, lnum) * lnum)

        return np.vstack(Ylab).sum(0)

    
    def getTotalRate(self, rateBaseName, r=None):
        '''
        retrieves the sum of all rates starting with 'rateName'
        :rateBaseName: rateBaseName
        :return: rate
        '''
        if r is None:
            r = range(len(self.results))

        rid = {v:k for k,v in enumerate(self.model.rateNames())}
        rsel = []
        for k,v in rid.items():
            if k.startswith(rateBaseName):
                rsel.append(v)

        V = self.getV()[:,rsel].sum(1)

        return V
            

        
        
   
    # these two do not belong here, should be part of model.py
    # they have been introduced in model.py but kept here for compatilibity reasons

    #def numericElasticities(self, y0, rate):
    #    '''
    #    y0: state vector
    #    rate: name of rate for which elasticities shall be determined
    #    '''
    #
    #    v0 = self.model.rates(y0)
    #
    #    def vi(y):
    #        v = self.model.rates(y)
    #        return v[rate]
    #
    #    jac = nd.Jacobian(vi,step=y0.min()/100)
    #
    #    epsilon = jac(y0)
    #
    #    return epsilon
    #
    #def numericJacobian(self, y0):
    #
    #    J = np.zeros([len(y0),len(y0)])
    #
    #    for i in range(len(y0)):
    #
    #        def fi(y):
    #            dydt = self.model.model(y)
    #            return dydt[i]
    #
    #        jac = nd.Jacobian(fi,step=y0.min()/100)
    #
    #        J[i,:] = jac(y0)
    #
    #    return np.matrix(J)

        


class ModelError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
