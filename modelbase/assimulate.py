#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 19:54:36 2018

@author: oliver
"""

from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

from .simulate import Simulate
from .simulate import LabelSimulate

import numpy as np


class Assimulate(Simulate):
    
    '''
    class Assimulate works analogously to class Simulate.
    The fundamental integration routines are
    - generate_integrator
    - set_initial_value
    - integrate
    - timeCourse
    All other methods are inherited from super class Simulate
    '''
    
    def __init__(self, model, **kwargs):
        
        #super(Assimulate, self).__init__(model, **kwargs)
        # the simple solution above does not work, because Simulate.__init__
        # calls generate_integrator, but this one already requires the 'f'
        # to be defined.
        
        # a more elegant way would be to have one base class that does 
        # everything except defining f and calling generate_integrator
        
        self.model = model

        def dydt(t, y, m):
            return m.model(y, t, **kwargs)
        def f(t,y):
            return(self.dydt(t,y,self.model))

        self.dydt = dydt
        self.f = f
        self._successful = True
        self._monitor = True
        self._warnings = False
        self.clearResults()
        if "verbosity" in kwargs:
            self.generate_integrator(verbosity=kwargs["verbosity"])
        else:
            self.generate_integrator()
            
    
    def generate_integrator(self, y0=None, name='---',verbosity=50):
 
        if y0 is None:
            y0 = np.zeros(len(self.model.cpdNames))
            
        self.problem = Explicit_Problem(self.f, y0=y0, name=name)
        self.integrator = CVode(self.problem)
        self.integrator.verbosity = verbosity

    
    def set_initial_value(self, y0, t0=0):
        
        self.integrator.y = y0
        self.integrator.t = t0
        

    def integrate(self, t, integrator=None, minstep=None, maxstep=None, nsteps=None):
        
        self._successful = True
        
        try:
            T,Y = self.integrator.simulate(t)
        except:
            print("Error while integrating with CVode")
            self._successful = False
            
        if len(Y.shape) == 1:
            Ylast = Y[-1]
        else:
            Ylast = Y[-1,:]
        return Ylast
    
    
    def timeCourse(self, Torig, y0, integrator=None, minstep=None, maxstep=None, nsteps=None):
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

        try:
            t,Y = self.integrator.simulate(T[-1],ncp_list=T)
        except:
            print("Error in timeCourse while integrating with CVode")
            self._successful = False
            
        if self.doesMonitor() and self.successful():
            self.results.append({'t': T, 'y': Y})

        return np.vstack(Y)
    
    


class LabelAssimulate(Assimulate,LabelSimulate):
    
    pass

