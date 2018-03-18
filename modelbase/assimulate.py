#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 19:54:36 2018

@author: oliver
"""

from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

from .simulate import Simulate

import numpy as np


class Assimulate(Simulate):
    
    def __init__(self, model, **kwargs):
        
        #super(Assimulate, self).__init__(model, **kwargs)
        
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
        self.generate_integrator()
            
    
    def generate_integrator(self, y0=None, name='---'):
 
        if y0 is None:
            y0 = np.zeros(len(self.model.cpdNames))
            
        self.problem = Explicit_Problem(self.f, y0=y0, name=name)
        self.integrator = CVode(self.problem)

    
    def set_initial_value(self, y0, t0=0):
        
        self.integrator.y = y0
        self.integrator.t = t0
        

    def integrate(self, t):
        
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
    
    
    
    
