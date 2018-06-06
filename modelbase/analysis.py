#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 21:30:15 2018

@author: oliver
"""

import scipy.optimize as opt

import numdifftools as nd

from .modelbase import Model

import numpy as np



class Analysis(Model):
    

    def numericElasticities(self, y0, rate):
        '''
        y0: state vector
        rate: name of rate for which elasticities shall be determined
        '''

        def vi(y):
            v = self.rates(y)
            return v[rate]

        jac = nd.Jacobian(vi,step=y0.min()/100)

        epsilon = jac(y0)

        return epsilon

    def allElasticities(self, y0, norm=False):
        '''
        calculates all elasticities:
        :param y0: state vector
        :return: all elasticities as np.matrix
        '''

        rateIds = self.rateNames()

        epsilon = np.zeros([len(rateIds), len(self.cpdNames)])

        for i in range(len(rateIds)):

            def vi(y):
                return self.rateFn[rateIds[i]](y)

            jac = nd.Jacobian(vi, step=y0.min()/100)

            epsilon[i,:] = jac(y0)

        if norm:
            v = np.array(self.rates(y0).values())
            epsilon = (1/v).reshape(len(v),1)*epsilon*y0

        return np.matrix(epsilon)


    def numericJacobian(self, y0, **kwargs):
        '''
        y0: state vector at which Jacobian is calculated
        '''
        J = np.zeros([len(y0),len(y0)])

        if np.isclose(y0.min(),0):
            jstep = None
        else:
            jstep = y0.min()/100

        for i in range(len(y0)):

            def fi(y):
                dydt = self.model(y, 0, **kwargs)
                return dydt[i]

            jac = nd.Jacobian(fi,step=jstep)

            J[i,:] = jac(y0)

        return np.matrix(J)


    def findSteadyState(self, y0, **kwargs):
        '''
        tries to find the steady-state by numerically solving the algebraic system dy/dt = 0.
        input: y0: initial guess
        TODO: this method can be improved. So far, it simply tries the standard solving method hybr
        '''

        def fn(x):
            return self.model(x, 0, **kwargs)
        sol = opt.root(fn, y0)

        if sol.success == True:
            return sol.x
        else:
            return False



    def concentrationControlCoefficients(self, y0, pname, norm=True, **kwargs):
        '''
        invokes findSteadyState to calculate the concentration control coefficients
        for parameter pname
        :input y0: initial guess for steady-state
        :input pname: parameter name to vary
        :input norm: if True (default), normalize coefficients
        :returns: response coefficients
        '''

        origValue = getattr(self.par, pname)

        def fn(x):
            self.par.update({pname: x})
            return self.findSteadyState(y0, **kwargs)

        jac = nd.Jacobian(fn, step=origValue/100.)

        cc = np.array(jac(origValue))

        self.par.update({pname: origValue})

        if norm:
            ss = self.findSteadyState(y0, **kwargs)
            cc = origValue * cc / ss.reshape(ss.shape[0],1)

        return cc
