__author__ = 'oliver'

import modelbase.parameters as param
import numpy as np

import numdifftools as nd

class AlgebraicModule(object):

    '''
    This class represents a module of a larger model. This module is characterised by
    the fact that several concentrations can be algebraicly calculated from a small number
    of variables. Examples are rapid equilibrium or quasi steady-state modules.

    
    '''

    def __init__(self,par,fn=None):
        '''
        initiation by passing parameters, a function and two lists:
        fn: function calculating concentration of compounds in module from the values of the variables describing the modules
        fn must accept two arguments fn(par,y). par: ParameterSet; y: np.array
        # amVars: names of variables that are dynamic in the embedding model,
        # amCpds: names of compounds that can be calculated from it
        '''

        self.par = param.ParameterSet(par)

        #self.varNames = amVars # actually not necessary here. Only when including in embedding model
        #self.cpdNames = amCpds

        if fn:
            self.convertFun = fn
        else:
            def dummy(x):
                return np.array([0.])
            self.convertFun = dummy


    def getConcentrations(self, y):

        if len(y.shape) == 1:
            return self.convertFun(self.par,y)
        
        else:
            return np.array([self.convertFun(self.par,y[i,:]) for i in range(y.shape[0])])

   

    def elasticities(self, y0):
        '''
        returns the derivatives of all compounds in the algebraic module
        versus the determining quantities at point y0. 
        Returned as np.matrix
        '''

        conc = self.getConcentrations(y0)

        theta = np.zeros([len(conc),len(y0)])

        for i in range(len(conc)):

            def si(y):
                conc = self.getConcentrations(y)
                return conc[i]

            jac = nd.Jacobian(si, step=y0.min()/100)

            theta[i,:] = jac(y0)

        return np.matrix(theta)

