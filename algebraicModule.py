__author__ = 'oliver'

import modelbase.parameters as param
import numpy as np

class AlgebraicModule(object):

    '''
    This class represents a module of a larger model. This module is characterised by
    the fact that several concentrations can be algebraicly calculated from a small number
    of variables. Examples are rapid equilibrium or quasi steady-state modules.

    
    '''

    def __init__(self,par,fn):
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

        self.convertFun = fn



    def getConcentrations(self, y):

        if len(y.shape) == 1:
            return self.convertFun(self.par,y)
        
        else:
            return np.array([self.convertFun(self.par,y[i,:]) for i in range(y.shape[0])])

   

