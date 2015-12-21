__author__ = 'oliver'

from modelbase.model import Model
import numpy as np


class AlgmModel(Model):
    def __init__(self, pars={}, defaultpars={}):
        super(AlgmModel,self).__init__(pars,defaultpars)
        self.algebraicModules = []

    def add_algebraicModule(self, am, amVars, amCpds):
        '''
        this adds a module in which several compound concentrations can be calculated algebraicly.
        am: modelbase.algebraicModule.AlgebraicModule
        amVars: list of names of variables used for module in embedding model
        amCpds: list of names of compounds which are calculated by the module from amVars
        '''
        
        self.algebraicModules.append({'am': am, 'amVars': amVars, 'amCpds': amCpds})


    def set_rate(self, rateName, fn, *args):
        '''
        sets a rate. Arguments:
        Input: rateName (string), fn (the function) and _names_ of compounds which are passed to the function.
        The function fn is called with the parameters self.par as first argument and the dynamic variables corresponding to the compounds as variable argument list.

        In contrast to class Model, args can contain a name from an algebraic module

        '''

        cpdids = {it: id for id, it in enumerate(self.cpdNames)}
        cnt = len(self.cpdNames)

        for ammod in self.algebraicModules:
            cpdids.update({it: id for id, it in enumerate(ammod['amCpds'], cnt)})
            cnt += len(ammod['amCpds'])

        argids = np.array([cpdids[x] for x in args])

        if len(argids) == 0:
            def v(y):
                return fn(self.par)
        else:
            def v(y):
                cpdarg = y[argids]
                return fn(self.par,*cpdarg)

        self.rateFn[rateName] = v
        


    def rates(self, y):

        arglist = [y]
        cpdids = {it: id for id, it in enumerate(self.cpdNames)}

        for ammod in self.algebraicModules:
            varids = np.array([cpdids[x] for x in ammod['amVars']])
            arglist.append(ammod['am'].getConcentrations(y[varids]))

        z = np.hstack(arglist)

        return {r:self.rateFn[r](z) for r in self.stoichiometries.keys()}

