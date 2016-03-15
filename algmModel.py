__author__ = 'oliver'

from modelbase.model import Model
import numpy as np

import numdifftools as nd


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


    
    def get_argids(self, *args):

        cpdids = {it: id for id, it in enumerate(self.cpdNames)}
        cnt = len(self.cpdNames)

        for ammod in self.algebraicModules:
            cpdids.update({it: id for id, it in enumerate(ammod['amCpds'], cnt)})
            cnt += len(ammod['amCpds'])

        return np.array([cpdids[x] for x in args])

    """
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
    """    

    def fullConcVec(self, y):
        '''
        returns the full concentration vector, including all concentrations from algebraic modules
        input: y - state vector of all dynamic variables
        output: z - state vector extended by all derived concentrations
        '''
        z = y.copy()
        #vlist = [y]
        cpdids = {it: id for id, it in enumerate(self.cpdNames)}

        for ammod in self.algebraicModules:
            varids = np.array([cpdids[x] for x in ammod['amVars']])
            zam = ammod['am'].getConcentrations(z[varids])

            cpdidsam = {it:id for id,it in enumerate(ammod['amCpds'], z.size)}

            z = np.hstack([z,zam])
            cpdids = dict(cpdids, **cpdidsam)
            #vlist.append(ammod['am'].getConcentrations(y[varids]))

        #z = np.hstack(vlist)

        return z


    def rates(self, y, **kwargs):

        z = self.fullConcVec(y)

        return {r:self.rateFn[r](z, **kwargs) for r in self.stoichiometries.keys()}



    def allCpdNames(self):
        ''' returns list of all compounds, including from algebraic modules '''
        names = []
        names.extend(self.cpdNames)
        for ammod in self.algebraicModules:
            names.extend(ammod['amCpds'])

        return names

 
    def allElasticities(self, y0):
        ''' 
        calculates all _direct_ elasticities:
        Rates usually depend on a concentration and not directly on a conserved equilbrium module variable.
        Therefore, the partial derivatives of the rate expression itself is zero wrt the equilibrium variable, but non-zero wrt to the concentration.
        Input: y0 - state vector
        '''

        rateIds = self.rateNames()

        epsilon = np.zeros([len(rateIds), len(self.allCpdNames())])

        z0 = self.fullConcVec(y0)

        for i in range(len(rateIds)):

            def vi(y):
                return self.rateFn[rateIds[i]](y)
                
            jac = nd.Jacobian(vi, step=z0.min()/100)

            epsilon[i,:] = jac(z0)

        return np.matrix(epsilon)

