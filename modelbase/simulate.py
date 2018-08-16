# This file is part of modelbase.
#
# modelbase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# modelbase is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with modelbase.  If not, see <http://www.gnu.org/licenses/>.

"""Simulate

Description of the module

"""
import numpy as np
import scipy.integrate as sci
import math
import itertools
import pickle

class Simulate(object):
    """Provides simulation routines

    Attributes
    ----------
    model : class
        modelbase.Model or modelbase.LabelModel instance
    **kwargs : dict
        User defined dynamic changes for rate functions

    Methods
    -------
    generate_integrator(integrator, max_step, nsteps)
        Generates a scipy.integrate.ode object
    set_initial_value(y0, t0)
        Sets initial values for the integration
    set_initial_value_to_last()
        Initialises the sci.ode integrator to last values stored in results
    integrate(t, minstep, maxstep, nsteps)
        Integrates model. Returns values at time point t
    timeCourse(Torig, y0, integrator, minstep, maxstep, nsteps)
        Numerical integration in given time range
    sim2SteadyState(y0, AbsTol, t0, step, maxstep)
        Simulation until numerical approximated steady state
    estimatePeriod(y0, t0, twait, tend, dt, osctol, varno)
        Estimates a period from a simulation running to a stable limit cycle
    getT
        Returns integration time vector over all simulations
    getY
        Returns integration values over all simulations with derived variables
    getVar
        Returns state variables over all simulations
    getVarByName
        Get compound concentrations by compound name
    getVarsByName
        Get compound concentrations by compound names
    getVarsByRegexp
        Get compound concentrations by regular expression
    getV
        Returns rate vector for simulation results
    getRate
        Returns rate vector for the simulation
    """

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
        """Returns integration success state

        Returns
        -------
        successul : bool
            Integration success state
        """
        return self._successful

    def doesMonitor(self, setMonitor=None):
        """Sets monitoring state

        Parameters
        ----------
        setMonitor : bool
            Monitoring state
        Returns
        -------
        monitor : bool
            Monitoring state
        """
        if setMonitor != None:
            self._monitor = setMonitor
        return self._monitor

    def clearResults(self):
        """Clears results array

        Returns
        -------
        None
        """
        self.results = []

    def storeResults(self, filename):
        """Stores integration results in pickle file

        Parameters
        ----------
        filename : str
            Name of the pickle file

        Returns
        -------
        None
        """
        f = open(filename,'wb')
        pickle.dump(self.results, f)
        f.close

    def loadResults(self, filename):
        """Loads results from a pickle file into the class

        Parameters
        ----------
        filename : str
            Name of pickle file

        Returns
        -------
        None
        """
        '''
        loads results from file and stores in results attribute
        NOTE: overrides old results
        '''
        res = pickle.load(open(filename,'rb'))
        self.results = res

    def generate_integrator(self, integrator='lsoda', max_step=0.1, nsteps=500):
        """Generates a scipy.integrate.ode object

        Parameters
        ----------
        integrator : str
            Sets the scipy integrator. Defaults to 'lsoda'
        max_step : float
            Maximal integration step size
        nsteps : int
            Number of integration steps

        Returns
        -------
        None
        """
        self.integrator = sci.ode(self.dydt).set_integrator(integrator, max_step=max_step, nsteps=nsteps)
        self._integrator = integrator
        self._max_step = max_step
        self._nsteps = nsteps

    def set_initial_value(self, y0, t0=0):
        """Sets initial values for the integration.

        Parameters
        ----------
        y0 : numpy.array or list
            Initial values
        t0 : int or float
            Initial time point. Defaults to 0.

        Returns
        -------
        None
        """
        self.integrator.set_initial_value(y0, t0)
        self.integrator.set_f_params(self.model)

    def set_initial_value_to_last(self):
        """Initialises the sci.ode integrator to last values stored in results

        Returns
        -------
        None
        """
        tlast = self.getT()[-1]
        ylast = self.getVarsByName(self.model.cpdNames)[-1,:]
        self.set_initial_value(ylast, tlast)


    def integrate(self, t, minstep=1e-8, maxstep=0.1, nsteps=500):
        """Integrates model. Returns values at time point t.

        Parameters
        ----------
        t : int or float
            Integration end point
        minstep : float
            Minimal step size
        maxstep : float
            Maximal step size
        nsteps : int
            Number of steps

        Returns
        -------
        y : list
            integration values
        """

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
        """Integration over time vector.

        Parameters
        ----------
        Torig : list or numpy.array
            Time array
        y0 : list or numpy.array
            Initial values
        integrator : str
            Scipy integrator. Defaults to "lsoda"
        minstep : float
            Minimal step size. Defaults to 1e-8
        maxstep : float
            Maximal step size. Defaults to 0.1
        nsteps : int
            Number of integration steps

        Returns
        -------
        Y : numpy.array
            Array of integrated values
        """
        self._successful = True
        T = Torig.copy()
        if y0 is not None:
            Y = [y0]
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
        """Simulation until numerical approximated steady state.

        Parameters
        ----------
        y0 : list or numpy.array
            Initial concentrations
        AbsTol : float
            Minimal difference between two simulation steps required
            as steady state cut-off. Defaults to 1e-7
        t0 : int or float
            Integration starting time point. Defaults to 0
        step : float
            Integration step size
        maxstep : int
            Maximal number of integration steps

        Returns
        -------
        Y : list
            Steady state concentrations
        """
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
        """Estimates a period from a simulation running to a stable limit cycle.
        Only works for a 'smooth' oscillation

        Parameters
        ----------
        y0 : list or numpy.array
            Initial concentrations
        t0 : int or float
            Integration time point. Defaults to 0.
        twait : int or float
            TODO - Fill in
        tend : int or float
            TODO - Fill in
        dt : float
            TODO - Fill in
        osctol : float
            TODO - Fill in
        varno : float
            TODO - Fill in
        Returns
        -------
        P : type
            TODO - Fill in
        ymax : float
            TODO - Fill in
        ymin : float
            TODO - Fill in
        """
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

    def getT(self, r=None):
        """Returns integration time vector over all simulations.

        Parameters
        ----------
        r : list
            Range of simulations

        Returns
        -------
        T : numpy.array
            Time of all results as one vector
        """
        if r is None:
            r = range(len(self.results))
        T = np.hstack([self.results[i]['t'] for i in r])
        return T

    def getY(self, r=None):
        """Returns integration values over all simulations with derived variables.

        Parameters
        ----------
        r : list
            Range of simulations

        Returns
        -------
        Y : numpy.array
            Integration values over all simulations with derived variables.
        """
        if r is None:
            r = range(len(self.results))
        Y = np.vstack([self.results[i]['y'] for i in r])
        return self.model.fullConcVec(Y)

    def getVar(self, j, r=None):
        """Returns state variables over all simulations

        Parameters
        ----------
        j : int or list
            IDs of the state variables
        r : list
            Range of simulations

        Returns
        -------
        X : numpy.array
            Concentrations of variables j
        """
        if type(j) == int:
            j = [j]
        if r is None:
            r = range(len(self.results))
        Y = self.getY(r)
        X = Y[:,j]
        if np.size(X,1) == 1:
            X = np.reshape(X,[np.size(X,0)])
        return X

    def getVarByName(self, cpdName, r=None):
        """Get state variable (compound) concentrations
        by compound name.

        Parameters
        ----------
        cpdName : str
            Name of the compound
        r : list
            Range of simulations

        Returns
        -------
        X : numpy.array
            Concentration vector
        """
        if r is None:
            r = range(len(self.results))
        ids = self.model.get_argids(cpdName)
        return self.getVar(ids, r)

    def getVarsByName(self, cpdNames, r=None):
        """Get state variables (compounds) concentrations
        by compound names.

        Parameters
        ----------
        cpdNames : list
            Names of the compound
        r : list
            Range of simulations

        Returns
        -------
        X : numpy.array
            Concentration vector
        """
        if r is None:
            r = range(len(self.results))
        ids = self.model.get_argids(*cpdNames)
        return self.getVar(ids, r)

    def getVarsByRegexp(self, regexp, r=None):
        """Get state variables (compounds) concentrations
        by regular expression.

        Parameters
        ----------
        regexp : str
            Compound name regular expression
        r : list
            Range of simulations

        Returns
        -------
        X : numpy.array
            Concentration vector
        """
        if r is None:
            r = range(len(self.results))
        ids = self.model.find_re_argids(regexp)
        return self.getVar(ids, r)



    def getV(self, r=None):
        """Returns rate vector for simulation results.

        Parameters
        ----------
        r : list
            Range of simulations

        Returns
        -------
        V : numpy.array
            Rate vector
        """
        if r is None:
            r = range(len(self.results))
        V = np.array([])

        for i in r:
            if not 'v' in self.results[i]:
                t = self.results[i]['t']
                y = self.results[i]['y']
                rlist = self.model.rateNames()
                Vnew = []
                for j in range(len(t)):
                    vd = self.model.rates(y[j], **{'t':t[j]})
                    vt = np.array([vd[k] for k in rlist])
                    Vnew.append(vt)
                self.results[i]['v'] = np.vstack(Vnew)
        V = np.vstack([self.results[i]['v'] for i in r])
        return V

    def getRate(self, rate, r=None):
        """Returns rate vector for the simulation.

        Parameters
        ----------
        rate : str
            Name of the rate function
        r : list
            Range of simulations

        Returns
        -------
        V : numpy.array
            Rate vector
        """
        if r is None:
            r = range(len(self.results))
        V = np.array([])
        for i in r:
            t = self.results[i]['t']
            y = self.results[i]['y']
            V = np.hstack([V,
                           np.array(
                           [self.model.rates(y[j], **{'t':t[j]})[rate] for j in range(len(t))]
                           )])
        return V


class LabelSimulate(Simulate):
    """Extends Simulate class with methods for label models.

    Methods
    -------
    getTotal(cpdBaseName, r)
        Retrieve sum of all labeled variants of a base compound
    getLabelAtPos(cpdBaseName, lpos, r)
        Retrieves total concentration of a compound where the given position is labeled
    getNumLabel(cpdBaseName, nlab, r)
        Retrieves total concentration of a compound with the given number of labeles
    getTotalLabel(cpdBaseName, r)
        Retrieves total concentration of a compound that is labeled
    getTotalRate(rateBaseName, r)
        Retrieves the sum of all rates starting with 'rateBaseName'
    """


    def getTotal(self, cpdBaseName, r=None):
        """Retrieve sum of all labeled variants of a base compound

        Parameters
        ----------
        cpdBaseName : str
            Name of the compound
        r : list
            Range of simulations

        Returns
        -------
        Y : numpy.array
            Sum of all labeled variants of a base compound
        """
        if r is None:
            r = range(len(self.results))
        c = self.model.cpdBaseNames[cpdBaseName]
        regexp = "\A" + cpdBaseName + '.' * c + "\Z"
        Y = self.getVarsByRegexp(regexp, r)
        return Y.sum(1)

    def getLabelAtPos(self, cpdBaseName, lpos, r=None):
        """Retrieves total concentration of a compound
        where the given position is labeled.

        Parameters
        ----------
        cpdBaseName : str
            Name of the compound
        lpos : int or list
            Labeled position(s)
        r : list
            Range of simulations

        Returns
        -------
        Y : numpy.array
            Total concentration
        """
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
        """Retrieves total concentration of a compound
        with the given number of labeles.

        Parameters
        ----------
        cpdBaseName : str
            Name of the compound
        nlab : int
            Number of labeled positions
        r : list
            Range of simulations

        Returns
        -------
        Y : numpy.array
            Total concentration
        """
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
        if len(Y.shape) > 1:
            return Y.sum(1)
        else:
            return Y


    def getTotalLabel(self, cpdBaseName, r=None):
        """Retrieves total concentration of a compound
        that is labeled.

        Parameters
        ----------
        cpdBaseName : str
            Name of the compound
        r : list
            Range of simulations

        Returns
        -------
        Y : numpy.array
            Total concentration
        """
        if r is None:
            r = range(len(self.results))
        c = self.model.cpdBaseNames[cpdBaseName]
        Ylab = []
        for lnum in range(1,c+1):
            Ylab.append(self.getNumLabel(cpdBaseName, lnum) * lnum)
        return np.vstack(Ylab).sum(0)


    def getTotalRate(self, rateBaseName, r=None):
        """Retrieves the sum of all rates starting with 'rateBaseName'

        Parameters
        ----------
        rateBaseName : str
            Name of the rate
        r : list
            Range of simulations

        Returns
        -------
        V : numpy.array
            Sum of rates
        """
        if r is None:
            r = range(len(self.results))
        rid = {v:k for k,v in enumerate(self.model.rateNames())}
        rsel = []
        for k,v in rid.items():
            if k.startswith(rateBaseName):
                rsel.append(v)
        V = self.getV()[:,rsel].sum(1)
        return V


class ModelError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
