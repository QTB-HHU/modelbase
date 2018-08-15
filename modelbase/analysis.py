import scipy.optimize as opt
import numdifftools as nd
import numpy as np

from .modelbase import Model


class Analysis(Model):
        """Class for model analyses, including common metabolic control
        analysis calculations

        Methods
        -------
        numericElasticities(y0, rate)
            Calculate numerical elasticities for a rate
        allElasticities(y0, norm)
            Numerically approximates elasticisties for all rates
        numericJacobian(y0, **kwargs)
            Calculate Jacobian
        findSteadyState(y0, **kwargs)
            Tries to find the steady-state by numerically solving the algebraic system dy/dt = 0
        concentrationControlCoefficients(y0, pname, norm, **kwargs)
            Calculates concentration control coefficients for a parameter
        """

    def numericElasticities(self, y0, rate):
        """Numerically approximates elasticisties for a rate

        Parameters
        ----------
        y0 : list or numpy.array
            State vector
        rate : str
            Name of the rate for which elasticities are calculated

        Returns
        -------
        epsilon : numdifftools.Jacobian
            elasticities
        """
        def vi(y):
            v = self.rates(y)
            return v[rate]
        jac = nd.Jacobian(vi,step=y0.min()/100)
        epsilon = jac(y0)
        return epsilon

    def allElasticities(self, y0, norm=False):
        """Numerically approximates elasticisties for all rates

        Parameters
        ----------
        y0 : list or numpy.array
            State vector
        norm : bool
            Normalization for elasticities

        Returns
        -------
        epsilon : numpy.matrix
            Matrix of elasticities
        """
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
        """Calculate Jacobian

        Parameters
        ----------
        y0 : list or numpy.array
            State vector for which Jacobian is calculated
        Returns
        -------
        J : numpy.matrix
            Jacobian
        """
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
        """Tries to find the steady-state by numerically solving the algebraic system dy/dt = 0.

        Parameters
        ----------
        y0 : list or numpy.array
            Initial guess

        Returns
        -------
        sol : list
            Possible solution
        """
        def fn(x):
            return self.model(x, 0, **kwargs)
        sol = opt.root(fn, y0)
        if sol.success == True:
            return sol.x
        else:
            return False

    def concentrationControlCoefficients(self, y0, pname, norm=True, **kwargs):
        """Calculates concentration control coefficients for a parameter.
        Uses findSteadyState.

        Parameters
        ----------
        y0 : list or numpy.array
            Initial steady-state guess
        pname : str
            Parameter name to vary
        norm : bool
            Whether to normalize coefficients. Defaults to True.

        Returns
        -------
        cc : list
            Response coefficients
        """
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
