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

"""Analysis

Functions for model analyses, including common metabolic control
analysis calculations

"""
import scipy.optimize as opt
import numdifftools as nd
import numpy as np

class Analysis:
    @staticmethod
    def numericElasticities(model, y0, rate):
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
            v = model.rates(y)
            return v[rate]
        jac = nd.Jacobian(vi, step=y0.min()/100)
        epsilon = jac(y0)
        return epsilon

    @staticmethod
    def allElasticities(model, y0, norm=False):
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
        rateIds = model.rateNames()
        epsilon = np.zeros([len(rateIds), len(model.cpdNames)])
        for i in range(len(rateIds)):
            def vi(y):
                return model.rateFn[rateIds[i]](y)
            jac = nd.Jacobian(vi, step=y0.min()/100)
            epsilon[i,:] = jac(y0)
        if norm:
            v = np.array(model.rates(y0).values())
            epsilon = (1/v).reshape(len(v),1)*epsilon*y0
        return np.matrix(epsilon)

    @staticmethod
    def numericJacobian(model, y0, **kwargs):
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
                dydt = model.model(y, 0, **kwargs)
                return dydt[i]
            jac = nd.Jacobian(fi,step=jstep)
            J[i,:] = jac(y0)
        return np.matrix(J)

    @staticmethod
    def findSteadyState(model, y0, **kwargs):
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
            return model.model(x, 0, **kwargs)
        sol = opt.root(fn, y0)
        if sol.success == True:
            return sol.x
        else:
            return False

    @staticmethod
    def concentrationControlCoefficients(model, y0, pname, norm=True, **kwargs):
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
        origValue = getattr(model.par, pname)
        def fn(x):
            model.par.update({pname: x})
            return findSteadyState(y0, **kwargs)
        jac = nd.Jacobian(fn, step=origValue/100.)
        cc = np.array(jac(origValue))
        model.par.update({pname: origValue})
        if norm:
            ss = findSteadyState(y0, **kwargs)
            cc = origValue * cc / ss.reshape(ss.shape[0],1)
        return cc
