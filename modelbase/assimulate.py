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

"""Assimulate

Description of the module

"""
import numpy as np

from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

from .simulate import Simulate
from .simulate import LabelSimulate


class Assimulate(Simulate):
    """Provides simulation routines

    Attributes
    ----------
    model : class
        modelbase.Model or modelbase.LabelModel instance
    **kwargs : dict
        User defined dynamic changes for rate functions

    Methods
    -------
    generate_integrator(y0, name, verbosity)
        Initializes CVODE integrator
    set_initial_value(y0, t0)
        Sets initial values
    integrate(t, integrator, minstep, maxstep, nsteps)
        Integrates model. Returns values at time point t
    timeCourse(Torig, y0, integrator, minstep, maxstep, nsteps)
        Integration over time vector
    """

    def __init__(self, model, **kwargs):
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
        """Initializes CVODE integrator

        Parameters
        ----------
        y0 : list or numpy.array
            Initial values
        name : str
            Name of the assimulo problem
        verbosity : int
            Verbosity of the integrator. Possible values = [10, 20, 30, 40, 50]
        Returns
        -------
        None
        """
        if y0 is None:
            y0 = np.zeros(len(self.model.cpdNames))
        self.problem = Explicit_Problem(self.f, y0=y0, name=name)
        self.integrator = CVode(self.problem)
        self.integrator.atol = 1e-8
        self.integrator.rtol = 1e-8
        self.integrator.verbosity = verbosity


    def set_initial_value(self, y0, t0=0):
        """Sets initial values

        Parameters
        ----------
        y0 : list or numpy.array
            Initial values
        t0 : int or float
            Initial time. Defaults to 0.

        Returns
        -------
        None
        """
        self.integrator.y = y0
        self.integrator.t = t0


    def integrate(self, t, integrator=None, minstep=None, maxstep=None, nsteps=None):
        """Integrates model. Returns values at time point t.

        .. deprecated:: 0.23
            `minstep` will be removed in version 1.0
            `maxstep` will be removed in version 1.0
            `nsteps` will be removed in version 1.0

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
        """Integration over time vector.

        .. deprecated:: 0.23
            `minstep` will be removed in version 1.0
            `maxstep` will be removed in version 1.0
            `nsteps` will be removed in version 1.0

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
        try:
            t,Y = self.integrator.simulate(T[-1],ncp_list=T)
        except:
            print("Error in timeCourse while integrating with CVode")
            self._successful = False
        if self.doesMonitor() and self.successful():
            self.results.append({'t': T, 'y': Y})
        return np.vstack(Y)

class LabelAssimulate(Assimulate, LabelSimulate):
    pass
