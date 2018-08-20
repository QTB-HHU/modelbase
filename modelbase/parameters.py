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

"""ParameterSet

Description of the module

"""


class ParameterSet(object):
    """Class containing model paramters

    Attributes
    ----------
    pars : dict
        Supplied parameters
    defaultpars : dict
        Default parameters overwriting supplied parameters

    Methods
    -------
    update(pars)
        Adds and updates parameters.
    """


    def __init__(self, parameters={}):
        if isinstance(parameters, dict):
            for k,v in parameters.items():
                setattr(self,k,v)
        elif isinstance(parameters, ParameterSet):
            for k, v in parameters.__dict__.items():
                setattr(self,k,v)
        else:
            raise TypeError("Function requires dict or ParameterSet input")

    def update(self, pars):
        """Adds and updates parameters

        Parameters
        ----------
        pars : dict or modelbase.parameters.ParameterSet
            Object containing new parameters

        Returns
        -------
        None

        Warns
        -----
        OverwritingKeys
            Prints warning if keys are overwritten
        """
        if isinstance(pars,dict):
            replaced_keys = [key for key in self.__dict__.keys() if key in pars]
            if replaced_keys:
                print("Warning: overwriting keys", replaced_keys)
            for k,v in pars.items():
                setattr(self,k,v)
        elif isinstance(pars,ParameterSet):
            self.update(pars.__dict__)
        else:
            raise TypeError("Function requires dict or ParameterSet input")
