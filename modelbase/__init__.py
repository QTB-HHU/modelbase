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

"""Init

Description of the module

"""
from .model import Model
from .model import LabelModel

try:
    from .assimulate import Assimulate
    from .assimulate import LabelAssimulate
    Simulate = Assimulate
    LabelSimulate = LabelAssimulate
except:
    print("Could not load modelbase.assimulate. Sundials support disabled.")
    from .simulate import Simulate
    from .simulate import LabelSimulate

from .analysis import Analysis


def Simulator(model):
    """ Chooses the simulator class according to the model type

    Parameters
    ----------
    model : modelbase.model
        The model instance

    Returns
    -------
    Simulate : object
        A simulate object according to the model type
    """
    if isinstance(model, LabelModel):
        return LabelSimulate(model)
    elif isinstance(model, Model):
        return Simulate(model)
    else:
        raise NotImplementedError
