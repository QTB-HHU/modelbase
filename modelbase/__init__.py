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
