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
    if isinstance(model,LabelModel):
        return LabelSimulate(model)
    elif isinstance(model,Model):
        return Simulate(model)
    else:
        raise NotImplementedError
        