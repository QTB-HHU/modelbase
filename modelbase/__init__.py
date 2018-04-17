from .model import Model
from .model import AlgmModel
from .model import LabelModel

try:
    from .assimulate import Assimulate
    from .assimulate import AlgmAssimulate
    from .assimulate import LabelAssimulate
    Simulate = Assimulate
    AlgmSimulate = AlgmAssimulate
    LabelSimulate = LabelAssimulate
except:
    print("Could not load modelbase.assimulate. Sundials support disabled.")
    from .simulate import Simulate
    from .simulate import AlgmSimulate
    from .simulate import LabelSimulate

