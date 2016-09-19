__author__ = 'oliver'

"""
contains some useful rate laws
"""

def massAction(p, *args):

    v = p
    for x in args:
        v = v * x

    return v


def MM1(Vmax, KM, X):
    """ returns Michaelis-Menten rate for irreversible reaction with one substrate """
    return Vmax * X / (KM + X)


def irreversibleMMUni(Vmax,KM):

    def _rateLaw(p,x):
        return getattr(p,Vmax)*x/(getattr(p,KM)+x)

    return _rateLaw


def reversibleMassActionUniUni(kf,eq):

    def _rateLaw(p,x,y):
        return getattr(p,kf)*(x-y/getattr(p,eq))

    return _rateLaw


def reversibleMassActionBiUni(kf,eq):

    def _rateLaw(p,x,y,z):
        return getattr(p,kf)*(x*y-z/getattr(p,eq))

    return _rateLaw

                            
