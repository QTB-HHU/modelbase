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

                            
def reversibleMassActionUniBi(kf,eq):

    def _rateLaw(p,x,y,z):
        return getattr(p,kf)*(x-y*z/getattr(p,eq))

    return _rateLaw

                            
def reversibleMassActionBiBi(kf,eq):

    def _rateLaw(p,a,b,c,d):
        return getattr(p,kf)*(a*b-c*d/getattr(p,eq))

    return _rateLaw

                            
def irrMMnoncompInh(k_Vmax,k_KM,k_KI):

    def _rateLaw(p,X,I):
        Vmax = getattr(p,k_Vmax)
        KM = getattr(p,k_KM)
        KI = getattr(p,k_KI)
        return Vmax*(X/(KM+X))/(1+I/KI)

    return _rateLaw
