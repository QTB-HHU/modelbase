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

