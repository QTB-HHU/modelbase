__author__ = 'oliver'

"""
contains some useful rate laws
"""

def massAction(p, *args):

    v = p
    for x in args:
        v = v * x

    return v
