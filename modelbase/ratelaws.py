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

"""RateLaws

Description of the module

"""
def massAction(p, *args):
    """Mass action law

    Parameters
    ----------
    p : int or float
        Parameter
    args : int or float
        Concentrations
    Returns
    -------
    v : float
        Rate law
    """
    v = p
    for x in args:
        v = v * x
    return v

def MM1(Vmax, KM, X):
    """Michaelis-Menten rate law

    Parameters
    ----------
    p : int or float
        Parameter
    args : int or float
        Concentrations
    Returns
    -------
    v : float
        Rate law
    """
    return Vmax * X / (KM + X)

def irreversibleMMUni(Vmax,KM):
    """Irreversible Michaelis-Menten rate law

    Parameters
    ----------
    p : int or float
        Parameter
    args : int or float
        Concentrations
    Returns
    -------
    v : float
        Rate law
    """
    def _rateLaw(p,x):
        return getattr(p,Vmax)*x/(getattr(p,KM)+x)
    return _rateLaw

def reversibleMassActionUniUni(kf,eq):
    """xyz law

    Parameters
    ----------
    p : int or float
        Parameter
    args : int or float
        Concentrations
    Returns
    -------
    v : float
        Rate law
    """
    def _rateLaw(p,x,y):
        return getattr(p,kf)*(x-y/getattr(p,eq))
    return _rateLaw

def reversibleMassActionBiUni(kf,eq):
    """xyz law

    Parameters
    ----------
    p : int or float
        Parameter
    args : int or float
        Concentrations
    Returns
    -------
    v : float
        Rate law
    """
    def _rateLaw(p,x,y,z):
        return getattr(p,kf)*(x*y-z/getattr(p,eq))
    return _rateLaw

def reversibleMassActionUniBi(kf,eq):
    """xyz law

    Parameters
    ----------
    p : int or float
        Parameter
    args : int or float
        Concentrations
    Returns
    -------
    v : float
        Rate law
    """
    def _rateLaw(p,x,y,z):
        return getattr(p,kf)*(x-y*z/getattr(p,eq))
    return _rateLaw

def reversibleMassActionBiBi(kf,eq):
    """xyz law

    Parameters
    ----------
    p : int or float
        Parameter
    args : int or float
        Concentrations
    Returns
    -------
    v : float
        Rate law
    """
    def _rateLaw(p,a,b,c,d):
        return getattr(p,kf)*(a*b-c*d/getattr(p,eq))
    return _rateLaw

def irrMMnoncompInh(k_Vmax,k_KM,k_KI):
    """xyz law

    Parameters
    ----------
    p : int or float
        Parameter
    args : int or float
        Concentrations
    Returns
    -------
    v : float
        Rate law
    """
    def _rateLaw(p,X,I):
        Vmax = getattr(p,k_Vmax)
        KM = getattr(p,k_KM)
        KI = getattr(p,k_KI)
        return Vmax*(X/(KM+X))/(1+I/KI)
    return _rateLaw
