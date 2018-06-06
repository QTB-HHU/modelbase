__author__ = 'oliver'

"""
This example simulates a simple protein kinase signalling cascade
It employs time-dependent external stimulus (as exp(l*t))
It illustrates how to use algebraic modules for conserved quantitites

Here: Three protein kinases X, Y, Z. 
    R(t)
     v
  Xi -> X
        v
     Yi -> Y
           v
        Zi -> Z
        
Reversion X/Y/Z -> Xi/Yi/Zi as mass-action
{XYZ}i is inactive form, {XYZ} active form

Three conserved quantitites:
    X + Xi = Xtot
    Y + Yi = Ytot
    Z + Zi = Ztot

Algebraic modules calculate {XYZ}i from {XYZ}
"""

if __name__ == '__main__':

    import modelbase

    import matplotlib.pyplot as plt

    import numpy as np

    print("Example 4 started...")


    # define dynamic variables
    cl = ['X','Y','Z']
    p = {'l':.5, 'k1':1., 'k2':1., 'k3':1., 'p':.5, 'tot':1}

    # instantiate model
    m = modelbase.Model(p)

    m.set_cpds(cl)

    # define the algebraic module exploiting conserved quantitites
    def conrel(par, y):
        return np.array([par.tot - y[0]])

    # add three algebraic modules, always returning the inactive form Xi, etc.
    m.add_algebraicModule(conrel, 'rapidEqx',['X'],['Xi'])
    m.add_algebraicModule(conrel, 'rapidEqy',['Y'],['Yi'])
    m.add_algebraicModule(conrel, 'rapidEqz',['Z'],['Zi'])

    # time-dependent stimulus (see timeDepExt.py)
    def v0(p,x,**kwargs):
        return x*np.exp(-p.l*kwargs['t'])

    # note that rate of conversion depends on Xi (made accessible by algebraic module)
    m.set_ratev('v0',v0,'Xi')
    m.set_stoichiometry('v0',{'X':1})

    # kinases
    def k1(p,x,y0):
        return p.k1*x*y0

    m.set_rate('k1',k1,'X','Yi')
    m.set_stoichiometry('k1',{'Y':1})

    def k2(p,x,y0):
        return p.k2*x*y0

    m.set_rate('k2',k2,'Y','Zi')
    m.set_stoichiometry('k2',{'Z':1})


    # define the constitutively active, unspecific phosphatase
    def phosph(p,x):
        return p.p * x

    m.set_rate('p1',phosph,'X')
    m.set_stoichiometry('p1',{'X':-1})

    m.set_rate('p2',phosph,'Y')
    m.set_stoichiometry('p2',{'Y':-1})

    m.set_rate('p3',phosph,'Z')
    m.set_stoichiometry('p3',{'Z':-1})

    s = modelbase.Simulate(m)
    s.timeCourse(np.linspace(0,50,500),np.zeros(3))

    plt.figure()
    plt.plot(s.getT(),s.getVarsByName(cl))
    plt.legend(cl)
    plt.draw_if_interactive()
    plt.show()

    print("OK!")

