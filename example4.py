__author__ = 'oliver'

"""
illustrates the usage of time dependent external conditions.
Simulates a simple kinase cascade.
Exploits algebraicModule for conserved quantitites.
"""

if __name__ == '__main__':

    import modelbase.model
    import modelbase.algmModel
    import modelbase.simulate
    import modelbase.results
    import modelbase.algebraicModule

    import matplotlib.pyplot as plt

    import numpy as np


    cl = ['X','Y','Z']
    p = {'l':.5, 'k1':1., 'k2':1., 'k3':1., 'p':.5}

    m = modelbase.algmModel.AlgmModel(p)

    m.set_cpds(cl)

    def conrel(par, y):
        return np.array([par.tot - y[0]])

    cr = modelbase.algebraicModule.AlgebraicModule({'tot':1}, conrel)
    m.add_algebraicModule(cr,['X'],['Xi'])
    m.add_algebraicModule(cr,['Y'],['Yi'])
    m.add_algebraicModule(cr,['Z'],['Zi'])

    # time-dependent stimulus
    def v0(p,x,**kwargs):
        return x*np.exp(-p.l*kwargs['t'])

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


    # def constitutively active, unspecific phosphatase
    def phosph(p,x):
        return p.p * x

    m.set_rate('p1',phosph,'X')
    m.set_stoichiometry('p1',{'X':-1})

    m.set_rate('p2',phosph,'Y')
    m.set_stoichiometry('p2',{'Y':-1})

    m.set_rate('p3',phosph,'Z')
    m.set_stoichiometry('p3',{'Z':-1})

    s = modelbase.simulate.Simulate(m)
    s.timeCourse(np.linspace(0,50,500),np.zeros(3))

    r = modelbase.results.Results(s)

    plt.figure()
    plt.plot(r.getT(),r.getVar([0,1,2]))
    plt.draw()

    print "OK!"

