__author__ = 'oliver'

if __name__ == '__main__':

    import modelbase.model
    import modelbase.algmModel
    import modelbase.simulate
    import modelbase.results
    import modelbase.algebraicModule

    import matplotlib.pyplot as plt

    import numpy as np


    cl = ['A']
    p = {'v0':1,'k2':0.1}

    m = modelbase.algmModel.AlgmModel(p)

    m.set_cpds(cl)

    def feq(par,y):
        return np.array([y[0]/(1+par.K),y[0]*par.K/(1+par.K)])

    eqm = modelbase.algebraicModule.AlgebraicModule({'K':5},feq)

    m.add_algebraicModule(eqm,['A'],['X','Y'])


    m.set_rate('v0',lambda p:p.v0)
    m.set_stoichiometry('v0',{'A':1})

    def v2(p,y):
        return p.k2*y

    m.set_rate('v2',v2,'Y')
    m.set_stoichiometry('v2',{'A':-1})


    s = modelbase.simulate.Simulate(m)
    s.timeCourse(np.linspace(0,100,1000),np.zeros(1))

    r = modelbase.results.Results(s)

    a = r.getVar([0])
    xy = np.array([eqm.getConcentrations(np.array([z])) for z in a])
    
    plt.figure()
    plt.plot(r.getT(),a)
    plt.plot(r.getT(),xy)
    plt.draw()

    print "OK!"

