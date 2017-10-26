__author__ = 'oliver'

if __name__ == '__main__':

    import modelbase.model
    import modelbase.simulate
    #import modelbase.results

    import matplotlib.pyplot as plt

    import numpy as np

    print("Example 1 started...")

    cl = ['X','Y']
    p = {'v0':1,'k1p':0.5,'k1m':1,'k2':0.1}

    m = modelbase.model.Model(p)

    m.set_cpds(cl)

    m.set_rate('v0',lambda p:p.v0)
    m.set_stoichiometry('v0',{'X':1})

    def v1(p,x,y):
        return p.k1p*x - p.k1m*y

    m.set_rate('v1',v1,'X','Y')
    m.set_stoichiometry('v1',{'X':-1,'Y':1})

    def v2(p,y):
        return p.k2*y

    m.set_rate('v2',v2,'Y')
    m.set_stoichiometry('v2',{'Y':-1})

    s = modelbase.simulate.Simulate(m)
    s.timeCourse(np.linspace(0,100,1000),np.zeros(3))

    #r = modelbase.results.Results(s)

    #plt.interactive(True)
    plt.figure()
    plt.plot(s.getT(),s.getVar([0,1]))
    plt.draw()
    plt.show()

    print("OK!")

