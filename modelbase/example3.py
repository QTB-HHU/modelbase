__author__ = 'oliver'

"""
illustrates the usage of time dependent external conditions 
"""

if __name__ == '__main__':

    import modelbase.model
    import modelbase.simulate

    import matplotlib.pyplot as plt

    import numpy as np

    print("Example 3 started...")


    cl = ['X']
    p = {'l':1,'k':0.1}

    m = modelbase.model.Model(p)

    m.set_cpds(cl)

    def v0(p,**kwargs):
        return np.exp(-p.l*kwargs['t'])

    m.set_ratev('v0',v0)
    m.set_stoichiometry('v0',{'X':1})

    def v1(p,x):
        return p.k*x

    m.set_rate('v1',v1,'X')
    m.set_stoichiometry('v1',{'X':-1})

    s = modelbase.simulate.Simulate(m)
    s.timeCourse(np.linspace(0,50,500),np.zeros(1))

    plt.figure()
    plt.plot(s.getT(),s.getVar([0]))
    plt.draw()

    print("OK!")

