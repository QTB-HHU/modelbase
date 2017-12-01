__author__ = 'oliver'

"""
This example illustrates the usage of time dependent external conditions 
Simulates is a simple 2-reaction chain -> X ->
Influx is time-dependent with exp(-p*t)
Outflux is mass-action
"""

if __name__ == '__main__':

    import modelbase

    import matplotlib.pyplot as plt

    import numpy as np

    print("Example 3 started...")


    # define variable
    cl = ['X']
    p = {'l':1,'k':0.1}

    # instantiate model
    m = modelbase.Model(p)

    m.set_cpds(cl)

    # define influx
    # if additional kwargs are used in the function, by default the key 't' holds the time
    def v0(p,**kwargs):
        return np.exp(-p.l*kwargs['t'])

    # this rate must be set with 'ratev'
    m.set_ratev('v0',v0)
    m.set_stoichiometry('v0',{'X':1})

    def v1(p,x):
        return p.k*x

    m.set_rate('v1',v1,'X')
    m.set_stoichiometry('v1',{'X':-1})

    s = modelbase.Simulate(m)
    s.timeCourse(np.linspace(0,50,500),np.zeros(1))

    plt.figure()
    plt.plot(s.getT(),s.getVar([0]))
    plt.draw_if_interactive()
    plt.show()

    print("OK!")

