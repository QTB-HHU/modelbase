__author__ = 'oliver'

'''
This example illustrates how to simulate a simple chain of three reactions with modelbase
Reactions:
    v0: influx ( -> X), constant rate
    v1: conversion (X -> Y), reversible mass-action
    v2: outflux (Y -> ), irreversible mass-action
'''

if __name__ == '__main__':

    import modelbase.model
    import modelbase.simulate
    #import modelbase.results

    import matplotlib.pyplot as plt

    import numpy as np

    print("Example 1 started...")

    # define metabolite species X and Y
    cl = ['X','Y']
    # define parameters
    p = {'v0':1,'k1p':0.5,'k1m':1,'k2':0.1}

    # instantiate model
    m = modelbase.model.Model(p)

    # tell it which metabolites it has
    m.set_cpds(cl)

    # define the reactions. Always define a rate function, then define the stoichiometries
    # remember: the rate functions always accept the model parameters as first argument,
    # the remaining arguments are metabolite concentrations as defined in the set_rate command
    
    # v0 is a constant function, creates one X
    m.set_rate('v0',lambda p:p.v0)
    m.set_stoichiometry('v0',{'X':1})

    # v1 is a reversible mass-action 
    def v1(p,x,y):
        return p.k1p*x - p.k1m*y

    m.set_rate('v1',v1,'X','Y') # tells the model to pass concentrations of X and Y
    m.set_stoichiometry('v1',{'X':-1,'Y':1}) # one X destroyed, one Y created

    # v2 is irreversible mass-action
    def v2(p,y):
        return p.k2*y

    m.set_rate('v2',v2,'Y')
    m.set_stoichiometry('v2',{'Y':-1})

    # do the simulation: create Simulate object, perform simulations
    
    # create modelbase.simulate.Simulate object
    s = modelbase.simulate.Simulate(m)
    # run a timecourse
    s.timeCourse(np.linspace(0,100,1000),np.zeros(2))

    #r = modelbase.results.Results(s)

    #plt.interactive(True)
    plt.figure()
    plt.plot(s.getT(),s.getY())
    plt.legend(m.cpdNames)
    plt.draw()
    plt.show()

    print("OK!")

