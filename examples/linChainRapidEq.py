__author__ = 'oliver'

'''
This example illustrates how to simulate a simple chain of three reactions with modelbase
using a rapid equilibrium approximation.
The reaction chain is -> X <=> Y ->
Here, X <=> Y is assumed to be very fast, i.e. in equilibrium

Approach: 
    1. We define one 'slow' variable A (= X+Y)
    2. We define the algebraic module, allowing to calculate X and Y from A
    3. We can use any quantity A, X, Y in our rate equations
'''

if __name__ == '__main__':

    import modelbase.model
    import modelbase.simulate
    import modelbase.algebraicModule

    import matplotlib.pyplot as plt

    import numpy as np

    print("Example 2 started...")


    # define slow variable A
    cl = ['A']
    p = {'v0':1,'k2':0.1}

    # instantiate model as AlgmModel (because it uses an algebraic module)
    m = modelbase.model.AlgmModel(p)

    m.set_cpds(cl)

    # this function defines the algebraic module. It accepts as first argument
    # the parameters, then the slow-changing variable
    # output are the two derived variables X and Y
    def feq(par,y):
        return np.array([y[0]/(1+par.K),y[0]*par.K/(1+par.K)])

    # define the algebraic module object
    eqm = modelbase.algebraicModule.AlgebraicModule({'K':5},feq)

    # add it to the model by specifying the names of the variables
    m.add_algebraicModule(eqm,['A'],['X','Y'])


    # constant influx to the pool A
    m.set_rate('v0',lambda p:p.v0)
    m.set_stoichiometry('v0',{'A':1})

    # mass-action outflux from the pool A
    # note that rate expression depends on variable Y!
    def v2(p,y):
        return p.k2*y

    m.set_rate('v2',v2,'Y')
    m.set_stoichiometry('v2',{'A':-1})


    # use the AlgmSimulate class to get access to the variables X and Y
    s = modelbase.simulate.AlgmSimulate(m)
    s.timeCourse(np.linspace(0,100,1000),np.zeros(1))

    #a = s.getVar([0])
    #xy = np.array([eqm.getConcentrations(np.array([z])) for z in a])
    
    #plt.figure()
    #plt.plot(s.getT(),a)
    #plt.plot(s.getT(),xy)
    #plt.draw()
    #plt.show()
    
    plt.figure()
    plt.plot(s.getT(),s.getY())
    plt.legend(m.allCpdNames())
    plt.draw_if_interactive()
    plt.show()

    print("OK!")

