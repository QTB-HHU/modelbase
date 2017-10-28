#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:23:32 2017

@author: oliver
"""

"""
This example illustrates how to construct an isotope-label specific model.

The first example simulates the equilibration of isotope distribution in the 
aldolase reaction.
Included reactions are triose-phosphate isomerase and fructose-bisphosphate aldolase:
    GAP <=> DHAP
    GAP + DHAP <=> FBP
    
The second example starts with unlabelled intermediates and simulates a 
constant influx of label into the 1/position of GAP, and a mass-action outflux
of FBP (all isotope variants)
"""

if __name__ == '__main__':

    import modelbase.model as mod
    import modelbase.simulate as sim
    import modelbase.ratelaws as rl
    
    import matplotlib.pyplot as plt

    import numpy as np

    print("Example 5 started...")

    # instantiate model
    m = mod.LabelModel()
    
    #define parameters
    m.par.update({'kf_TPI': 1.0,
                  'Keq_TPI': 21.0,
                  'kf_Ald': 2000.0,
                  'Keq_Ald': 7000.0})
    
    # set the 'base' compounds. Second argument defines numbers of carbons
    # these will be automatically expanded into 8 (=2^3), 8, and 64 (=2^6) 
    # isotope variants, respectively
    m.add_base_cpd('GAP',3)
    m.add_base_cpd('DHAP',3)
    m.add_base_cpd('FBP',6)

    # define a simple mass-action rate-law for the forward TPI reaction
    # GAP -> DHAP
    def v1f(p,y):
        return rl.massAction(p.kf_TPI,y)

    # tell the model this is a carbon map reaction
    # arguments:
    # name, rate function, carbon map (here the numbers are reversed!),
    # list of substrates, list of products, variables to be passed
    m.add_carbonmap_reaction('TPIf',v1f,[2,1,0],['GAP'],['DHAP'],'GAP')

    def v1r(p,y):
        return rl.massAction(p.kf_TPI/p.Keq_TPI,y)

    m.add_carbonmap_reaction('TPIr',v1r,[2,1,0],['DHAP'],['GAP'],'DHAP')

    def v2f(p,y,z):
        return rl.massAction(p.kf_Ald,y,z)

    m.add_carbonmap_reaction('Aldf',v2f,[0,1,2,3,4,5],['DHAP','GAP'],['FBP'],'DHAP','GAP')

    def v2r(p,y):
        return rl.massAction(p.kf_Ald/p.Keq_Ald,y)

    m.add_carbonmap_reaction('Aldr',v2r,[0,1,2,3,4,5],['FBP'],['DHAP','GAP'],'FBP')

    # set initial concentrations
    GAP0 = 2.5e-5
    DHAP0 = GAP0 * m.par.Keq_TPI
    FBP0 = GAP0 * DHAP0 * m.par.Keq_Ald
    y0d = {'GAP': GAP0,
           'DHAP': DHAP0,
           'FBP': FBP0}
    
    # simulate equilibration of the labels
    y0 = m.set_initconc_cpd_labelpos(y0d,labelpos={'GAP':0})
    s = sim.LabelSimulate(m)
    T = np.linspace(0,20,1000)
    s.timeCourse(T,y0)
    
    plt.figure()
    plt.plot(s.getT(),np.vstack([s.getLabelAtPos('FBP',i) for i in range(6)]).transpose())
    plt.legend([str(i+1) for i in range(6)])
    plt.title("Label in FBP positions - equilibration")
    plt.show()
    plt.draw_if_interactive()
    
    # now simulate steady influx of label into GAP, position 0
    kout = 0.05
    vin = 2*9e-5*kout
    m.par.update({'vin':vin,'kout':kout})
    
    # influx reaction is simulated as a 'normal' reaction 
    # generating a GAP with label in 1-position (GAP100)
    m.set_rate('vin',lambda p:p.vin)
    m.set_stoichiometry('vin',{'GAP100':1})
    
    def vout(p,y):
        return rl.massAction(p.kout,y)
    m.add_carbonmap_reaction('vout',vout,[0,1,2,3,4,5],['FBP'],[],'FBP')
    
    y0 = m.set_initconc_cpd_labelpos(y0d)
    T = np.linspace(0,100,1000)
    s2 = sim.LabelSimulate(m)
    s2.timeCourse(T,y0)

    plt.plot(s2.getT(),np.vstack([s2.getLabelAtPos('FBP',i)/s2.getTotal('FBP') for i in range(6)]).transpose())
    plt.legend([str(i+1) for i in range(6)])
    plt.title("Relative label in FBP - dynamic influx of label, steady state")
    plt.draw_if_interactive()
       
    print("OK!")
    