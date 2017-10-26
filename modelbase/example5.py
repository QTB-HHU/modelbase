#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:23:32 2017

@author: oliver
"""

if __name__ == '__main__':

    import modelbase.model as mod
    import modelbase.simulate as sim
    import modelbase.ratelaws as rl
    
    import matplotlib.pyplot as plt

    import numpy as np

    print("Example 5 started...")

    m = mod.LabelModel()
    
    m.par.update({'kf_TPI': 1.0,
                  'Keq_TPI': 21.0,
                  'kf_Ald': 2000.0,
                  'Keq_Ald': 7000.0})
    
    m.add_base_cpd('GAP',3)
    m.add_base_cpd('DHAP',3)
    m.add_base_cpd('FBP',6)

    def v1f(p,y):
        return rl.massAction(p.kf_TPI,y)

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
    y0 = m.set_initconc_cpd_labelpos(y0d,labelpos={'GAP':0})
    
    s = sim.LabelSimulate(m)
    T = np.linspace(0,20,1000)
    s.timeCourse(T,y0)
    
    plt.figure()
    plt.plot(s.getT(),np.vstack([s.getLabelAtPos('FBP',i) for i in range(6)]).transpose())
    plt.legend([str(i+1) for i in range(6)])
    plt.title("Label in FBP")
    plt.draw_if_interactive()
    
    print("OK!")
    