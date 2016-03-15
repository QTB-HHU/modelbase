# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 22:19:33 2015

@author: oliver
"""


# TODO: allow dictionary or ParameterSet as pars / defaultpars

class ParameterSet(object):
    
    def __init__(self, pars={}, defaultpars={}):  # -- Anna changed here for pars to be optional
        ''' 
        sets parameters to defaultpars, 
        overwrites these with pars if provided
        '''
        
        mypars = pars.copy()

        for k in defaultpars.keys():
            mypars.setdefault(k,defaultpars[k])

        for k,v in mypars.items():
            setattr(self,k,v)



    def update(self, pars):
        '''
        updates parameters with keys/values provided in pars
        '''

        for k,v in pars.items():
            setattr(self,k,v)

