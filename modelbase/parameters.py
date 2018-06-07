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
        Updates parameters
        Input:
            pars: dictionary or ParameterSet object
        '''

        if isinstance(pars,dict):
            replaced_keys = [key for key in self.__dict__.keys() if key in pars]
            if replaced_keys:
                print("Warning: overwriting keys",replaced_keys)
            for k,v in pars.items():
                setattr(self,k,v)
        elif isinstance(pars,ParameterSet):
            self.update(pars.__dict__)
        else:
            raise TypeError("Function requires dict or ParameterSet input")
