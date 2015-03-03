# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 22:19:33 2015

@author: oliver
"""

class ParameterSet:
    
    def __init__(self, pars={}):  # -- Anna changed here for pars to be optional

        for k,v in pars.items():
            setattr(self,k,v)

