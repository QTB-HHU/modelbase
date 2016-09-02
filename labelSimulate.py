__author__ = 'oliver'

from modelbase.simulate import Simulate

import itertools

import numpy as np

class LabelSimulate(Simulate):
    '''
    subclass of Simulate, including several access methods for labels
    '''

    def getTotal(self, cpdBaseName, r=None):
        '''
        retrieves total concentration for compound cpdBaseName
        :param cpdBaseName: base name of compound
        :return: vector with concentrations
        '''

        if r is None:
            r = range(len(self.results))
        
        c = self.model.cpdBaseNames[cpdBaseName]

        regexp = "\A" + cpdBaseName + '.' * c + "\Z"

        Y = self.getVarsByRegexp(regexp, r)

        return Y.sum(1)

    def getLabelAtPos(self, cpdBaseName, lpos, r=None):
        '''
        retrieves total of cpdBaseName where label is at lpos
        :param cpdBaseName: base name of compound
        :param lpos: position of label
        :return: vector with concentrations
        '''

        if r is None:
            r = range(len(self.results))
        
        if type(lpos) == int:
            lpos = [lpos]

        c = self.model.cpdBaseNames[cpdBaseName]

        l = ['.'] * c
        for p in lpos:
            l[p] = '1'

        regexp = "\A" + cpdBaseName + ''.join(l) + "\Z"

        Y = self.getVarsByRegexp(regexp, r)

        return Y.sum(1)

    def getNumLabel(self, cpdBaseName, nlab, r=None):
        '''
        retrieves total of cpdBaseName with exactly nlab labels
        :param cpdBaseName: base name of compound
        :param nlab: numbers of labels
        :return: vector with concentrations
        '''

        if r is None:
            r = range(len(self.results))

        c = self.model.cpdBaseNames[cpdBaseName]
        
        lcom = itertools.combinations(range(c),nlab)

        cpdNames = []
        for i in lcom:
            l = ['0'] * c
            for p in i:
                l[p] = '1'
            cpdNames.append(cpdBaseName + ''.join(l))

        Y = self.getVarsByName(cpdNames, r)

        return Y.sum(1)


    def getTotalLabel(self, cpdBaseName, r=None):
        '''
        retrieves total labels of cpdBaseName
        :param cpdBaseName: base name of compound
        :return: vector with concentrations
        '''

        if r is None:
            r = range(len(self.results))

        c = self.model.cpdBaseNames[cpdBaseName]
        
        Ylab = []
        for lnum in range(1,c):
            Ylab.append(self.getNumLabel(cpdBaseName, lnum) * lnum)

        return np.vstack(Ylab).sum(0)
