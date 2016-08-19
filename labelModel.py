__author__ = 'oliver'

from modelbase.model import Model
from modelbase.algmModel import AlgmModel
from modelbase.algebraicModule import AlgebraicModule

import itertools

import numpy as np


def generateLabelCpds(cpdName, c):
    '''
    generates label versions of a compound.
    input: string cpdName, int c (number of carbon atoms)
    output: list of compounds names with all labeling patterns accroding to Name000, Name001 etc
    '''

    cpdList = [cpdName+''.join(i) for i in itertools.product(('0','1'), repeat = c)]

    return cpdList

def mapCarbons(sublabels, carbonmap):
    '''
    generates a redistributed string for the substrates (sublabels) according to carbonmap
    '''
    prodlabels = ''.join([sublabels[carbonmap[i]] for i in range(len(carbonmap))])
    return prodlabels

def splitLabel(label, numc):
    '''
    splits the label string according to the lengths given in the list/vector numc
    '''
    splitlabels = []
    cnt = 0
    for i in range(len(numc)):
        splitlabels.append(label[cnt:cnt+numc[i]])
        cnt += numc[i]
    return splitlabels


class LabelModel(AlgmModel):
    '''
    LabelModel allows to define a model with carbon labelling pattern information
    
    Important information on usage:
    -------------------------------
    Compounds must be added with the add_cpd method, which here takes two arguments:
    cpdName (string) and c (int) specifying number of carbon atoms
    '''

    def __init__(self, pars={}, defaultpars={}):
        super(LabelModel,self).__init__(pars,defaultpars)
        self.cpdBaseNames = {}


    def add_cpd(self, cpdName, c):
        '''
        adds compound to model, generating all possible labelling patterns
        '''
        self.cpdBaseNames[cpdName] = c
        labelNames = generateLabelCpds(cpdName,c)
        super(LabelModel,self).add_cpds(labelNames) # add all labelled names

        # now define an algebraic module for the sum of all labels
        # e.g. if CO20, CO21 are the unlabelled and labelled CO2's,
        # the total can be accessed by 'CO2' (likewise for any other more complicated compound)
        def totalconc(par, y):
            return np.array([y.sum()])
        tc = AlgebraicModule({}, totalconc)
        self.add_algebraicModule(tc,labelNames,[cpdName])

    def set_base_rate(self, rateBaseName, fn, numsubs, *args):
        '''
        sets an identical rate expression for all isotope labelling patters of the substrates
        numsubs: int defining the number of substrates in *args
        '''
        
        # first collect the lengths (num of C) of the substrates
        cs = np.array([self.cpdBaseNames[i] for i in args[0:numsubs]])
        otherargs = [] # collect the rest
        if numsubs < len(args):
            otherargs = args[numsubs:len(args)]


        # get all possible combinations of label patterns
        rateLabels = generateLabelCpds('',cs.sum())

        for l in rateLabels:
            rateName = rateBaseName+l
            cnt=0
            sublabels = []
            for i in range(numsubs):
                sublabels.append(args[i]+l[cnt:cnt+cs[i]])
                cnt += cs[i]
                # print i, sublabels[i]

            newargs = sublabels+otherargs
            # print newargs
            self.set_rate(rateName, fn, *newargs)

    def add_carbonmap_reaction(self, rateBaseName, fn, carbonmap, subList, prodList, *args):
        '''
        sets all rates for reactions for all isotope labelling patterns of the substrates.
        Sets all stoichiometries for these reactions.
        requires additionally
        - carbonmap: a list defining how the carbons appear in the products
          (of course, number of Cs must be the same for substrates and products)
        - subList: list of substrates
        - prodList: list of products
        - *args: list of arguments required to calculate rate using function fn 
          (including substrates and possibly allosteric effectors). 
          In this list, substrate names MUST come first

        examples for carbon maps:
        TPI: GAP [0,1,2] -> DHAP [2,1,0] (order changes here), carbonmap = [2,1,0]
        Ald: DHAP [0,1,2] + GAP [3,4,5] -> FBP, carbonmap = [0,1,2,3,4,5]
        TK: E4P [0,1,2,3] + X5P [4,5,6,7,8] -> GAP [6,7,8] + F6P [4,5,0,1,2,3], carbonmap = [6,7,8,4,5,0,1,2,3]
        '''

        # first collect the lengths (num of C) of the substrates and products
        cs = np.array([self.cpdBaseNames[s] for s in subList])
        cp = np.array([self.cpdBaseNames[p] for p in prodList])

        # get all args from *args that are not substrates (can be passed directly)
        otherargs = list(args[len(cs):len(args)])
        print "otherargs:", otherargs

        # get all possible combinations of label patterns for substrates
        rateLabels = generateLabelCpds('',cs.sum())

        for l in rateLabels: # loop through all patterns
            print l
            pl = mapCarbons(l, carbonmap) # get product labels
            sublabels = splitLabel(l, cs)
            prodlabels = splitLabel(pl, cp)

            subargs = [args[i]+sublabels[i] for i in range(len(cs))]
            print subargs
            prodargs = [prodList[i]+prodlabels[i] for i in range(len(cp))]
            print prodargs

            rateName = rateBaseName+l

            # set rate
            rateargs = subargs+otherargs
            print rateargs
            self.set_rate(rateName, fn, *rateargs)

            # set stoichiometry dictionary
            # FIXME think about the possibility that a stoichiometry is not +/-1...
            stDict = {k:-1 for k in subargs}
            stDict.update({k:1 for k in prodargs})
            print stDict
            self.set_stoichiometry(rateName, stDict)
            




if __name__ == '__main__':

    m=LabelModel({'k':1.})
    m.add_cpd('GAP',3)
    m.add_cpd('DHAP',3)
    def fn(par,x):
        return par.k*x
    m.add_carbonmap_reaction('TPI',fn,[2,1,0],['GAP'],['DHAP'],'GAP')
    y=np.array(range(8))*.1
    y2=np.hstack([y,np.zeros(8)])
    print m.rates(y2)
