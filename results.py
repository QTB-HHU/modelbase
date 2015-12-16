__author__ = 'oliver'

import numpy as np


class Results():

    def __init__(self, s):
        """
        class of methods inheriting from Simulation useful for plotting and analysis
        """
        self.model = s.model
        self.results = s.results


    def getT(self, r=None):
        """
        :param r: range of steps of simulation for which results we are interested in
        :return: time of all results as one vector
        """
        if r is None:
            r = range(len(self.results))

        T = np.hstack([self.results[i]['t'] for i in r])

        return T

    def getVar(self, j, r=None):
        """
        :param j: int of the state variable [0:PQred, ..., 5:ATPsynth]
        :param r: range of steps of simulation for which results we are interested in
        :return: concentrations of variable j as one vector
        """

        if type(j) == int:
            j = [j]

        if r is None:
            r = range(len(self.results))

        X = np.vstack([self.results[i]['y'][:,j] for i in r])
        #X = np.vstack([np.reshape(self.results[i]['y'][:, j],np.size(self.results[i]['y'],0),len(j)) for i in r])

        if np.size(X,1) == 1:
            X = np.reshape(X,[np.size(X,0)])

        return X

    def getRate(self, rate, r=None):
        """
        :param rate: name of the rate
        :param r: range of steps of simulation for which results we are interested in
        :return: rate with name 'rate' of all results as one vector
        """
        if r is None:
            r = range(len(self.results))

        V = np.array([])

        for i in r:
            t = self.results[i]['t']
            y = self.results[i]['y']

            V = np.hstack([V,
                           np.array(
                           [self.model.rates(y[j])[rate] for j in range(len(t))]
                           )])

        return V
