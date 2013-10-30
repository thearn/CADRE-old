import numpy as np

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array


class KSfunction(object): 
    """Helper class that can be used inside other components to aggregate constraint 
    vectors with a KS function"""

    def compute(self, g, rho=50): 
        """gets the value of the KS function for the given array of constraints"""

        self.rho = rho
        self.g_max = np.max(g)
        self.g_diff = g-self.g_max
        self.exponents = np.exp(rho * self.g_diff)
        self.summation = np.sum(self.exponents)
        self.KS = self.g_max + 1.0/rho * np.log(self.summation)

        return self.KS

    def derivatives(self):
        """returns a row vector of [dKS_gd, dKS_drho]"""
        dsum_dg = self.rho*self.exponents
        dKS_dsum = 1.0/self.rho/self.summation
        self.dKS_dg = dKS_dsum * dsum_dg

        dsum_drho = np.sum(self.g_diff*self.exponents)
        self.dKS_drho = dKS_dsum * dsum_drho

        return self.dKS_dg, self.dKS_drho



class KSComp(Component): 
    """Aggregates a number of functions to a single value via the 
    Kreisselmeier-Steinhauser Function""" 

    rho = Float(.1, iotype="in", desc="hyperparameter for the KS function")
    KS = Float(0, iotype="out", desc="value of the aggregate KS function")

    def __init__(self, n=2): 
        super(KS, self).__init__()

        self.n = n

        self.add('g',Array(zeros((n,)), size=(n,1), dtype=Float, iotype="in", 
            desc="array of function values to be aggregated"))

        self._ks = KSfunction()

    def execute(self): 
        self.KS = self._ks.compute(self.g, self.rho)

    def linearize(self): 
        """linearize around the last executed point""" 

        #use g_max, exponsnte, summation from last executed point
        self.J = np.hstack(self._ks.derivatives())


    def provideDer(self): 
        ins = ('g','rho')
        outs = ('KS', )
        return ins, outs, self.J

