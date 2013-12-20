''' Bspline module for CADRE '''

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

from MBI.MBI import MBI

import numpy as np


class BsplineParameters(Component):

    '''Creates a Bspline interpolant for several CADRE variables
       so that their time histories can be shaped with m control points
       instead of n time points.'''

    def __init__(self, n, m):
        super(BsplineParameters, self).__init__()

        self.n = n
        self.m = m
        self.add('t1', Float(0.,
                             units='s',
                             desc='Start time',
                             iotype='in'))

        self.add('t2',  Float(43200.,
                              units='s',
                              desc='End time',
                              iotype='in'))

        self.B = MBI(np.zeros(n),
                     [np.linspace(self.t1,self.t2,n)], [self.m], [4]).getJacobian(0,0)

        self.Bdot = MBI(np.zeros(n),
                        [np.linspace(self.t1,self.t2,n)], [self.m], [4]).getJacobian(1,0)

        self.BT = self.B.transpose()
        self.BdotT = self.Bdot.transpose()

        self.add('CP_P_comm', Array(np.zeros((self.m,)),
                                    size=(self.m,),
                                    dtype=float,
                                    units='W',
                                    desc='Communication power at the control points',
                                    iotype='in'))

        self.add('CP_gamma', Array(np.zeros((self.m,)),
                                   size=(self.m,),
                                   dtype=float,
                                   units='rad',
                                   desc='Satellite roll angle at control points',
                                   iotype='in'))

        self.add('CP_Isetpt', Array(np.zeros((12,self.m)),
                                    size=(12,self.m),
                                    dtype=float,
                                    units='A',
                                    desc='Currents of the solar panels at the control points',
                                    iotype='in'))

        self.add('P_comm', Array(np.ones((n,)),
                                 size=(n,), dtype=float,
                                 units='W',
                                 desc='Communication power over time',
                                 iotype='out'))

        self.add('Gamma', Array(0.1*np.ones((n,)),
                                size=(n,),
                                dtype=float,
                                units='rad',
                                desc='Satellite roll angle over time',
                                iotype='out'))

        self.add('Isetpt',Array(0.2*np.ones((12,n)),
                                size=(12,n),
                                dtype=float,
                                units="A",
                                desc="Currents of the solar panels over time",
                                iotype='out'))

    def linearize(self):
        """ Calculate and save derivatives (i.e., Jacobian). """
        # Derivatives are simple
        return

    def execute(self):
        """ Calculate output. """

        self.P_comm = self.B.dot(self.CP_P_comm)
        self.Gamma = self.B.dot(self.CP_gamma)
        for k in range(12):
            self.Isetpt[k, :] = self.B.dot(self.CP_Isetpt[k, :])

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian """

        if 'CP_P_comm' in arg:
            result['P_comm'] += self.B.dot(arg['CP_P_comm'])

        if 'CP_gamma' in arg:
            result['Gamma'] += self.B.dot(arg['CP_gamma'])

        if 'CP_Isetpt' in arg:
            for k in range(12):
                result['Isetpt'][k, :] += self.B.dot(arg['CP_Isetpt'][k, :])

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian """

        if 'P_comm' in arg and 'CP_P_comm' in result:
            result['CP_P_comm'] += self.BT.dot(arg['P_comm'])

        if 'Gamma' in arg and 'CP_gamma' in result:
            result['CP_gamma'] += self.BT.dot(arg['Gamma'])

        if 'Isetpt' in arg and 'CP_Isetpt' in result:
            for k in range(12):
                result['CP_Isetpt'][k, :] += self.BT.dot(arg['Isetpt'][k, :])
