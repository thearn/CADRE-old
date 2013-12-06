''' Communications Discpline for CADRE '''

import numpy as np
import scipy.sparse
import MBI
import os
from openmdao.lib.datatypes.api import Float, Array
from openmdao.main.api import Component

from CADRE.kinematics import fixangles, computepositionspherical, \
     computepositionsphericaljacobian, computepositionrotd,\
     computepositionrotdjacobian

import rk4
7
# Allow non-standard variable names for scientific calc
# pylint: disable-msg=C0103


class Comm_DataDownloaded(rk4.RK4):

    """ Integrate the incoming data rate to compute the time history of data
    downloaded from the satelite."""

    def __init__(self, n_times):
        super(Comm_DataDownloaded, self).__init__()

        # Inputs
        self.add(
            'Dr', 
            Array(
                np.zeros(n_times),
                iotype='in',
                shape=(n_times,),
                units="Gibyte/s",
                desc="Download rate over time"
            )
        )

        # Initial State
        self.add(
            'Data0',
            Array(
                [0.0],
                iotype='in',
                shape=(1,),
                units="Gibyte",
                desc="Initial downloaded data state"
            )
        )

        # States
        self.add(
            'Data',
            Array(
                np.zeros((1, n_times)),
                iotype='out',
                shape=(1, n_times),
                units="Gibyte",
                desc="Downloaded data state over time"
            )
        )

        self.state_var = "Data"
        self.init_state_var = "Data0"
        self.external_vars = ["Dr"]

        self.dfdy = np.array([[0.]])
        self.dfdx = np.array([[1.]])

    def f_dot(self, external, state):
        return external[0]

    def df_dy(self, external, state):
        return self.dfdy

    def df_dx(self, external, state):
        return self.dfdx


class Comm_AntRotation(Component):

    ''' Fixed antenna angle to time history of the quaternion '''

    # Inputs
    antAngle = Float(0., iotype="in", copy=None)

    def __init__(self, n):
        super(Comm_AntRotation, self).__init__()

        # Outputs
        self.add(
            'q_A',
            Array(
                np.zeros((4, n)),
                iotype='out',
                shape=(4, n),
                units="unitless",
                desc="Quarternion matrix in antenna angle frame over time"
            )
        )

        self.dq_dt = np.zeros(4)

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        rt2 = np.sqrt(2)
        self.dq_dt[0] = - np.sin(self.antAngle / 2.) / 2.
        self.dq_dt[1] = np.cos(self.antAngle / 2.) / rt2 / 2.
        self.dq_dt[2] = - np.cos(self.antAngle / 2.) / rt2 / 2.
        self.dq_dt[3] = 0.0

    def execute(self):
        """ Calculate output. """

        rt2 = np.sqrt(2)
        self.q_A[0,:] = np.cos(self.antAngle/2.)
        self.q_A[1,:] = np.sin(self.antAngle/2.) / rt2
        self.q_A[2,:] = - np.sin(self.antAngle/2.) / rt2
        self.q_A[3,:] = 0.0

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'antAngle' in arg and 'q_A' in result:
            for k in xrange(4):
                result['q_A'][k,:] += self.dq_dt[k] * arg['antAngle']

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'q_A' in arg and 'antAngle' in result:
            for k in xrange(4):
                result['antAngle'] += self.dq_dt[k] * np.sum(arg['q_A'][k,:])


class Comm_AntRotationMtx(Component):

    """ Translate antenna angle into the body frame. """

    def __init__(self, n):
        super(Comm_AntRotationMtx, self).__init__()

        self.n = n

        # Inputs
        self.add(
            'q_A',
            Array(
                np.zeros((4, self.n)),
                iotype='in',
                shape=(4, self.n),
                desc="Quarternion matrix in antenna angle frame over time"
            )
        )

        # Outputs
        self.add(
            'O_AB', 
            Array(
                np.zeros((3, 3, self.n)),
                iotype='out',
                shape=(3, 3, self.n),
                units="unitless",
                desc="Rotation matrix from antenna angle to body-fixed frame over time"
            )
        )

        self.J = np.empty((self.n, 3, 3, 4))

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        A = np.zeros((4, 3))
        B = np.zeros((4, 3))
        dA_dq = np.zeros((4, 3, 4))
        dB_dq = np.zeros((4, 3, 4))

        dA_dq[0,:, 0] = (1, 0, 0)
        dA_dq[1,:, 0] = (0, 1, 0)
        dA_dq[2,:, 0] = (0, 0, 1)
        dA_dq[3,:, 0] = (0, 0, 0)

        dA_dq[0,:, 1] = (0, 0, 0)
        dA_dq[1,:, 1] = (0, 0, -1)
        dA_dq[2,:, 1] = (0, 1, 0)
        dA_dq[3,:, 1] = (1, 0, 0)

        dA_dq[0,:, 2] = (0, 0, 1)
        dA_dq[1,:, 2] = (0, 0, 0)
        dA_dq[2,:, 2] = (-1, 0, 0)
        dA_dq[3,:, 2] = (0, 1, 0)

        dA_dq[0,:, 3] = (0, -1, 0)
        dA_dq[1,:, 3] = (1, 0, 0)
        dA_dq[2,:, 3] = (0, 0, 0)
        dA_dq[3,:, 3] = (0, 0, 1)


        dB_dq[0,:, 0] = (1, 0, 0)
        dB_dq[1,:, 0] = (0, 1, 0)
        dB_dq[2,:, 0] = (0, 0, 1)
        dB_dq[3,:, 0] = (0, 0, 0)

        dB_dq[0,:, 1] = (0, 0, 0)
        dB_dq[1,:, 1] = (0, 0, 1)
        dB_dq[2,:, 1] = (0, -1, 0)
        dB_dq[3,:, 1] = (1, 0, 0)

        dB_dq[0,:, 2] = (0, 0, -1)
        dB_dq[1,:, 2] = (0, 0, 0)
        dB_dq[2,:, 2] = (1, 0, 0)
        dB_dq[3,:, 2] = (0, 1, 0)

        dB_dq[0,:, 3] = (0, 1, 0)
        dB_dq[1,:, 3] = (-1, 0, 0)
        dB_dq[2,:, 3] = (0, 0, 0)
        dB_dq[3,:, 3] = (0, 0, 1)

        for i in range(0, self.n):
            A[0,:] = ( self.q_A[0, i], -self.q_A[3, i],  self.q_A[2, i])
            A[1,:] = ( self.q_A[3, i],  self.q_A[0, i], -self.q_A[1, i])
            A[2,:] = (-self.q_A[2, i],  self.q_A[1, i],  self.q_A[0, i])
            A[3,:] = ( self.q_A[1, i],  self.q_A[2, i],  self.q_A[3, i])

            B[0,:] = ( self.q_A[0, i],  self.q_A[3, i], -self.q_A[2, i])
            B[1,:] = (-self.q_A[3, i],  self.q_A[0, i],  self.q_A[1, i])
            B[2,:] = ( self.q_A[2, i], -self.q_A[1, i],  self.q_A[0, i])
            B[3,:] = ( self.q_A[1, i],  self.q_A[2, i],  self.q_A[3, i])

            for k in range(0, 4):
                self.J[i,:,:, k] = np.dot(dA_dq[:,:, k].T, B) + \
                    np.dot(A.T, dB_dq[:,:, k])

    def execute(self):
        """ Calculate output. """

        A = np.zeros((4, 3))
        B = np.zeros((4, 3))

        for i in range(0, self.n):
            A[0,:] = ( self.q_A[0, i], -self.q_A[3, i],  self.q_A[2, i])
            A[1,:] = ( self.q_A[3, i],  self.q_A[0, i], -self.q_A[1, i])
            A[2,:] = (-self.q_A[2, i],  self.q_A[1, i],  self.q_A[0, i])
            A[3,:] = ( self.q_A[1, i],  self.q_A[2, i],  self.q_A[3, i])

            B[0,:] = ( self.q_A[0, i],  self.q_A[3, i], -self.q_A[2, i])
            B[1,:] = (-self.q_A[3, i],  self.q_A[0, i],  self.q_A[1, i])
            B[2,:] = ( self.q_A[2, i], -self.q_A[1, i],  self.q_A[0, i])
            B[3,:] = ( self.q_A[1, i],  self.q_A[2, i],  self.q_A[3, i])

            self.O_AB[:,:, i] = np.dot(A.T, B)

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'q_A' in arg and 'O_AB' in result:
            for u in xrange(3):
                for v in xrange(3):
                    for k in xrange(4):
                        result['O_AB'][u, v,:] += \
                            self.J[:, u, v, k] * arg['q_A'][k,:]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'O_AB' in arg and 'q_A' in result:
            for u in range(3):
                for v in range(3):
                    for k in range(4):
                        result['q_A'][k,:] += self.J[:, u, v, k] * \
                            arg['O_AB'][u, v,:]


class Comm_BitRate(Component):

    ''' Compute the data rate the satellite receives. '''

    # constants
    pi = 2 * np.arccos(0.)
    c = 299792458
    Gr = 10 ** (12.9 / 10.)
    Ll = 10 ** (-2.0 / 10.)
    f = 437e6
    k = 1.3806503e-23
    SNR = 10 ** (5.0 / 10.)
    T = 500.
    alpha = c ** 2 * Gr * Ll / 16.0 / pi ** 2 / f ** 2 / k / SNR / T / 1e6

    def __init__(self, n):
        super(Comm_BitRate, self).__init__()
        self.n = n

        # Inputs
        self.add(
            'P_comm',
            Array(
                np.zeros(self.n),
                iotype='in',
                shape=(self.n, ),
                units="W",
                desc="Communication power over time"
            )
        )

        self.add(
            'gain',
            Array(
                np.zeros(self.n), 
                iotype='in',
                shape=(self.n, ),
                units="unitless",
                desc="Transmitter gain over time"
            )
        )

        self.add(
            'GSdist',
            Array(
                np.zeros(self.n),
                iotype='in',
                shape=(self.n, ),
                units="km",
                desc="Distance from ground station to satellite over time"
            )
        )

        self.add(
            'CommLOS',
            Array(
                np.zeros(self.n),
                iotype='in',
                shape=(self.n, ),
                units="unitless",
                desc="Satellite to ground station line of sight over time"
            )
        )

        # Outputs
        self.add(
            'Dr',
            Array(
                np.zeros(self.n),
                iotype='out',
                shape=(self.n, ),
                units="Gibyte/s",
                desc="Download rate over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        S2 = 0.
        self.dD_dP = np.zeros(self.n)
        self.dD_dGt = np.zeros(self.n)
        self.dD_dS = np.zeros(self.n)
        self.dD_dLOS = np.zeros(self.n)

        for i in range(0, self.n):

            if np.abs(self.GSdist[i]) > 1e-10:
                S2 = self.GSdist[i] * 1e3
            else:
                S2 = 1e-10

            self.dD_dP[i] = self.alpha * self.gain[i] * \
                self.CommLOS[i] / S2 ** 2
            self.dD_dGt[i] = self.alpha * self.P_comm[i] * \
                self.CommLOS[i] / S2 ** 2
            self.dD_dS[i] = -2.0 * 1e3 * self.alpha * self.P_comm[i] * \
                self.gain[i] * self.CommLOS[i] / S2 ** 3
            self.dD_dLOS[i] = self.alpha * \
                self.P_comm[i] * self.gain[i] / S2 ** 2

    def execute(self):
        """ Calculate output. """

        for i in range(0, self.n):
            if np.abs(self.GSdist[i]) > 1e-10:
                S2 = self.GSdist[i] * 1e3
            else:
                S2 = 1e-10
            self.Dr[i] = self.alpha * self.P_comm[i] * self.gain[i] * \
                self.CommLOS[i] / S2 ** 2

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'Dr' in result:
            if 'P_comm' in arg:
                result['Dr'] += self.dD_dP * arg['P_comm']
            if 'gain' in arg:
                result['Dr'] += self.dD_dGt * arg['gain']
            if 'GSdist' in arg:
                result['Dr'] += self.dD_dS * arg['GSdist']
            if 'CommLOS' in arg:
                result['Dr'] += self.dD_dLOS * arg['CommLOS']

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'Dr' in arg:
            if 'P_comm' in result:
                result['P_comm'] += self.dD_dP.T * arg['Dr']
            if 'gain' in result:
                result['gain'] += self.dD_dGt.T * arg['Dr']
            if 'GSdist' in result:
                result['GSdist'] += self.dD_dS.T * arg['Dr']
            if 'CommLOS' in result:
                result['CommLOS'] += self.dD_dLOS.T * arg['Dr']


class Comm_Distance(Component):

    '''Calculates distance from ground station to satellitle'''

    def __init__(self, n):
        super(Comm_Distance, self).__init__()

        self.n = n

        # Inputs
        self.add(
            'r_b2g_A', 
            Array(
                np.zeros((3, self.n)),
                iotype='in',
                shape=(3, self.n),
                units="km",
                desc="Position vector from satellite to ground station in antenna angle frame over time"
            )
        )

        # Outputs
        self.add(
            'GSdist',
            Array(
                np.zeros(self.n),
                iotype='out',
                shape=(self.n,),
                units="km",
                desc="Distance from ground station to satellite over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.J = np.zeros((self.n, 3))

        for i in range(0, self.n):
            norm = np.dot(self.r_b2g_A[:, i], self.r_b2g_A[:, i]) ** 0.5
            if norm > 1e-10:
                self.J[i,:] = self.r_b2g_A[:, i] / norm
            else:
                self.J[i,:] = 0.

    def execute(self):
        """ Calculate output. """

        for i in range(0, self.n):
            self.GSdist[i] = np.dot(
                self.r_b2g_A[:, i], self.r_b2g_A[:, i]) ** 0.5

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'r_b2g_A' in arg and 'GSdist' in result:
            for k in xrange(3):
                result['GSdist'] += self.J[:, k] * arg['r_b2g_A'][k,:]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'GSdist' in arg and 'r_b2g_A' in result:
            for k in xrange(3):
                result['r_b2g_A'][k,:] += self.J[:, k] * arg['GSdist']


class Comm_EarthsSpin(Component):

    ''' Returns the Earth quaternion as a function of time. '''

    def __init__(self, n):
        super(Comm_EarthsSpin, self).__init__()
        self.n = n

        # Inputs
        self.add('t', Array(np.zeros(self.n),
                            iotype='in',
                            shape=(self.n, ),
                            units="s",
                            desc="Time"))

        # Outputs
        self.add('q_E', Array(
                np.zeros((4, self.n)),
                iotype='out',
                shape=(4, self.n),
                units="unitless",
                desc="Quarternion matrix in Earth-fixed frame over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        ntime = self.n
        self.dq_dt = np.zeros((ntime, 4))

        fact = np.pi / 3600.0 / 24.0
        theta = fact * self.t

        self.dq_dt[:, 0] = -np.sin(theta) * fact
        self.dq_dt[:, 3] = -np.cos(theta) * fact

    def execute(self):
        """ Calculate output. """

        fact = np.pi / 3600.0 / 24.0
        theta = fact * self.t

        self.q_E[0,:] = np.cos(theta)
        self.q_E[3,:] = -np.sin(theta)

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 't' in arg and 'q_E' in result:
            for k in range(4):
                result['q_E'][k,:] += self.dq_dt[:, k] * arg['t']

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'q_E' in arg and 't' in result:
            for k in range(4):
                result['t'] += self.dq_dt[:, k] * arg['q_E'][k,:]


class Comm_EarthsSpinMtx(Component):

    ''' Quaternion to rotation matrix for the earth spin. '''

    def __init__(self, n):
        super(Comm_EarthsSpinMtx, self).__init__()

        self.n = n

        # Inputs
        self.add(
            'q_E',
            Array(
                np.zeros((4, self.n)),
                iotype='in',
                shape=(4, self.n),
                units="unitless",
                desc="Quarternion matrix in Earth-fixed frame over time" 
            )
        )

        # Outputs
        self.add(
            'O_IE',
            Array(
                np.zeros((3, 3, self.n)),
                iotype='out',
                shape=(3, 3, self.n),
                units="unitless",
                desc="Rotation matrix from Earth-centered inertial frame to Earth-fixed frame over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        A = np.zeros((4, 3))
        B = np.zeros((4, 3))

        dA_dq = np.zeros((4, 3, 4))
        dB_dq = np.zeros((4, 3, 4))

        self.J = np.zeros((self.n, 3, 3, 4))

        dA_dq[0,:, 0] = (1, 0, 0)
        dA_dq[1,:, 0] = (0, 1, 0)
        dA_dq[2,:, 0] = (0, 0, 1)
        dA_dq[3,:, 0] = (0, 0, 0)

        dA_dq[0,:, 1] = (0, 0, 0)
        dA_dq[1,:, 1] = (0, 0, -1)
        dA_dq[2,:, 1] = (0, 1, 0)
        dA_dq[3,:, 1] = (1, 0, 0)

        dA_dq[0,:, 2] = (0, 0, 1)
        dA_dq[1,:, 2] = (0, 0, 0)
        dA_dq[2,:, 2] = (-1, 0, 0)
        dA_dq[3,:, 2] = (0, 1, 0)

        dA_dq[0,:, 3] = (0, -1, 0)
        dA_dq[1,:, 3] = (1, 0, 0)
        dA_dq[2,:, 3] = (0, 0, 0)
        dA_dq[3,:, 3] = (0, 0, 1)


        dB_dq[0,:, 0] = (1, 0, 0)
        dB_dq[1,:, 0] = (0, 1, 0)
        dB_dq[2,:, 0] = (0, 0, 1)
        dB_dq[3,:, 0] = (0, 0, 0)

        dB_dq[0,:, 1] = (0, 0, 0)
        dB_dq[1,:, 1] = (0, 0, 1)
        dB_dq[2,:, 1] = (0, -1, 0)
        dB_dq[3,:, 1] = (1, 0, 0)

        dB_dq[0,:, 2] = (0, 0, -1)
        dB_dq[1,:, 2] = (0, 0, 0)
        dB_dq[2,:, 2] = (1, 0, 0)
        dB_dq[3,:, 2] = (0, 1, 0)

        dB_dq[0,:, 3] = (0, 1, 0)
        dB_dq[1,:, 3] = (-1, 0, 0)
        dB_dq[2,:, 3] = (0, 0, 0)
        dB_dq[3,:, 3] = (0, 0, 1)

        for i in range(0, self.n):
            A[0,:] = ( self.q_E[0, i], -self.q_E[3, i],  self.q_E[2, i])
            A[1,:] = ( self.q_E[3, i],  self.q_E[0, i], -self.q_E[1, i])
            A[2,:] = (-self.q_E[2, i],  self.q_E[1, i],  self.q_E[0, i])
            A[3,:] = ( self.q_E[1, i],  self.q_E[2, i],  self.q_E[3, i])

            B[0,:] = ( self.q_E[0, i],  self.q_E[3, i], -self.q_E[2, i])
            B[1,:] = (-self.q_E[3, i],  self.q_E[0, i],  self.q_E[1, i])
            B[2,:] = ( self.q_E[2, i], -self.q_E[1, i],  self.q_E[0, i])
            B[3,:] = ( self.q_E[1, i],  self.q_E[2, i],  self.q_E[3, i])

            for k in range(0, 4):
                self.J[i,:,:, k] = np.dot(dA_dq[:,:, k].T, B) + \
                    np.dot(A.T, dB_dq[:,:, k])

    def execute(self):
        """ Calculate output. """

        A = np.zeros((4, 3))
        B = np.zeros((4, 3))

        for i in range(0, self.n):
            A[0,:] = ( self.q_E[0, i], -self.q_E[3, i],  self.q_E[2, i])
            A[1,:] = ( self.q_E[3, i],  self.q_E[0, i], -self.q_E[1, i])
            A[2,:] = (-self.q_E[2, i],  self.q_E[1, i],  self.q_E[0, i])
            A[3,:] = ( self.q_E[1, i],  self.q_E[2, i],  self.q_E[3, i])

            B[0,:] = ( self.q_E[0, i],  self.q_E[3, i], -self.q_E[2, i])
            B[1,:] = (-self.q_E[3, i],  self.q_E[0, i],  self.q_E[1, i])
            B[2,:] = ( self.q_E[2, i], -self.q_E[1, i],  self.q_E[0, i])
            B[3,:] = ( self.q_E[1, i],  self.q_E[2, i],  self.q_E[3, i])

            self.O_IE[:,:, i] = np.dot(A.T, B)

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'q_E' in arg and 'O_IE' in result:
            for u in range(3):
                for v in range(3):
                    for k in range(4):
                        result['O_IE'][u, v,:] += self.J[:, u, v, k] * \
                            arg['q_E'][k,:]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'O_IE' in arg and 'q_E' in result:
            for u in range(3):
                for v in range(3):
                    for k in range(4):
                        result['q_E'][k,:] += self.J[:, u, v, k] * \
                            arg['O_IE'][u, v,:]


class Comm_GainPattern(Component):

    ''' Determines transmitter gain based on an external az-el map. '''

    def __init__(self, n, rawG=None):
        super(Comm_GainPattern, self).__init__()
        self.n = n

        if rawG is None:
            fpath = os.path.dirname(os.path.realpath(__file__))
            rawGdata = np.genfromtxt(fpath + '/data/Comm/Gain.txt')
            rawG = (10 ** (rawGdata / 10.0)).reshape((361, 361), order='F')

        # Inputs
        self.add(
            'azimuthGS',
            Array(
                np.zeros(n),
                iotype='in',
                shape=(n,),
                units="rad",
                desc="Azimuth angle from satellite to ground station in Earth-fixed frame over time"
            )
        )

        self.add(
            'elevationGS',
            Array(
                np.zeros(n),
                iotype='in',
                shape=(self.n,),
                units="rad",
                desc="Elevation angle from satellite to ground station in Earth-fixed frame over time"
            )
        )

        # Outputs
        self.add('gain', Array(np.zeros(n),
                               iotype='out',
                               shape=(n,),
                               units="unitless",
                               desc="Transmitter gain over time"))

        pi = np.pi
        az = np.linspace(0, 2 * pi, 361)
        el = np.linspace(0, 2 * pi, 361)

        self.MBI = MBI.MBI(rawG, [az, el], [15, 15], [4, 4])
        self.x = np.zeros((self.n, 2), order='F')

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.dg_daz = self.MBI.evaluate(self.x, 1)[:, 0]
        self.dg_del = self.MBI.evaluate(self.x, 2)[:, 0]

    def execute(self):
        """ Calculate output. """

        result = fixangles(self.n, self.azimuthGS, self.elevationGS)
        self.x[:, 0] = result[0]
        self.x[:, 1] = result[1]
        self.gain = self.MBI.evaluate(self.x)[:, 0]

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'azimuthGS' in arg and 'gain' in result:
            result['gain'] += self.dg_daz * arg['azimuthGS']
        if 'elevationGS' in arg and 'gain' in result:
            result['gain'] += self.dg_del * arg['elevationGS']

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'azimuthGS' in result and 'gain' in arg:
            result['azimuthGS'] += self.dg_daz * arg['gain']
        if 'elevationGS' in result and 'gain' in arg:
            result['elevationGS'] += self.dg_del * arg['gain']


class Comm_GSposEarth(Component):

    ''' Returns position of the ground station in Earth frame. '''

    # constants
    Re = 6378.137
    d2r = np.pi / 180.

    # Inputs
    lon = Float(0.0, iotype="in", units="rad", desc="Longitude of ground station in Earth-fixed frame")
    lat = Float(0.0, iotype="in", units="rad", desc="Latitude of ground station in Earth-fixed frame")
    alt = Float(0.0, iotype="in", units="rad", desc="Altitude of ground station in Earth-fixed frame")

    def __init__(self, n):
        super(Comm_GSposEarth, self).__init__()

        self.n = n

        # Outputs
        self.add(
            'r_e2g_E',
            Array(
                np.zeros((3, self.n)),
                iotype='out',
                shape=(3, self.n),
                units="km",
                desc="Position vector from earth to ground station in Earth-fixed frame over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.dr_dlon = np.zeros(3)
        self.dr_dlat = np.zeros(3)
        self.dr_dalt = np.zeros(3)

        cos_lat = np.cos(self.d2r * self.lat)
        sin_lat = np.sin(self.d2r * self.lat)
        cos_lon = np.cos(self.d2r * self.lon)
        sin_lon = np.sin(self.d2r * self.lon)

        r_GS = (self.Re + self.alt)

        self.dr_dlon[0] = -self.d2r * r_GS * cos_lat * sin_lon
        self.dr_dlat[0] = -self.d2r * r_GS * sin_lat * cos_lon
        self.dr_dalt[0] = cos_lat * cos_lon

        self.dr_dlon[1] = self.d2r * r_GS * cos_lat * cos_lon
        self.dr_dlat[1] = -self.d2r * r_GS * sin_lat * sin_lon
        self.dr_dalt[1] = cos_lat * sin_lon

        self.dr_dlon[2] = 0.
        self.dr_dlat[2] = self.d2r * r_GS * cos_lat
        self.dr_dalt[2] = sin_lat

    def execute(self):
        """ Calculate output. """

        cos_lat = np.cos(self.d2r * self.lat)
        r_GS = (self.Re + self.alt)

        self.r_e2g_E[0,:] = r_GS * cos_lat * np.cos(self.d2r*self.lon)
        self.r_e2g_E[1,:] = r_GS * cos_lat * np.sin(self.d2r*self.lon)
        self.r_e2g_E[2,:] = r_GS * np.sin(self.d2r*self.lat)

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'lon' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] += self.dr_dlon[k] * arg['lon']
        if 'lat' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] += self.dr_dlat[k] * arg['lat']
        if 'alt' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] += self.dr_dalt[k] * arg['alt']

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'r_e2g_E' in arg:
            for k in xrange(3):
                if 'lon' in result:
                    result['lon'] += self.dr_dlon[k] * np.sum(arg['r_e2g_E'][k,:])
                if 'lat' in result:
                    result['lat'] += self.dr_dlat[k] * np.sum(arg['r_e2g_E'][k,:])
                if 'alt' in result:
                    result['alt'] += self.dr_dalt[k] * np.sum(arg['r_e2g_E'][k,:])


class Comm_GSposECI(Component):

    ''' Convert time history of ground station position from earth frame
    to inertial frame
    '''

    def __init__(self, n):
        super(Comm_GSposECI, self).__init__()
        self.n = n

        # Inputs
        self.add(
            'O_IE', 
            Array(
                np.zeros((3, 3, self.n)),
                iotype='in',
                shape=(3, 3, self.n),
                units="unitless",
                desc="Rotation matrix from Earth-centered inertial frame to Earth-fixed frame over time"
            )
        )

        self.add(
            'r_e2g_E',
            Array(
                np.zeros((3, self.n)),
                iotype='in',
                shape=(3, self.n),
                units="km",
                desc="Position vector from earth to ground station in Earth-fixed frame over time"
            )
        )

        # Outputs
        self.add(
            'r_e2g_I',
            Array(
                np.zeros((3, self.n)),
                iotype='out',
                shape=(3, self.n),
                units="km",
                desc="Position vector from earth to ground station in Earth-centered inertial frame over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.J1 = np.zeros((self.n, 3, 3, 3))

        for k in range(0, 3):
            for v in range(0, 3):
                self.J1[:, k, k, v] = self.r_e2g_E[v,:]

        self.J2 = np.transpose(self.O_IE, (2, 0, 1))

    def execute(self):
        """ Calculate output. """

        for i in range(0, self.n):
            self.r_e2g_I[:, i] = np.dot(self.O_IE[:,:, i],
                                        self.r_e2g_E[:, i])

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'r_e2g_I' in result:
            for k in xrange(3):
                for u in xrange(3):
                    if 'O_IE' in arg:
                        for v in xrange(3):
                            result['r_e2g_I'][k,:] += self.J1[:, k, u, v] * \
                                arg['O_IE'][u, v,:]

                    if 'r_e2g_E' in arg:
                        result['r_e2g_I'][k,:] += self.J2[:, k, u] * \
                            arg['r_e2g_E'][u,:]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'r_e2g_I' in arg:
            for k in xrange(3):
                if 'O_IE' in result:
                    for u in xrange(3):
                        for v in xrange(3):
                            result['O_IE'][u, v,:] += self.J1[:, k, u, v] * \
                                arg['r_e2g_I'][k,:]
                if 'r_e2g_E' in result:
                    for j in xrange(3):
                        result['r_e2g_E'][j,:] += self.J2[:, k, j] * \
                            arg['r_e2g_I'][k,:]


class Comm_LOS(Component):

    ''' Determines if the Satellite has line of sight with the ground
    stations. '''

    # constants
    Re = 6378.137

    def __init__(self, n):
        super(Comm_LOS, self).__init__()
        self.n = n

        # Inputs
        self.add(
            'r_b2g_I', 
            Array(
                np.zeros((3, n)),
                iotype='in',
                shape=(3, self.n),
                units="km",
                desc="Position vector from satellite to ground station in Earth-centered inertial frame over time"
            )
        )

        self.add(
            'r_e2g_I',
            Array(
                np.zeros((3, n)),
                iotype='in',
                shape=(3, self.n),
                units="km",
                desc="Position vector from earth to ground station in Earth-centered inertial frame over time"
            )
        )

        # Outputs
        self.add(
            'CommLOS',
            Array(
                np.zeros(n),
                iotype='out',
                shape=(self.n, ),
                units="unitless",
                desc="Satellite to ground station line of sight over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.dLOS_drb = np.zeros((self.n, 3))
        self.dLOS_dre = np.zeros((self.n, 3))

        Rb = 10.0
        for i in range(0, self.n):

            proj = np.dot(self.r_b2g_I[:, i], self.r_e2g_I[:, i]) / self.Re

            if proj > 0:
                self.dLOS_drb[i,:] = 0.
                self.dLOS_dre[i,:] = 0.
            elif proj < -Rb:
                self.dLOS_drb[i,:] = 0.
                self.dLOS_dre[i,:] = 0.
            else:
                x = (proj - 0) / (-Rb - 0)
                dx_dproj = -1. / Rb
                dLOS_dx = 6 * x - 6 * x ** 2
                dproj_drb = self.r_e2g_I[:, i]
                dproj_dre = self.r_b2g_I[:, i]

                self.dLOS_drb[i,:] = dLOS_dx * dx_dproj * dproj_drb
                self.dLOS_dre[i,:] = dLOS_dx * dx_dproj * dproj_dre

    def execute(self):
        """ Calculate output. """

        Rb = 100.0
        for i in range(0, self.n):
            proj = np.dot(self.r_b2g_I[:, i], self.r_e2g_I[:, i]) / self.Re

            if proj > 0:
                self.CommLOS[i] = 0.
            elif proj < -Rb:
                self.CommLOS[i] = 1.
            else:
                x = (proj - 0) / (-Rb - 0)
                self.CommLOS[i] = 3 * x ** 2 - 2 * x ** 3

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'CommLOS' in result:
            for k in xrange(3):
                if 'r_b2g_I' in arg:
                    result['CommLOS'] += self.dLOS_drb[:, k] * arg['r_b2g_I'][k,:]
                if 'r_e2g_I' in arg:
                    result['CommLOS'] += self.dLOS_dre[:, k] * arg['r_e2g_I'][k,:]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'CommLOS' in arg:
            for k in xrange(3):
                if 'r_b2g_I' in result:
                    result['r_b2g_I'][k,:] += self.dLOS_drb[:, k] * arg['CommLOS']
                if 'r_e2g_I' in result:
                    result['r_e2g_I'][k,:] += self.dLOS_dre[:, k] * arg['CommLOS']


class Comm_VectorAnt(Component):

    '''Transform from antenna to body frame'''

    def __init__(self, n):
        super(Comm_VectorAnt, self).__init__()
        self.n = n

        # Inputs
        self.add(
            'r_b2g_B', 
            Array(
                np.zeros((3, n)),
                iotype='in',
                shape=(3, n),
                units="km",
                desc="Position vector from satellite to ground station in body-fixed frame over time"
            )
        )

        self.add(
            'O_AB',
            Array(
                np.zeros((3, 3, n)),
                iotype='in',
                shape=(3, 3, n),
                units="unitless",
                desc="Rotation matrix from antenna angle to body-fixed frame over time"
            )
        )

        # Outputs
        self.add(
            'r_b2g_A',
            Array(
                np.zeros((3, n)),
                iotype='out',
                shape=(3, n),
                units="km",
                desc="Position vector from satellite to ground station in antenna angle frame over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.J1, self.J2 = computepositionrotdjacobian(self.n, self.r_b2g_B,
                                                       self.O_AB)

    def execute(self):
        """ Calculate output. """

        self.r_b2g_A = computepositionrotd(self.n, self.r_b2g_B, self.O_AB)

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'r_b2g_A' in result:
            for k in xrange(3):
                if 'O_AB' in arg:
                    for u in xrange(3):
                        for v in xrange(3):
                            result['r_b2g_A'][k,:] += self.J1[:, k, u, v] * \
                                arg['O_AB'][u, v,:]
                if 'r_b2g_B' in arg:
                    for j in xrange(3):
                        result['r_b2g_A'][k,:] += self.J2[:, k, j] * \
                            arg['r_b2g_B'][j,:]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'r_b2g_A' in arg:
            for k in xrange(3):
                if 'O_AB' in result:
                    for u in xrange(3):
                        for v in xrange(3):
                            result['O_AB'][u, v,:] += self.J1[:, k, u, v] * \
                                arg['r_b2g_A'][k,:]
                if 'r_b2g_B' in result:
                    for j in xrange(3):
                        result['r_b2g_B'][j,:] += self.J2[:, k, j] * \
                            arg['r_b2g_A'][k,:]


class Comm_VectorBody(Component):

    '''Transform from body to inertial frame'''

    def __init__(self, n):
        super(Comm_VectorBody, self).__init__()
        self.n = n

        # Inputs
        self.add(
            'r_b2g_I',
            Array(
                np.zeros((3, n)),
                iotype='in',
                shape=(3, n),
                units="km",
                desc="Position vector from satellite to ground station in Earth-centered inertial frame over time"
            )
        )

        self.add(
            'O_BI',
            Array(
                np.zeros((3, 3, n)),
                iotype='in',
                shape=(3, 3, n),
                units="unitless",
                desc="Rotation matrix from body-fixed frame to Earth-centered inertial frame over time"
            )
        )

        # Outputs
        self.add(
            'r_b2g_B',
            Array(
                np.zeros((3, n)),
                iotype='out',
                shape=(3, n),
                units="km",
                desc="Position vector from satellite to ground station in body-fixed frame over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.J1 = np.zeros((self.n, 3, 3, 3))

        for k in range(0, 3):
            for v in range(0, 3):
                self.J1[:, k, k, v] = self.r_b2g_I[v,:]

        self.J2 = np.transpose(self.O_BI, (2, 0, 1))

    def execute(self):
        """ Calculate output. """

        for i in range(0, self.n):
            self.r_b2g_B[:, i] = np.dot(self.O_BI[:,:, i], self.r_b2g_I[:, i])

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'r_b2g_B' in result:
            for k in range(3):
                if 'O_BI' in arg:
                    for u in range(3):
                        for v in range(3):
                            result['r_b2g_B'][k,:] += self.J1[:, k, u, v] * \
                                arg['O_BI'][u, v,:]
                if 'r_b2g_I' in arg:
                    for j in range(3):
                        result['r_b2g_B'][k,:] += self.J2[:, k, j] * \
                            arg['r_b2g_I'][j,:]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'r_b2g_B' in arg:
            for k in range(3):
                if 'O_BI' in result:
                    for u in range(3):
                        for v in range(3):
                            result['O_BI'][u, v,:] += self.J1[:, k, u, v] * \
                                arg['r_b2g_B'][k,:]
                if 'r_b2g_I' in result:
                    for j in range(3):
                        result['r_b2g_I'][j,:] += self.J2[:, k, j] * \
                            arg['r_b2g_B'][k,:]


class Comm_VectorECI(Component):

    '''Determine vector between satellite and ground station.'''

    def __init__(self, n):
        super(Comm_VectorECI, self).__init__()
        self.n = n

        # Inputs
        self.add(
            'r_e2g_I',
            Array(
                np.zeros((3, n)),
                iotype='in',
                shape=(3, n),
                units="km",
                desc="Position vector from earth to ground station in Earth-centered inertial frame over time"
            )
        )

        self.add(
            'r_e2b_I',
            Array(
                np.zeros((6, n)),
                iotype='in',
                shape=(6, n),
                units="unitless",
                desc="Position and velocity vector from earth to satellite in Earth-centered inertial frame over time"
            )
        )

        # Outputs
        self.add(
            'r_b2g_I',
            Array(
                np.zeros((3, n)),
                iotype='out',
                shape=(3, n),
                units="km",
                desc="Position vector from satellite to ground station in Earth-centered inertial frame over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """
        # Derivatives are simple
        return

    def execute(self):
        """ Calculate output. """

        self.r_b2g_I = self.r_e2g_I - self.r_e2b_I[:3,:]

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'r_e2g_I' in arg:
            result['r_b2g_I'] += arg['r_e2g_I']
        if 'r_e2b_I' in arg:
            result['r_b2g_I'] += -arg['r_e2b_I'][:3,:]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'r_b2g_I' in arg:
            if 'r_e2g_I' in result:
                result['r_e2g_I'] += arg['r_b2g_I']
            if 'r_e2b_I' in result:
                result['r_e2b_I'][:3,:] += -arg['r_b2g_I']


class Comm_VectorSpherical(Component):

    '''Convert satellite-ground vector into Az-El.'''

    def __init__(self, n):
        super(Comm_VectorSpherical, self).__init__()
        self.n = n

        # Inputs
        self.add(
            'r_b2g_A',
            Array(
                np.zeros((3, n)),
                iotype='in',
                shape=(3, self.n),
                units="km",
                desc="Position vector from satellite to ground station in antenna angle frame over time"
            )
        )

        # Outputs
        self.add(
            'azimuthGS',
            Array(
                np.zeros(n), 
                iotype='out',
                shape=(n,),
                units="rad",
                desc="Azimuth angle from satellite to ground station in Earth-fixed frame over time"
            )
        )

        self.add(
            'elevationGS',
            Array(
                np.zeros(n),
                iotype='out',
                shape=(n,),
                units="rad",
                desc="Elevation angle from satellite to ground station in Earth-fixed frame over time"
            )
        )

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.Ja1, self.Ji1, self.Jj1, self.Ja2, self.Ji2, self.Jj2 = \
            computepositionsphericaljacobian(self.n, 3 * self.n, self.r_b2g_A)

        self.J1 = scipy.sparse.csc_matrix((self.Ja1, (self.Ji1, self.Jj1)),
                                          shape=(self.n, 3 * self.n))
        self.J2 = scipy.sparse.csc_matrix((self.Ja2, (self.Ji2, self.Jj2)),
                                          shape=(self.n, 3 * self.n))
        self.J1T = self.J1.transpose()
        self.J2T = self.J2.transpose()

    def execute(self):
        """ Calculate output. """

        azimuthGS, elevationGS = computepositionspherical(self.n, self.r_b2g_A)
        self.azimuthGS = azimuthGS
        self.elevationGS = elevationGS

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'r_b2g_A' in arg:
            r_b2g_A = arg['r_b2g_A'].reshape((3 * self.n), order='F')
            if 'azimuthGS' in result:
                result['azimuthGS'] += self.J1.dot(r_b2g_A)
            if 'elevationGS' in result:
                result['elevationGS'] += self.J2.dot(r_b2g_A)

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'azimuthGS' in arg:
            az_GS = arg['azimuthGS']
            result['r_b2g_A'] += (self.J1T.dot(az_GS)).reshape((3, self.n),
                                                               order='F')
        if 'elevationGS' in arg:
            el_GS = arg['elevationGS']
            result['r_b2g_A'] += (self.J2T.dot(el_GS)).reshape((3, self.n),
                                                               order='F')