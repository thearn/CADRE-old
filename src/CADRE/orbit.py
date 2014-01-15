''' Orbit discipline for CADRE '''

from math import sqrt

import numpy as np

from openmdao.lib.datatypes.api import Float, Array
from openmdao.main.api import Component

import rk4

# Allow non-standard variable names for scientific calc
# pylint: disable-msg=C0103

# Constants
mu = 398600.44
Re = 6378.137
J2 = 1.08264e-3
J3 = -2.51e-6
J4 = -1.60e-6

C1 = -mu
C2 = -1.5*mu*J2*Re**2
C3 = -2.5*mu*J3*Re**3
C4 = 1.875*mu*J4*Re**4

class Orbit_Dynamics(rk4.RK4):
    """Computes the Earth to body position vector in Earth-centered intertial frame."""

    def __init__(self, n_times):
        super(Orbit_Dynamics, self).__init__()

        # Inputs
        self.add('r_e2b_I0',Array(np.zeros((6,)),
                            size=(6,), iotype="in",
                            dtype=np.float,
                            fd_step=1e-2,
                            units="unitless",
                            desc="Initial position and velocity vectors from earth to satellite in Earth-centered inertial frame"))

        # Outputs
        self.add('r_e2b_I', Array(1000*np.ones((6, n_times)),
                                  size=(6, n_times),
                                  dtype=np.float,
                                  iotype="out",
                                  units="unitless",
                                  desc="Position and velocity vectors from earth to satellite in Earth-centered inertial frame over time"))

        self.state_var = 'r_e2b_I'
        self.init_state_var = 'r_e2b_I0'

        self.dfdx = np.zeros((6, 1))

    def list_deriv_vars(self):
        input_keys = ('r_e2b_I0',)
        output_keys = ('r_e2b_I',)
        return input_keys, output_keys

    def f_dot(self, external, state):

        x = state[0]
        y = state[1]
        z = state[2] if abs(state[2]) > 1e-15 else 1e-5

        z2 = z*z
        z3 = z2*z
        z4 = z3*z

        r = sqrt(x*x + y*y + z2)

        r2 = r*r
        r3 = r2*r
        r4 = r3*r
        r5 = r4*r
        r7 = r5*r*r

        T2 = 1 - 5*z2/r2
        T3 = 3*z - 7*z3/r2
        T4 = 1 - 14*z2/r2 + 21*z4/r4
        T3z = 3*z - 0.6*r2/z
        T4z = 4 - 28.0/3.0*z2/r2

        f_dot = np.zeros((6,))
        f_dot[0:3] = state[3:]
        f_dot[3:] = state[0:3]*(C1/r3 + C2/r5*T2 + C3/r7*T3 + C4/r7*T4)
        f_dot[5] += z*(2.0*C2/r5 + C3/r7*T3z + C4/r7*T4z)

        return f_dot

    def df_dy(self, external, state):

        x = state[0]
        y = state[1]
        z = state[2] if abs(state[2]) > 1e-15 else 1e-5

        z2 = z*z
        z3 = z2*z
        z4 = z3*z

        r = sqrt(x*x + y*y + z2)

        r2 = r*r
        r3 = r2*r
        r4 = r3*r
        r5 = r4*r
        r6 = r5*r
        r7 = r6*r
        r8 = r7*r

        drdx = x/r
        drdy = y/r
        drdz = z/r

        T2 = 1 - 5*z2/r2
        T3 = 3*z - 7*z3/r2
        T4 = 1 - 14*z2/r2 + 21*z4/r4
        T3z = 3*z - 0.6*r2/z
        T4z = 4 - 28.0/3.0*z2/r2

        dT2_dx = (10*z2)/(r3)*drdx
        dT2_dy = (10*z2)/(r3)*drdy
        dT2_dz = (10*z2)/(r3)*drdz - 10.*z/r2

        dT3_dx = 14*z3/r3*drdx
        dT3_dy = 14*z3/r3*drdy
        dT3_dz = 14*z3/r3*drdz - 21.*z2/r2 + 3

        dT4_dx = 28*z2/r3*drdx - 84.*z4/r5*drdx
        dT4_dy = 28*z2/r3*drdy - 84.*z4/r5*drdy
        dT4_dz = 28*z2/r3*drdz - 84.*z4/r5*drdz - 28*z/r2 + 84*z3/r4

        dT3z_dx = -1.2*r/z*drdx
        dT3z_dy = -1.2*r/z*drdy
        dT3z_dz = -1.2*r/z*drdz + 0.6*r2/z2 + 3

        dT4z_dx = 56.0/3.0*z2/r3*drdx
        dT4z_dy = 56.0/3.0*z2/r3*drdy
        dT4z_dz = 56.0/3.0*z2/r3*drdz - 56.0/3.0*z/r2

        eye = np.identity(3)

        dfdy = np.zeros((6, 6))
        dfdy[0:3, 3:] += eye

        dfdy[3:, :3] = dfdy[3:, :3] + eye*(C1/r3 + C2/r5*T2 + C3/r7*T3 + C4/r7*T4)
        dfdy[3:, 0] = dfdy[3:, 0] + drdx*state[:3]*(-3*C1/r4 - 5*C2/r6*T2 - 7*C3/r8*T3 - 7*C4/r8*T4)
        dfdy[3:, 1] = dfdy[3:, 1] + drdy*state[:3]*(-3*C1/r4 - 5*C2/r6*T2 - 7*C3/r8*T3 - 7*C4/r8*T4)
        dfdy[3:, 2] = dfdy[3:, 2] + drdz*state[:3]*(-3*C1/r4 - 5*C2/r6*T2 - 7*C3/r8*T3 - 7*C4/r8*T4)
        dfdy[3:, 0] = dfdy[3:, 0] + state[:3]*(C2/r5*dT2_dx + C3/r7*dT3_dx + C4/r7*dT4_dx)
        dfdy[3:, 1] = dfdy[3:, 1] + state[:3]*(C2/r5*dT2_dy + C3/r7*dT3_dy + C4/r7*dT4_dy)
        dfdy[3:, 2] = dfdy[3:, 2] + state[:3]*(C2/r5*dT2_dz + C3/r7*dT3_dz + C4/r7*dT4_dz)
        dfdy[5, 0] = dfdy[5, 0] + drdx*z*(-5*C2/r6*2 - 7*C3/r8*T3z - 7*C4/r8*T4z)
        dfdy[5, 1] = dfdy[5, 1] + drdy*z*(-5*C2/r6*2 - 7*C3/r8*T3z - 7*C4/r8*T4z)
        dfdy[5, 2] = dfdy[5, 2] + drdz*z*(-5*C2/r6*2 - 7*C3/r8*T3z - 7*C4/r8*T4z)
        dfdy[5, 0] = dfdy[5, 0] + z*(C3/r7*dT3z_dx + C4/r7*dT4z_dx)
        dfdy[5, 1] = dfdy[5, 1] + z*(C3/r7*dT3z_dy + C4/r7*dT4z_dy)
        dfdy[5, 2] = dfdy[5, 2] + z*(C3/r7*dT3z_dz + C4/r7*dT4z_dz)
        dfdy[5, 2] = dfdy[5, 2] + (C2/r5*2 + C3/r7*T3z + C4/r7*T4z)
        #dfdy[5, 2] = dfdy[5, 2] + (2.0*C2 + (C3*T3z + C4*T4z)/r2)/r5

        #print dfdy
        return dfdy

    def df_dx(self, external, state):

        return self.dfdx


class Orbit_Initial(Component):
    """Computes initial position and velocity vectors of Earth to body position"""

    # Inputs
    altPerigee = Float(500., iotype="in", copy=None)
    altApogee = Float(500., iotype="in", copy=None)
    RAAN = Float(66.279, iotype="in", copy=None)
    Inc = Float(82.072, iotype="in", copy=None)
    argPerigee = Float(0, iotype="in", copy=None)
    trueAnomaly = Float(337.987, iotype="in", copy=None)

    def __init__(self):
        super(Orbit_Initial, self).__init__()

        #Outputs
        self.add('r_e2b_I0', Array(np.ones((6,)),
                                   size=(6,),
                                   dtype=np.float,
                                   iotype='out',
                                   units="unitless",
                                   desc="Initial position and velocity vectors from Earth to satellite in Earth-centered inertial frame"))

    def list_deriv_vars(self):
        input_keys = ('altPerigee','altApogee','RAAN','Inc','argPerigee',
                      'trueAnomaly')
        output_keys = ('r_e2b_I0',)
        return input_keys, output_keys

    def compute(self, altPerigee, altApogee, RAAN, Inc, argPerigee, trueAnomaly):
        ''' Compute position and velocity from orbital elements '''

        Re = 6378.137
        mu = 398600.44

        def S(v):
            S = np.zeros((3,3),complex)
            S[0,:] = [0, -v[2], v[1]]
            S[1,:] = [v[2], 0, -v[0]]
            S[2,:] = [-v[1], v[0], 0]
            return S

        def getRotation(axis, angle):
            R = np.eye(3,dtype=complex) + S(axis)*np.sin(angle) + \
                (1 - np.cos(angle)) * (np.outer(axis, axis) - np.eye(3, dtype=complex))
            return R

        d2r = np.pi/180.0
        r_perigee = Re + altPerigee
        r_apogee = Re + altApogee
        e = (r_apogee-r_perigee)/(r_apogee+r_perigee)
        a = (r_perigee+r_apogee)/2
        p = a*(1-e**2)
        h = np.sqrt(p*mu)

        rmag0 = p/(1+e*np.cos(d2r*trueAnomaly))
        r0_P = np.array([rmag0*np.cos(d2r*trueAnomaly), rmag0*np.sin(d2r*trueAnomaly), 0], complex)
        v0_P = np.array([-np.sqrt(mu/p)*np.sin(d2r*trueAnomaly), np.sqrt(mu/p)*(e+np.cos(d2r*trueAnomaly)), 0], complex)

        O_IP = np.eye(3, dtype=complex)
        O_IP = np.dot(O_IP, getRotation(np.array([0,0,1]),RAAN*d2r))
        O_IP = np.dot(O_IP, getRotation(np.array([1,0,0]),Inc*d2r))
        O_IP = np.dot(O_IP, getRotation(np.array([0,0,1]),argPerigee*d2r))

        r0_ECI = np.dot(O_IP, r0_P)
        v0_ECI = np.dot(O_IP, v0_P)

        return r0_ECI, v0_ECI

    def linearize(self):
        """ Calculate and save derivatives, (i.e., Jacobian) """

        h = 1e-16
        ih = complex(0, h)
        v = np.zeros(6, complex)
        v[:] = [self.altPerigee, self.altApogee, self.RAAN, self.Inc,
                self.argPerigee, self.trueAnomaly]
        self.J = np.zeros((6,6))

        for i in range(6):
            v[i] += ih
            r0_ECI, v0_ECI = self.compute(v[0], v[1], v[2], v[3], v[4], v[5])
            v[i] -= ih
            self.J[:3,i] = r0_ECI.imag/h
            self.J[3:,i] = v0_ECI.imag/h

    def execute(self):
        """ Calculate output. """

        r0_ECI, v0_ECI = self.compute(self.altPerigee, self.altApogee,
                                      self.RAAN, self.Inc, self.argPerigee,
                                      self.trueAnomaly)
        self.r_e2b_I0[:3] = r0_ECI.real
        self.r_e2b_I0[3:] = v0_ECI.real

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian """

        J = self.J

        if 'r_e2b_I0' in result:
            if 'altPerigee' in arg:
                result['r_e2b_I0'][:3] += J[:3, 0]*arg['altPerigee']
                result['r_e2b_I0'][3:] += J[3:, 0]*arg['altPerigee']
            if 'altApogee' in arg:
                result['r_e2b_I0'][:3] += J[:3, 1]*arg['altApogee']
                result['r_e2b_I0'][3:] += J[3:, 1]*arg['altApogee']
            if 'RAAN' in arg:
                result['r_e2b_I0'][:3] += J[:3, 2]*arg['RAAN']
                result['r_e2b_I0'][3:] += J[3:, 2]*arg['RAAN']
            if 'Inc' in arg:
                result['r_e2b_I0'][:3] += J[:3, 3]*arg['Inc']
                result['r_e2b_I0'][3:] += J[3:, 3]*arg['Inc']
            if 'argPerigee' in arg:
                result['r_e2b_I0'][:3] += J[:3, 4]*arg['argPerigee']
                result['r_e2b_I0'][3:] += J[3:, 4]*arg['argPerigee']
            if 'trueAnomaly' in arg:
                result['r_e2b_I0'][:3] += J[:3, 5]*arg['trueAnomaly']
                result['r_e2b_I0'][3:] += J[3:, 5]*arg['trueAnomaly']


    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian """

        J = self.J

        if 'r_e2b_I0' in arg:
            r_e2b_I0 = arg['r_e2b_I0']
            r_e2b_I0_0 = r_e2b_I0[:3]
            r_e2b_I0_3 = r_e2b_I0[3:]


            if 'altPerigee' in result:
                result['altPerigee'] += sum(J[:3, 0]*r_e2b_I0_0) + \
                                        sum(J[3:, 0]*r_e2b_I0_3)
            if 'altApogee' in result:
                result['altApogee'] += sum(J[:3, 1]*r_e2b_I0_0) + \
                                       sum(J[3:, 1]*r_e2b_I0_3)
            if 'RAAN' in result:
                result['RAAN'] += sum(J[:3, 2]*r_e2b_I0_0) + \
                                  sum(J[3:, 2]*r_e2b_I0_3)
            if 'Inc' in result:
                result['Inc'] += sum(J[:3, 3]*r_e2b_I0_0) + \
                                 sum(J[3:, 3]*r_e2b_I0_3)
            if 'argPerigee' in result:
                result['argPerigee'] += sum(J[:3, 4]*r_e2b_I0_0) + \
                                        sum(J[3:, 4]*r_e2b_I0_3)
            if 'trueAnomaly' in result:
                result['trueAnomaly'] += sum(J[:3, 5]*r_e2b_I0_0) + \
                                         sum(J[3:, 5]*r_e2b_I0_3)
