from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import numpy as np

import rk4


class ReactionWheel_Motor(Component):
    '''Compute motor torque'''
    
    def __init__(self, n):
        super(ReactionWheel_Motor, self).__init__()
        self.n = n

        self.add('J_RW', 2.8e-5)

        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n),
                               units="N*m",
                               desc="Torque vector of reaction wheel over time",
                               dtype=np.float, iotype='in'))
        self.add('w_B', Array(np.zeros((3,n)), size=(3,n),
                              units="1/s",
                              desc="Angular velocity vector in body-fixed frame over time",
                              dtype=np.float, iotype='in'))
        self.add('w_RW', Array(np.zeros((3,n)), size=(3,n),
                               units="1/s",
                               desc="Angular velocity vector of reaction wheel over time",
                               dtype=np.float, iotype='in'))

        self.add('T_m', Array(np.ones((3,n)), size=(3,n),
                              units="N*m",
                              desc="Torque vector of motor over time",
                              dtype=np.float, iotype='out'))

    def linearize(self):
        w_Bx = np.zeros((3,3))
        self.dT_dTm = np.zeros((self.n,3,3))
        self.dT_dwb = np.zeros((self.n,3,3))
        self.dT_dh = np.zeros((self.n,3,3))
        dwx_dwb = np.zeros((3,3,3))
        h_RW = self.J_RW * self.w_RW[:]
        for i in range(0,self.n):
            w_Bx[0,:] = (0., -self.w_B[2,i] , self.w_B[1,i])
            w_Bx[1,:] = (self.w_B[2,i], 0., -self.w_B[0,i])
            w_Bx[2,:] = (-self.w_B[1,i], self.w_B[0,i], 0.)

            dwx_dwb[0,:,0] = (0., 0., 0.)
            dwx_dwb[1,:,0] = (0., 0., -1.)
            dwx_dwb[2,:,0] = (0., 1., 0.)

            dwx_dwb[0,:,1] = (0., 0., 1.)
            dwx_dwb[1,:,1] = (0., 0., 0.)
            dwx_dwb[2,:,1] = (-1., 0., 0.)

            dwx_dwb[0,:,2] = (0., -1., 0.)
            dwx_dwb[1,:,2] = (1., 0., 0.)
            dwx_dwb[2,:,2] = (0., 0., 0.)

            for k in range(0,3):
                self.dT_dTm[i,k,k] = -1.
                self.dT_dwb[i,:,k] = -np.dot(dwx_dwb[:,:,k] , h_RW[:,i])
            self.dT_dh[i,:,:] = -w_Bx

    def execute(self):
        w_Bx = np.zeros((3,3))
        h_RW = self.J_RW * self.w_RW[:]
        for i in range(0,self.n):
            w_Bx[0,:] = (0., -self.w_B[2,i], self.w_B[1,i])
            w_Bx[1,:] = (self.w_B[2,i], 0., -self.w_B[0,i])
            w_Bx[2,:] = (-self.w_B[1,i], self.w_B[0,i], 0.)

            self.T_m[:,i] = -self.T_RW[:,i] - np.dot(w_Bx , h_RW[:,i])

    def apply_deriv(self, arg, result):

        if 'T_m' in result:
            for k in range(3):
                for j in range(3):
                    if 'T_RW' in arg:
                        result['T_m'][k,:] += self.dT_dTm[:,k,j] * arg['T_RW'][j,:]
                    if 'w_B' in arg:
                        result['T_m'][k,:] += self.dT_dwb[:,k,j] * arg['w_B'][j,:]
                    if 'w_RW' in arg:
                        result['T_m'][k,:] += self.dT_dh[:,k,j] * arg['w_RW'][j,:] * self.J_RW

    def apply_derivT(self, arg, result):

        if 'T_m' in arg:
            for k in range(3):
                for j in range(3):
                    if 'T_RW' in result:
                        result['T_RW'][j,:] += self.dT_dTm[:,k,j] * arg['T_m'][k,:]
                    if 'w_B' in result:
                        result['w_B'][j,:] += self.dT_dwb[:,k,j] * arg['T_m'][k,:]
                    if 'w_RW' in result:
                        result['w_RW'][j,:] += self.dT_dh[:,k,j] * arg['T_m'][k,:] * self.J_RW


class ReactionWheel_Power(Component):
    '''Compute reaction wheel power'''
    
    #constants
    V = 4.0
    a = 4.9e-4
    b = 4.5e2
    I0 = 0.017

    def __init__(self, n):
        super(ReactionWheel_Power, self).__init__()
        self.n = n

        self.add('w_RW', Array(np.zeros((3,n)), size=(3,n),
                               units="1/s",
                               desc="Angular velocity vector of reaction wheel over time",
                               dtype=np.float, iotype='in'))
        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n),
                               units="N*m",
                               desc="Torque vector of reaction wheel over time",
                               dtype=np.float, iotype='in'))

        self.add('P_RW', Array(np.ones((3,n)), size=(3,n),
                               units='W',
                               desc='Reaction wheel power over time',
                               dtype=np.float, iotype='out'))

    def linearize(self):
        self.dP_dw = np.zeros((self.n,3))
        self.dP_dT = np.zeros((self.n,3))
        for i in range(self.n):
            for k in range(3):
                prod = 2 * self.V * (self.a * self.w_RW[k,i] + self.b * self.T_RW[k,i])
                self.dP_dw[i,k] = self.a * prod
                self.dP_dT[i,k] = self.b * prod

    def execute(self):
        for i in range(self.n):
            for k in range(3):
                self.P_RW[k,i] = self.V * (self.a * self.w_RW[k,i] + self.b * self.T_RW[k,i])**2 + self.V * self.I0

    def apply_deriv(self, arg, result):
        
        if 'P_RW' in result:
            for k in range(3):
                if 'w_RW' in arg:
                    result['P_RW'][k,:] += self.dP_dw[:,k] * arg['w_RW'][k,:]
                if 'T_RW' in arg:
                    result['P_RW'][k,:] += self.dP_dT[:,k] * arg['T_RW'][k,:]

    def apply_derivT(self, arg, result):

        if 'P_RW' in arg:
            for k in range(3):
                if 'w_RW' in result:
                    result['w_RW'][k,:] += self.dP_dw[:,k] * arg['P_RW'][k,:]
                if 'T_RW' in result:
                    result['T_RW'][k,:] += self.dP_dT[:,k] * arg['P_RW'][k,:]


class ReactionWheel_Torque(Component):
    '''Compute reaction wheel torque'''
    
    def __init__(self, n):
        super(ReactionWheel_Torque, self).__init__()
        self.n = n

        self.add('T_tot', Array(np.zeros((3,n)), size=(3,n),
                                units='N*m',
                                desc='Total reaction wheel torque over time',
                                dtype=np.float, iotype='in'))

        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n),
                               units="N*m",
                               desc="Torque vector of reaction wheel over time",
                               dtype=np.float, iotype='out'))

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """
        # Derivatives are simple
        return

    def execute(self):
        self.T_RW[:] = self.T_tot[:]

    def apply_deriv(self, arg, result):

        if 'T_tot' in arg and 'T_RW' in result:
            result['T_RW'][:] += arg['T_tot'][:]

    def apply_derivT(self, arg, result):

        if 'T_RW' in arg and 'T_tot' in result:
            result['T_tot'][:] += arg['T_RW'][:]


class ReactionWheel_Dynamics(rk4.RK4):

    def __init__(self, n_times):
        super(ReactionWheel_Dynamics, self).__init__()

        self.add('w_B', Array(np.zeros((3, n_times)), size=(3, n_times),
                              units="1/s",
                              desc="Angular velocity vector in body-fixed frame over time",
                              dtype=np.float, iotype='in'))
        self.add('T_RW', Array(np.zeros((3, n_times)), size=(3, n_times),
                               units="N*m",
                               desc="Torque vector of reaction wheel over time",
                               dtype=np.float, iotype='in'))

        self.add('w_RW', Array(np.zeros((3, n_times)), size=(3, n_times),
                               units="1/s",
                               desc="Angular velocity vector of reaction wheel over time",
                               dtype=np.float, iotype='out'))
        self.add('w_RW0', Array(np.zeros((3,)), size=(3,),
                                units="1/s",
                                desc="Initial angular velocity vector of reaction wheel",
                                dtype=np.float, iotype='in'))

        self.state_var = 'w_RW'
        self.init_state_var = 'w_RW0'
        self.external_vars = ['w_B', 'T_RW']

        self.jy = np.zeros((3, 3))

        self.djy_dx = np.zeros((3, 3, 3))
        self.djy_dx[:,:,0] = [[0,0,0],[0,0,-1],[0,1,0]]
        self.djy_dx[:,:,1] = [[0,0,1],[0,0,0],[-1,0,0]]
        self.djy_dx[:,:,2] = [[0,-1,0],[1,0,0],[0,0,0]]

        self.J_RW = 2.8e-5 #unit conversion of some kind

    def f_dot(self, external, state):

        self.jy[0, :] = [0., -external[2], external[1]]
        self.jy[1, :] = [external[2], 0., -external[0]]
        self.jy[2, :] = [-external[1], external[0], 0.]

        #TODO sort out unit conversion here with T_RW
        return (-external[3:]/2.8e-5 - self.jy.dot(state))


    def df_dy(self, external, state):

        self.jy[0, :] = [0., -external[2], external[1]]
        self.jy[1, :] = [external[2], 0., -external[0]]
        self.jy[2, :] = [-external[1], external[0], 0.]
        return -self.jy

    def df_dx(self, external, state):
        self.jx = np.zeros((3, 6))

        for i in xrange(3):
            self.jx[i, 0:3] = -self.djy_dx[:,:,i].dot(state)
            self.jx[i, i+3] = -1.0 / self.J_RW
        return self.jx

