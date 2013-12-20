''' Sun discipline for CADRE '''

import numpy as np
import scipy.sparse

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

from kinematics import computepositionrotd, computepositionrotdjacobian
from kinematics import computepositionspherical, computepositionsphericaljacobian

class Sun_LOS( Component ):

    '''Compute the Satellite to sun line of sight.'''
    
    def __init__(self, n=2):
        super(Sun_LOS, self).__init__()

        self.n = n

        self.r1 = 6378.137*0.85 # Earth's radius is 6378 km. 0.85 is the alpha in John Hwang's paper
        self.r2 = 6378.137

        self.add('r_e2b_I', Array(np.zeros((6, n), order='F'),
                                  size=(6,n, ), dtype=np.float,
                                  units = "unitless",
                                  desc="Position and velocity vectors from " +
                                  "Earth to satellite in Earth-centered " +
                                  "inertial frame over time.",
                                  iotype="in"))
        
        self.add('r_e2s_I', Array(np.zeros((3, n), order='F'), size=(3,n, ),
                                  dtype=np.float,
                                  units="km", desc="Position vector from " +
                                  "Earth to sun in Earth-centered inertial " +
                                  "frame over time.",
                                  iotype="in"))

        self.add('LOS', Array(np.zeros((n, ), order='F'), size=(n, ),
                              dtype=np.float,
                              units="unitless",
                              desc="Satellite to sun " +
                              "line of sight over time",
                              iotype="out"
                              ))

    def execute(self):

        for i in range( self.n ):
            r_b = self.r_e2b_I[:3,i]
            r_s = self.r_e2s_I[:3,i]
            dot = np.dot( r_b, r_s )
            cross = np.cross( r_b, r_s )
            dist = np.sqrt( cross.dot(cross) )

            if dot >= 0.0 :
                self.LOS[i] = 1.0
            elif dist <= self.r1 :
                self.LOS[i] = 0.0
            elif dist >= self.r2 :
                self.LOS[i] = 1.0
            else :
                x = ( dist - self.r1 ) / ( self.r2 - self.r1 )
                self.LOS[i] = 3 *x ** 2 - 2 * x**3

    def linearize(self):

        nj = 3*self.n

        Jab = np.zeros(shape=(nj,), dtype = np.float)
        Jib = np.zeros(shape=(nj,), dtype = np.int)
        Jjb = np.zeros(shape=(nj,), dtype = np.int)
        Jas = np.zeros(shape=(nj,), dtype = np.float)
        Jis = np.zeros(shape=(nj,), dtype = np.int)
        Jjs = np.zeros(shape=(nj,), dtype = np.int)

        r_b = np.zeros(shape=(3,), dtype = np.int)
        r_s = np.zeros(shape=(3,), dtype = np.int)
        Bx = np.zeros(shape=(3,3,), dtype = np.int)
        Sx = np.zeros(shape=(3,3,), dtype = np.int)
        cross = np.zeros(shape=(3,), dtype = np.int)
        #ddist_cross = np.zeros(shape=(3,), dtype = np.int)
        dcross_drb = np.zeros(shape=(3,3,), dtype = np.int)
        dcross_drs = np.zeros(shape=(3,3,), dtype = np.int)
        dLOS_dx = np.zeros(shape=(3,), dtype = np.int)
        dLOS_drs = np.zeros(shape=(3,), dtype = np.int)
        dLOS_drb = np.zeros(shape=(3,), dtype = np.int)

        for i in range(self.n):
            r_b = self.r_e2b_I[:3,i]
            r_s = self.r_e2s_I[:3,i]
            Bx = crossMatrix(r_b)
            Sx = crossMatrix(-r_s)
            dot = np.dot(r_b,r_s)
            cross = np.cross(r_b,r_s)
            dist = np.sqrt(np.dot(cross,cross))

            if dot >= 0.0 :
                dLOS_drb[:] = 0.0
                dLOS_drs[:] = 0.0
            elif dist <= self.r1 :
                dLOS_drb[:] = 0.0
                dLOS_drs[:] = 0.0
            elif dist >= self.r2 :
                dLOS_drb[:] = 0.0
                dLOS_drs[:] = 0.0
            else:
                x = (dist-self.r1)/(self.r2-self.r1)
                #LOS = 3*x**2 - 2*x**3
                ddist_dcross = cross/dist
                dcross_drb = Sx
                dcross_drs = Bx
                dx_ddist = 1.0/(self.r2-self.r1)
                dLOS_dx = 6*x - 6*x**2
                dLOS_drb = dLOS_dx*dx_ddist*np.dot(ddist_dcross,dcross_drb)
                dLOS_drs = dLOS_dx*dx_ddist*np.dot(ddist_dcross,dcross_drs)

            for k in range(3) :
                iJ = i*3 + k
                Jab[iJ] = dLOS_drb[k]
                Jib[iJ] = i
                Jjb[iJ] = (i)*6 + k
                Jas[iJ] = dLOS_drs[k]
                Jis[iJ] = i
                Jjs[iJ] = (i)*3 + k

        self.Jb = scipy.sparse.csc_matrix((Jab,(Jib,Jjb)),shape=(self.n,6*self.n))
        self.Js = scipy.sparse.csc_matrix((Jas,(Jis,Jjs)),shape=(self.n,3*self.n))
        self.JbT = self.Jb.transpose()
        self.JsT = self.Js.transpose()

    def apply_deriv(self, arg, result):

        if 'r_e2b_I' in arg:
            r_e2b_I = arg['r_e2b_I'][:].reshape((6*self.n),order='F')
            result['LOS'] += self.Jb.dot(r_e2b_I)

        if 'r_e2s_I' in arg:
            r_e2s_I = arg['r_e2s_I'][:].reshape((3*self.n),order='F')
            result['LOS'] += self.Js.dot(r_e2s_I)

    def apply_derivT(self, arg, result):

        if 'LOS' in arg:
            LOS = arg['LOS']

            if 'r_e2b_I' in result:
                result['r_e2b_I'] += self.JbT.dot(LOS).reshape((6,self.n),order='F')
            if 'r_e2s_I' in result:
                result['r_e2s_I'] += self.JsT.dot(LOS).reshape((3,self.n),order='F')

def crossMatrix(v):

        # so m[1,0] is v[2], for example
        m = np.array( [
                         [ 0.0, -v[2], v[1] ],
                         [ v[2],  0.0, -v[0]],
                         [ -v[1], v[0], 0.0 ]
                       ] )
        return m

class Sun_PositionBody( Component ):

    '''Position vector from earth to sun in body-fixed frame'''
    
    def __init__(self, n=2):
        super(Sun_PositionBody, self).__init__()

        self.n = n

        self.add('O_BI', Array(np.zeros((3, 3, n), order='F'),
                               size=(3,3,n, ), dtype=np.float,
                               units="unitless",
                               desc="Rotation matrix from the " +
                               "Earth-centered inertial frame " +
                               "to the satellite frame.",
                               iotype="in"))

        self.add('r_e2s_I', Array(np.zeros((3, n), order='F'),
                                  size=(3,n, ), dtype=np.float,
                                  units="km", desc="Position vector " +
                                  "from Earth to Sun in Earth-centered " +
                                  "inertial frame over time.",
                                  iotype="in"))

        self.add('r_e2s_B', Array(np.zeros((3,n, ), order='F'),
                                  size=(3,n, ), dtype=np.float,
                                  iotype="out",
                                  units = "km", desc="Position vector " +
                                  "from Earth to Sun in body-fixed " +
                                  "frame over time." ))

    def execute(self):
        self.r_e2s_B = computepositionrotd(self.n, self.r_e2s_I, self.O_BI)

    def linearize(self):
        self.J1, self.J2 = computepositionrotdjacobian(self.n, self.r_e2s_I, self.O_BI )

    def apply_deriv(self, arg, result):

        if 'O_BI' in arg:
            for k in range(3):
                for u in range(3):
                    for v in range(3):
                        result['r_e2s_B'][k,:] += self.J1[:,k,u,v] * arg['O_BI'][u,v,:]

        if 'r_e2s_I' in arg:
            for k in range(3):
                for j in range(3):
                    result['r_e2s_B'][k,:] += self.J2[:,k,j] * arg['r_e2s_I'][j,:]

    def apply_derivT(self, arg, result):

        if 'r_e2s_B' in arg:
            for k in range(3):
                
                if 'O_BI' in result:
                    for u in range(3):
                        for v in range(3):
                            result['O_BI'][u,v,:] += self.J1[:,k,u,v] * arg['r_e2s_B'][k,:]
                if 'r_e2s_I' in result:
                    for j in range(3):
                        result['r_e2s_I'][j,:] += self.J2[:,k,j] * arg['r_e2s_B'][k,:]


class Sun_PositionECI( Component ):

    '''
    Compute the position vector from Earth to Sun in
    Earth-centered inertial frame.
    '''

    #constants
    d2r = np.pi/180.
    LD = Float(0., iotype="in", units="unitless", copy=None)

    def __init__(self, n=2):
        super(Sun_PositionECI, self).__init__()

        self.n = n

        #self.add('LD', Array(np.zeros((1,), order='F'), size=(1,), dtype=np.float, iotype="in"))

        self.add(
		    't', 
		    Array(
		        np.zeros((n,), order='F'),
			size=(n,),
			dtype=np.float,
                        units="s",
			desc="Time",
			iotype="in"
		    )
	)

        self.add('r_e2s_I', Array(np.zeros((3,n, ), order='F'),
                                  size=(3,n, ),
                                  dtype=np.float,
                                  units="km",
                                  desc="Position vector from Earth " +
                                  "to Sun in Earth-centered inertial " +
                                  "frame over time.",
                                  iotype="out"))

        self.Ja = np.zeros(3*self.n)
        self.Ji = np.zeros(3*self.n)
        self.Jj = np.zeros(3*self.n)

    def execute(self):
        T = self.LD + self.t[:]/3600./24.
        for i in range(0,self.n):
            L = self.d2r*280.460 + self.d2r*0.9856474*T[i]
            g = self.d2r*357.528 + self.d2r*0.9856003*T[i]
            Lambda = L + self.d2r*1.914666*np.sin(g) + self.d2r*0.01999464*np.sin(2*g)
            eps = self.d2r*23.439 - self.d2r*3.56e-7*T[i]
            self.r_e2s_I[0,i] = np.cos(Lambda)
            self.r_e2s_I[1,i] = np.sin(Lambda)*np.cos(eps)
            self.r_e2s_I[2,i] = np.sin(Lambda)*np.sin(eps)

    def linearize(self):
        T = self.LD + self.t[:]/3600./24.
        dr_dt = np.empty(3)
        for i in range(0,self.n):
            L = self.d2r*280.460 + self.d2r*0.9856474*T[i]
            g = self.d2r*357.528 + self.d2r*0.9856003*T[i]
            Lambda = L + self.d2r*1.914666*np.sin(g) + self.d2r*0.01999464*np.sin(2*g)
            eps = self.d2r*23.439 - self.d2r*3.56e-7*T[i]

            dL_dt = self.d2r*0.9856474
            dg_dt = self.d2r*0.9856003
            dlambda_dt = dL_dt + self.d2r*1.914666*np.cos(g)*dg_dt + self.d2r*0.01999464*np.cos(2*g)*2*dg_dt
            deps_dt = -self.d2r*3.56e-7

            dr_dt[0] = -np.sin(Lambda)*dlambda_dt
            dr_dt[1] = np.cos(Lambda)*np.cos(eps)*dlambda_dt - np.sin(Lambda)*np.sin(eps)*deps_dt
            dr_dt[2] = np.cos(Lambda)*np.sin(eps)*dlambda_dt + np.sin(Lambda)*np.cos(eps)*deps_dt

            for k in range(0,3):
                #iJ = (i-1)*3 + k #This is the original implementation, but it may not work because of index differences between fortran and python
                #self.Ja[iJ] = dr_dt[k]
                #self.Ji[iJ] = iJ - 1
                #self.Jj[iJ] = i - 1
                iJ = i*3 + k #This should resolve the index issues
                self.Ja[iJ] = dr_dt[k]
                self.Ji[iJ] = iJ
                self.Jj[iJ] = i

        self.J = scipy.sparse.csc_matrix((self.Ja,(self.Ji,self.Jj)),shape=(3*self.n,self.n))
        self.JT = self.J.transpose()

    def apply_deriv(self, arg, result):

        if 'LD' in arg and 't' in arg:
            result['r_e2s_I'][:] += self.J.dot( arg['LD'] + arg['t']/3600./24. ).reshape((3,self.n),order='F')

    def apply_derivT(self, arg, result):

        if 'r_e2s_I' in arg:
            r_e2s_I = arg['r_e2s_I'][:].reshape((3*self.n),order='F')
            if 'LD' in result:
                result['LD'] += sum(self.JT.dot(r_e2s_I))
            if 't' in result:
                result['t'][:] += self.JT.dot(r_e2s_I)/3600.0/24.0


class Sun_PositionSpherical(Component):

    '''Compute the elevation angle of the Sun in the body-fixed frame.'''

    def __init__(self, n=2):
        super(Sun_PositionSpherical, self).__init__()

        self.n = n

        self.add('r_e2s_B', Array(np.zeros((3, n)), size=(3, n),
                                  units = "km",
                                  desc="Position vector from " +
                                      "Earth to Sun in body-fixed " +
                                      "frame over time.",
                                  dtype=np.float, iotype="in"))

        self.add('azimuth', Array(np.zeros((n,)), size=(n,), dtype=np.float,
                                  units='rad',
                                  desc='Ezimuth angle of the Sun ' +
                                      'in the body-fixed frame over time.',
                                  iotype="out"))
        self.add('elevation', Array(np.zeros((n,)), size=(n,), dtype=np.float,
                                    units='rad',
                                    desc="Elevation angle of the " +
                                        "Sun in the body-fixed frame " +
                                        "over time.",
                                    iotype="out"))

    def execute(self):
        azimuth, elevation = computepositionspherical(self.n, self.r_e2s_B)

        self.azimuth = azimuth
        self.elevation = elevation

    def linearize(self):
        self.Ja1, self.Ji1, self.Jj1, self.Ja2, self.Ji2, self.Jj2 = \
                  computepositionsphericaljacobian(self.n, 3*self.n, self.r_e2s_B)
        self.J1 = scipy.sparse.csc_matrix((self.Ja1, (self.Ji1, self.Jj1)),
                                          shape=(self.n, 3*self.n))
        self.J2 = scipy.sparse.csc_matrix((self.Ja2, (self.Ji2, self.Jj2)),
                                          shape=(self.n, 3*self.n))
        self.J1T = self.J1.transpose()
        self.J2T = self.J2.transpose()

        return

    def apply_deriv(self, arg, result):

        if 'r_e2s_B' in arg:
            r_e2s_B = arg['r_e2s_B'].reshape((3*self.n), order='F')
            result['azimuth'] += self.J1.dot(r_e2s_B)
            result['elevation'] += self.J2.dot(r_e2s_B)


    def apply_derivT(self, arg, result):

        if 'azimuth' in arg and 'r_e2s_B' in result:
            azimuth = arg['azimuth'][:]
            result['r_e2s_B'][:] += self.J1T.dot(azimuth).reshape((3,self.n), order='F')

        if 'elevation' in arg and 'r_e2s_B' in result:
            elevation = arg['elevation'][:]
            result['r_e2s_B'][:] += self.J2T.dot(elevation).reshape((3,self.n), order='F')

        return result
