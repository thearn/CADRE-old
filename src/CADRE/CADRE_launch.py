import numpy as np

from openmdao.main.api import Assembly, Component

from openmdao.lib.drivers.api import CONMINdriver
from openmdao.lib.drivers.api import SLSQPdriver
from pyopt_driver import pyopt_driver

from openmdao.main.datatypes.api import Float, Array, Int

from comm import Comm_EarthsSpin, Comm_EarthsSpinMtx
from orbit import Orbit_Initial, Orbit_Dynamics
from scipy.stats import kurtosis, moment


class Uniformity(Component):

    """ Measures the departure of a sample from a uniform distribution
        based on excess kurtosis (4th order sample moment).
    """

    def __init__(self, n):
        super(Uniformity, self).__init__()
        self.n = n
        self.add('sample', Array(iotype='in'))
        self.add('k', Float(0., iotype='out'))

    def execute(self):
        self.k_ = kurtosis(self.sample.flatten())
        self.k = (self.k_ + 6. / 5) ** 2

    def linearize(self):
        self.J = 2 * (self.k_ + 6. / 5) * self.d_kurtosis(self.sample)

    def apply_deriv(self, arg, result):
        if "sample" in arg and "k" in result:
            result['k'] += np.dot(self.J, arg['sample'])

    def apply_derivT(self, arg, result):
        if "sample" in result and "k" in arg:
            result['sample'] += arg['k'] * self.J

    def d_kurtosis(self, x):
        return (moment(x, 2) * self.d_moment(x, 4) - 2 * moment(x, 4)
                * self.d_moment(x, 2)) / moment(x, 2) ** 3

    def d_moment(self, x, p=2):
        n = float(len(x))
        return self.d_norm_p(x - np.mean(x), p) / n

    def d_norm_p(self, x, p=2):
        d = p * np.linalg.norm(x, p) ** (p - 1)
        d = d * self.d_norm(x, p)
        return d

    def d_norm(self, x, p=2):
        return np.abs(x) ** (p - 2) * (x) / (np.linalg.norm(x, p) ** (p - 1))


class GroundLOC(Component):

    """ Gives the lat and lon location of the ground beneath a satellite
    """
    Re = 6378.137
    r2d = 180 / np.pi

    def __init__(self, n):
        super(GroundLOC, self).__init__()
        self.n = n
        self.add('O_IE', Array(np.zeros((3, 3, self.n)), iotype='in'))

        self.add('r_e2b_I', Array(np.zeros((6, self.n)), iotype='in'))

        self.add('lats', Array(np.zeros(self.n), iotype='out'))
        self.add('lons', Array(np.zeros(self.n), iotype='out'))

    def linearize(self):
        self.J = np.zeros((self.n, 2, 3))
        self.J_O_IE = np.zeros((self.n, 2, 3, 3))

        for i in xrange(self.n):

            O = self.O_IE[:,:, i].T
            r = self.r_e2b_I[:3, i]
            d = self.r2d

            self.J[i, 0, 0] = self.d_lat(
                O[2, 0], O[2, 1], O[2, 2], d, r[0], r[1], r[2])
            self.J[i, 0, 1] = self.d_lat(
                O[2, 1], O[2, 0], O[2, 2], d, r[1], r[0], r[2])
            self.J[i, 0, 2] = self.d_lat(
                O[2, 2], O[2, 1], O[2, 0], d, r[2], r[1], r[0])

            self.J[i, 1, 0] = self.d_lon(O[1, 0], O[1, 1], O[1, 2], d,
                                         O[0, 0], O[0, 1], O[0, 2], r[0], r[1], r[2])
            self.J[i, 1, 1] = self.d_lon(O[1, 1], O[1, 0], O[1, 2], d,
                                         O[0, 1], O[0, 0], O[0, 2], r[1], r[0], r[2])
            self.J[i, 1, 2] = self.d_lon(O[1, 2], O[1, 1], O[1, 0], d,
                                         O[0, 2], O[0, 1], O[0, 0], r[2], r[1], r[0])

            # self.J_O_IE[i, 0, 2, 0] = self.d_lat_O_IE(
            #     r[0], r[1], r[2], d, O[2, 0], O[2, 1], O[2, 2])
            # self.J_O_IE[i, 0, 2, 1] = self.d_lat_O_IE(
            #     r[1], r[0], r[2], d, O[2, 1], O[2, 0], O[2, 2])
            # self.J_O_IE[i, 0, 2, 2] = self.d_lat_O_IE(
            #     r[2], r[1], r[0], d, O[2, 2], O[2, 1], O[2, 0])

            # self.J_O_IE[i, 1, 0, 0] = self.d_lon_O_IE(r[0], r[1], r[2], d,
            #                                           O[1, 0], O[1, 1],
            #                                           O[1, 2], O[0, 0],
            #                                           O[0, 1], O[0, 2])

            # self.J_O_IE[i, 1, 0, 1] = self.d_lon_O_IE(r[1], r[0], r[2], d,
            #                                           O[1, 0], O[1, 1],
            #                                           O[1, 2], O[0, 1],
            #                                           O[1, 0], O[0, 2])

    def apply_deriv(self, arg, result):
        if 'r_e2b_I' in arg:
            result['lats'] += np.dot(self.J[:, 0,:], arg['r_e2b_I'][:3])
            result['lons'] += np.dot(self.J[:, 1,:], arg['r_e2b_I'][:3])

    def apply_derivT(self, arg, result):

        if 'lats' in arg:
            result['r_e2b_I'][0,:] += arg['lats'] * self.J[:, 0, 0]
            result['r_e2b_I'][1,:] += arg['lats'] * self.J[:, 0, 1]
            result['r_e2b_I'][2,:] += arg['lats'] * self.J[:, 0, 2]
        if 'lons' in arg:
            result['r_e2b_I'][0,:] += arg['lons'] * self.J[:, 1, 0]
            result['r_e2b_I'][1,:] += arg['lons'] * self.J[:, 1, 1]
            result['r_e2b_I'][2,:] += arg['lons'] * self.J[:, 1, 2]

    def d_lat(self, a, b, c, d, x, y, z):
        return -(d * (-a * (y ** 2 + z ** 2) + b * x * y + c * x * z)) / (np.sqrt((x ** 2 + y ** 2 + z ** 2) ** 3) *
                                                                          np.sqrt(-((a ** 2 - 1) * x ** 2 + 2 * a * x * (b * y + c * z) + (b ** 2 - 1) * y ** 2 + 2 * b * c * y * z + (c ** 2 - 1) * z ** 2) / (x ** 2 + y ** 2 + z ** 2)))

    def d_lon(self, a, b, c, d, e, f, g, x, y, z):
        num = -(d * (-a * (f * y + g * z) + b * e * y + c * e * z))
        denom = (a ** 2 * x ** 2 + 2 * a * x * (b * y + c * z) + b ** 2 * y ** 2 + 2 * b * c * y * z + c ** 2 *
                 z ** 2 + e ** 2 * x ** 2 + 2 * e * x * (f * y + g * z) + f ** 2 * y ** 2 + 2 * f * g * y * z + g ** 2 * z ** 2)
        if not denom:
            return 0.
        else:
            return num / denom

    def d_lat_O_IE(self, a, b, c, d, x, y, z):
        return d * a / (np.sqrt(a ** 2 + b ** 2 + c ** 2) * np.sqrt(1 - (a * x + b * y + c * z) ** 2 / (a ** 2 + b ** 2 + c ** 2)))

    def d_lon_O_IE(self, a, b, c, d, w, v, u, x, y, z,):
        return -d * (a * (a * w + b * v + c * u)) / (a ** 2 * (w ** 2 + x ** 2) + 2 * c * (a * u * w + a * x * z + b * u * v + b * y * z) + 2 * a * b * (v * w + x * y) + b ** 2 * (v ** 2 + y ** 2) + c ** 2 * (u ** 2 + z ** 2))

    def execute(self):
        self.npos = np.zeros((self.n, 3))
        self.g_pos = np.zeros((self.n, 3))
        for i in xrange(self.n):
            self.pos = self.r_e2b_I[:3, i]
            self.npos[i,:] = self.pos / np.linalg.norm(self.pos, 2) * self.Re
            self.g_pos[i,:] = np.dot(self.O_IE[:,:, i].T, self.npos[i,:])
            self.lats[i] = np.arcsin(self.g_pos[i, 2] / self.Re) * self.r2d
            self.lons[i] = np.arctan2(
                self.g_pos[i, 1], self.g_pos[i, 0]) * self.r2d


class CADRE_Launch(Assembly):

    """ Allows for analysis of the launch parameters of CADRE.

        Considers multiple launch parameters and their effects on
        coverage of the Earth's thermosphere.

        Ultimately, a launch that provides the most uniform sampling is
        favorable, which is expected to be given by a polar orbit
        (Inclination near 90). Other launch parameters are probably not
        very useful in comparison.
    """

    def __init__(self, n=200):

        super(CADRE_Launch, self).__init__()

        # Analysis parameters
        self.n = n
        self.add('t', Array(np.zeros((n,)), size=(n,),
                            dtype=np.float, iotype="in"))
        self.add('t1', Float(0., iotype='in'))
        self.add('t2', Float(43200., iotype='in'))
        h = (self.t2 - self.t1) / (self.n - 1)
        self.add("h", Float(h, iotype="in", copy=None))

        self.t = np.array(range(0, n)) * h

        #self.add("driver", pyopt_driver.pyOptDriver())
        #self.driver.optimizer = "SNOPT"

        #self.add('driver', SLSQPdriver())

        self.add('driver', CONMINdriver())

        # Orbit components
        self.add("Orbit_Initial", Orbit_Initial())
        self.driver.workflow.add("Orbit_Initial")

        self.add("Orbit_Dynamics", Orbit_Dynamics(n))
        self.driver.workflow.add("Orbit_Dynamics")

        self.add("Comm_EarthsSpin", Comm_EarthsSpin(n))
        self.driver.workflow.add("Comm_EarthsSpin")

        self.add("Comm_EarthsSpinMtx", Comm_EarthsSpinMtx(n))
        self.driver.workflow.add("Comm_EarthsSpinMtx")

        self.add("GroundLOC", GroundLOC(n))
        self.driver.workflow.add("GroundLOC")

        self.add("Lon_uniform", Uniformity(n))
        self.driver.workflow.add("Lon_uniform")

        self.add("Lat_uniform", Uniformity(n))
        self.driver.workflow.add("Lat_uniform")

        self.connect("t", "Comm_EarthsSpin.t")
        self.connect("h", "Orbit_Dynamics.h")
        self.connect("Comm_EarthsSpin.q_E", "Comm_EarthsSpinMtx.q_E")
        self.connect("Comm_EarthsSpinMtx.O_IE", "GroundLOC.O_IE")

        self.connect("Orbit_Initial.r_e2b_I0", "Orbit_Dynamics.r_e2b_I0")
        self.connect("Orbit_Dynamics.r_e2b_I", "GroundLOC.r_e2b_I")

        self.connect("GroundLOC.lats", "Lat_uniform.sample")
        self.connect("GroundLOC.lons", "Lon_uniform.sample")

        self.driver.add_objective("Lat_uniform.k + Lon_uniform.k")
        self.driver.add_parameter(
            "Orbit_Initial.altPerigee", low=500, high=1000)
        self.driver.add_parameter(
            "Orbit_Initial.altApogee", low=500, high=1000)
        self.driver.add_parameter(
            "Orbit_Initial.RAAN", low=-180, high=180)
        self.driver.add_parameter(
            "Orbit_Initial.Inc", low=0, high=90)
        self.driver.add_parameter(
            "Orbit_Initial.argPerigee", low=0, high=90)

if __name__ == "__main__":
    import time
    from scipy.optimize import fmin, fmin_slsqp
    print 30 * "-"
    print "with OpenMDAO optimizer:"
    a = CADRE_Launch()
    #a.add('driver', SLSQPdriver())
    #a.driver.add_objective("Lat_uniform.k + Lon_uniform.k")
    #a.driver.add_parameter(
    #    ["Orbit_Initial.altPerigee", "Orbit_Initial.altApogee"],
    #    low=500, high=1000)
    #a.driver.add_parameter(
    #    "Orbit_Initial.RAAN", low=-180, high=180)
    #a.driver.add_parameter(
    #    "Orbit_Initial.Inc", low=0, high=90)
    #a.driver.add_parameter(
    #    "Orbit_Initial.argPerigee", low=0, high=90)
    tt = time.time()
    a.run()
    a.driver.workflow.check_gradient(outputs=['Orbit_Dynamics.r_e2b_I'])
    exit()
    l1, l2 = a.GroundLOC.lats, a.GroundLOC.lons
    print "min/max lats:", min(l1), max(l1)
    print "min/max lons:", min(l2), max(l2)
    print "objective:", a.Lat_uniform.k + a.Lon_uniform.k
    print(a.Orbit_Initial.altPerigee,
          a.Orbit_Initial.altApogee,
          a.Orbit_Initial.RAAN,
          a.Orbit_Initial.Inc,
          a.Orbit_Initial.argPerigee)
    print "Elapsed time: ", time.time() - tt, "seconds"
    print 30 * "-"
    print "without OpenMDAO optimizer:"
    a = CADRE_Launch()
    tt = time.time()

    def f(orbit):
        a.Orbit_Initial.altPerigee = orbit[0]
        a.Orbit_Initial.altApogee = orbit[0]
        a.Orbit_Initial.RAAN = orbit[1]
        a.Orbit_Initial.Inc = orbit[2]
        a.run()
        return a.Lat_uniform.k + a.Lon_uniform.k
    fmin_slsqp(f, [600, 0, 45], bounds=[
               (500, 1000), (-180, 180), (0, 90)],
               iprint=1)
    l1, l2 = a.GroundLOC.lats, a.GroundLOC.lons
    print "min/max lats:", min(l1), max(l1)
    print "min/max lons:", min(l2), max(l2)
    print "objective:", a.Lat_uniform.k + a.Lon_uniform.k
    print(a.Orbit_Initial.altPerigee,
          a.Orbit_Initial.altApogee,
          a.Orbit_Initial.RAAN,
          a.Orbit_Initial.Inc,
          a.Orbit_Initial.argPerigee)
    print "Elapsed time: ", time.time() - tt, "seconds"
