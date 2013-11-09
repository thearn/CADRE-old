import numpy as np

from openmdao.main.api import Assembly, Component

from openmdao.lib.drivers.api import CONMINdriver
from openmdao.lib.drivers.api import SLSQPdriver
from pyopt_driver import pyopt_driver

from openmdao.main.datatypes.api import Float, Array, Int

from comm import Comm_EarthsSpin, Comm_EarthsSpinMtx
from orbit import Orbit_Initial, Orbit_Dynamics
from scipy.stats import kurtosis
import sys


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
        self.k = (kurtosis(self.sample.flatten()) + 6. / 5) ** 2
        sys.stdout.flush()


class GroundLOC(Component):

    """ Gives the lat and lon location of a satellite
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

    def execute(self):
        for i in xrange(self.n):
            self.pos = self.r_e2b_I[:3, i]
            self.npos = self.pos / np.linalg.norm(self.pos) * self.Re
            self.g_pos = np.dot(self.O_IE[:, :, i].T, self.npos)
            self.lats[i] = np.arcsin(self.g_pos[2] / self.Re) * self.r2d
            self.lons[i] = np.arctan2(self.g_pos[1], self.g_pos[0]) * self.r2d


class CADRE_Launch(Assembly):

    """ Optimization of the launch parameters of CADRE
    """

    def __init__(self, n=500):

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

        #self.add('driver', CONMINdriver())

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

        #self.driver.add_objective("Lat_uniform.k + Lon_uniform.k")
        # self.driver.add_parameter(
        #     "Orbit_Initial.altPerigee", low=500, high=1000)
        # self.driver.add_parameter(
        #     "Orbit_Initial.altApogee", low=500, high=1000)
        # self.driver.add_parameter(
        #     "Orbit_Initial.RAAN", low=-180, high=180)
        # self.driver.add_parameter(
        #     "Orbit_Initial.Inc", low=0, high=90)
        # self.driver.add_parameter(
        #     "Orbit_Initial.argPerigee", low=0, high=90)

if __name__ == "__main__":
    import time
    from scipy.optimize import fmin, fmin_slsqp
    a = CADRE_Launch()
    a.add('driver', SLSQPdriver())
    a.driver.add_objective("Lat_uniform.k + Lon_uniform.k")
    a.driver.add_parameter(
        "Orbit_Initial.altPerigee", low=500, high=1000)
    a.driver.add_parameter(
        "Orbit_Initial.altApogee", low=500, high=1000)
    a.driver.add_parameter(
        "Orbit_Initial.RAAN", low=-180, high=180)
    a.driver.add_parameter(
        "Orbit_Initial.Inc", low=0, high=90)
    a.driver.add_parameter(
        "Orbit_Initial.argPerigee", low=0, high=90)
    tt = time.time()
    print 30 * "-"
    print "with OpenMDAO optimizer:"
    l1, l2 = a.GroundLOC.lats, a.GroundLOC.lons
    print "min/max lats:", min(l1), max(l1)
    print "min/max lons:", min(l2), max(l2)
    print "\n"
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
    c1 = lambda x: x[0] - 500
    c2 = lambda x: 1000 - x[0]
    c3 = lambda x: x[1] - 500
    c4 = lambda x: 1000 - x[1]
    c5 = lambda x: x[2] + 180
    c6 = lambda x: 180 - x[2]
    c7 = lambda x: x[3]
    c8 = lambda x: 90 - x[3]

    def f(orbit):
        a.Orbit_Initial.altPerigee = orbit[0]
        a.Orbit_Initial.altApogee = orbit[1]
        a.Orbit_Initial.RAAN = orbit[2]
        a.Orbit_Initial.Inc = orbit[3]
        a.run()
        return a.Lat_uniform.k + a.Lon_uniform.k
    fmin_slsqp(f, [600, 600, 0, 45], ieqcons=[c1, c2, c3, c4, c5, c6, c7, c8])
    l1, l2 = a.GroundLOC.lats, a.GroundLOC.lons
    print "min/max lats:", min(l1), max(l1)
    print "min/max lons:", min(l2), max(l2)
    print "\n"
    print(a.Orbit_Initial.altPerigee,
          a.Orbit_Initial.altApogee,
          a.Orbit_Initial.RAAN,
          a.Orbit_Initial.Inc,
          a.Orbit_Initial.argPerigee)
    print "Elapsed time: ", time.time() - tt, "seconds"
