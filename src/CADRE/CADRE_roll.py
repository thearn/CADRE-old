from openmdao.main.api import Assembly, Component
from openmdao.main.datatypes.api import Float, Array, Int
import numpy as np
from CADRE import CADRE
from pyopt_driver import pyopt_driver
from openmdao.lib.drivers.api import CONMINdriver
from openmdao.lib.drivers.api import SLSQPdriver
import os


class Deflection(Component):

    def __init__(self, n=1500):
        super(Deflection, self).__init__()
        self.add('r_b2g_A', Array(np.zeros((3, n)), iotype='in',
                                  shape=(3, n)))
        #self.add("d", Array(np.zeros(n), iotype='out'))
        self.add("d", Float(1., iotype='out'))
        self.n = n
        self.nd = np.zeros(n)

    def linearize(self):
        self.J = np.zeros((3, self.n))

        self.J[0] = self.r_b2g_A[0] / self.nd
        self.J[1] = self.r_b2g_A[1] / self.nd
        self.J[2] = self.r_b2g_A[2] / self.nd

    def apply_deriv(self, arg, result):
        if 'r_b2g_A' in arg:
            result['d'] += np.sum(self.J * arg['r_b2g_A'])

    def apply_derivT(self, arg, result):
        if 'd' in arg:
            result['r_b2g_A'][0, :] += arg['d'] * self.J[0]
            result['r_b2g_A'][1, :] += arg['d'] * self.J[1]
            result['r_b2g_A'][2, :] += arg['d'] * self.J[2]

    def execute(self):
        self.nd = np.sqrt(np.sum(self.r_b2g_A ** 2, axis=0))
        self.d = sum(self.nd)


class CADRE_Roll(Assembly):

    def __init__(self, n=1500, m=300):
        super(CADRE_Roll, self).__init__()

        # add SNOPT driver
        #self.add("driver", pyopt_driver.pyOptDriver())
        #self.driver.optimizer = "SNOPT"
        # self.driver.options = {'Major optimality tolerance': 1e-8,
        #                       'Iterations limit': 500000000}

        #self.add("driver", SLSQPdriver())

        # orbit position and velocity data for each design point
        r_e2b_I0 = [-4969.91222,  4624.84149,
                    1135.9414,  0.1874654, -1.62801666,  7.4302362]

        # number of days since launch for each design point
        LD = 5417.5

        self.add("CADRE", CADRE(n, m))
        self.driver.workflow.add("CADRE")
        self.CADRE.set("LD", LD)
        self.CADRE.set("r_e2b_I0", r_e2b_I0)

        self.add("Deflection", Deflection(n))
        self.driver.workflow.add("Deflection")

        self.CADRE.create_passthrough("Comm_VectorBody.r_b2g_B")
        self.connect("CADRE.r_b2g_B", "Deflection.r_b2g_A")

        # add parameters to driver

        # self.driver.add_parameter("CADRE.CP_gamma", low=0,
        #                          high=np.pi / 2.)

        # add objective
        # self.driver.add_objective("Deflection.d")

if __name__ == "__main__":
    import pylab
    import pickle
    a = CADRE_Roll()
    a.run()
    pylab.figure()
    pylab.subplot(211)
    pylab.plot(a.Deflection.nd)
    # pylab.plot(a.CADRE.BsplineParameters.CP_gamma)
    print a.Deflection.d

    a.CADRE.CP_gamma = pickle.load(
        open("src/CADRE/test/data1346.pkl"))["3:CP_gamma"]
    a.run()
    print a.Deflection.d
    pylab.subplot(212)
    pylab.plot(a.Deflection.nd)
    # pylab.plot(a.CADRE.BsplineParameters.CP_gamma)
    pylab.show()
