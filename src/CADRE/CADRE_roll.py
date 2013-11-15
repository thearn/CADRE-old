from openmdao.main.api import Assembly, Component
from openmdao.main.datatypes.api import Float, Array, Int
import numpy as np
from CADRE import CADRE
from pyopt_driver import pyopt_driver
from openmdao.lib.drivers.api import CONMINdriver
import os


class Deflection(Component):

    def __init__(self, n=1500):
        super(Deflection, self).__init__()
        self.add('r_b2g_A', Array(np.zeros((3, n)), iotype='in',
                                  shape=(3, n)))
        self.add("d", Array(np.zeros(n), iotype='in'))

    def linearize(self):
        self.J = np.zeros((3, n))

        self.J[0] = self.r_b2g_A[0] / d
        self.J[1] = self.r_b2g_A[1] / d
        self.J[2] = self.r_b2g_A[2] / d

    def apply_deriv(self, arg, result):
        if 'r_b2g_A' in arg:
            result['d'] += np.sum(self.J * arg['r_b2g_A'])

    def apply_derivT(self, arg, result):
        if 'd' in arg:
            result['r_b2g_A'][0, :] += arg['d'] * self.J[0]
            result['r_b2g_A'][1, :] += arg['d'] * self.J[1]
            result['r_b2g_A'][2, :] += arg['d'] * self.J[2]

    def execute(self):
        self.d = np.sqrt(np.sum(self.r_b2g_A ** 2, axis=0))


class CADRE_Roll(Assembly):

    def __init__(self, n=1500, m=300):
        super(CADRE_Roll, self).__init__()

        # add SNOPT driver
        self.add("driver", pyopt_driver.pyOptDriver())
        self.driver.optimizer = "SNOPT"
        self.driver.options = {'Major optimality tolerance': 1e-8,
                               'Iterations limit': 500000000}

        #self.add("driver", CONMINdriver())

        # orbit position and velocity data for each design point
        r_e2b_I0 = [-4969.91222,  4624.84149,
                    1135.9414,  0.1874654, -1.62801666,  7.4302362]

        # number of days since launch for each design point
        LD = 5417.5

        self.add("CADRE", CADRE(n, m))
        self.CADRE.set("LD", LD)
        self.CADRE.set("r_e2b_I0", r_e2b_I0)

        self.add("Deflection", Deflection(n))
        self.driver.workflow.add("Deflection")

        self.CADRE.create_passthrough("Comm_VectorAnt.r_b2g_A")
        self.connect("CADRE.r_b2g_A", "Deflection.r_b2g_A")

        # add parameters to driver

        self.driver.add_parameter("CADRE.CP_gamma", low=0,
                                  high=np.pi / 2.)

        # add objective
        self.driver.add_objective("sum(Deflection.d)")

if __name__ == "__main__":
    a = CADRE_Roll()
    a.run()
