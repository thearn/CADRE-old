from openmdao.main.api import Assembly, Component
from openmdao.main.datatypes.api import Float, Array, Int
import numpy as np
from CADRE import CADRE
from pyopt_driver import pyopt_driver
from openmdao.lib.drivers.api import CONMINdriver
from openmdao.lib.drivers.api import SLSQPdriver
import os


class NetGain(Component):

    gain = Array(iotype="in")
    net = Float(iotype="out")

    def execute(self):
        self.net = sum(self.gain)


class CADRE_Roll(Assembly):

    def __init__(self, n=1500, m=300):
        super(CADRE_Roll, self).__init__()

        # add SNOPT driver
        self.add("driver", pyopt_driver.pyOptDriver())
        self.driver.optimizer = "SNOPT"
        self.driver.options = {'Major optimality tolerance': 1e-2,
                               'Iterations limit': 5000}

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
        self.CADRE.create_passthrough("Comm_GainPattern.gain")
        self.create_passthrough("CADRE.gain")

        self.add("NetGain", NetGain())
        self.driver.workflow.add("NetGain")
        self.connect("gain", "NetGain.gain")

        self.driver.add_parameter("CADRE.CP_gamma", low=0,
                                  high=np.pi / 2.)

        # add objective
        self.driver.add_objective("-NetGain.net")

if __name__ == "__main__":
    import pylab
    import pickle
    a = CADRE_Roll()
    print sum(a.gain)
    a.run()
    print sum(a.gain)
    quit()
    # a.run()
    # pylab.figure()
    # pylab.subplot(211)
    # pylab.plot(a.Deflection.nd)
    # pylab.plot(a.CADRE.Comm_GainPattern.gain)

    # a.CADRE.CP_gamma = pickle.load(
    #     open("src/CADRE/test/data1346.pkl"))["3:CP_gamma"]
    # a.run()
    # pylab.subplot(212)
    # pylab.plot(a.Deflection.nd)
    # pylab.plot(a.CADRE.Comm_GainPattern.gain)
    # pylab.show()
