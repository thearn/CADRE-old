import os.path
import numpy as np

from openmdao.main.api import Assembly
from openmdao.lib.drivers.api import CONMINdriver

from pyopt_driver import pyopt_driver

from .CADRE_assembly import CADRE


class CADRE_Optimization(Assembly):

    def __init__(self, n=1500, m=300, npts=6):
        super(CADRE_Optimization, self).__init__()

        # add SNOPT driver
        driver = self.add("driver", pyopt_driver.pyOptDriver())
        driver.optimizer = "SNOPT"
        driver.options = {'Major optimality tolerance': 1e-4,
                          'Iterations limit': 500000000,
                          'New basis file': 10}
        if os.path.exists("fort.10"):
            driver.options["Old basis file"] = 10

        #driver = self.add("driver", CONMINdriver())

        # Raw data to load
        fpath = os.path.dirname(os.path.realpath(__file__))
        fpath = os.path.join(fpath, 'data')
        solar_raw1 = np.genfromtxt(fpath + '/Solar/Area10.txt')
        solar_raw2 = np.loadtxt(fpath + '/Solar/Area_all.txt')
        comm_rawGdata = np.genfromtxt(fpath + '/Comm/Gain.txt')
        comm_raw = (10 ** (comm_rawGdata / 10.0)).reshape((361, 361), order='F')
        power_raw = np.genfromtxt(fpath + '/Power/curve.dat')

        # Load launch data
        launch_data = np.loadtxt(fpath + '/Launch/launch1.dat')

        # orbit position and velocity data for each design point
        r_e2b_I0s = launch_data[1::2, 1:]

        # number of days since launch for each design point
        LDs = launch_data[1::2, 0] - 2451545

        # LDs = [5233.5, 5294.5, 5356.5, 5417.5, 5478.5, 5537.5]

        # r_e2b_I0s = [np.array([4505.29362, -3402.16069, -3943.74582,
        #                        4.1923899, -1.56280012,  6.14347427]),
        #              np.array(
        #                  [-1005.46693,  -597.205348, -6772.86532, -0.61047858,
        #                   -7.54623146,  0.75907455]),
        #              np.array(
        #                  [4401.10539,  2275.95053, -4784.13188, -5.26605537,
        #                   -1.08194926, -5.37013745]),
        #              np.array(
        #                  [-4969.91222,  4624.84149,  1135.9414,  0.1874654,
        #                   -1.62801666,  7.4302362]),
        #              np.array(
        #                  [-235.021232,  2195.72976,  6499.79919, -2.55956031,
        #                   -6.82743519,  2.21628099]),
        #              np.array(
        #                  [-690.314375, -1081.78239, -6762.90367,  7.44316722,
        #                   1.19745345, -0.96035904])]

        # build design points
        names = ['pt%s' % i for i in range(npts)]
        for i, name in enumerate(names):
            comp = self.add(name, CADRE(n, m, solar_raw1, solar_raw2,
                                        comm_raw, power_raw))
            comp.set("LD", LDs[i])
            comp.set("r_e2b_I0", r_e2b_I0s[i])

            # add parameters to driver
            driver.add_parameter("%s.CP_Isetpt" % name, low=0., high=0.4)
            driver.add_parameter("%s.CP_gamma" % name, low=0, high=np.pi / 2.)
            driver.add_parameter("%s.CP_P_comm" % name, low=0., high=25.)
            driver.add_parameter("%s.iSOC[0]" % name, low=0.2, high=1.)

            # add constraints
            driver.add_constraint("%s.ConCh <= 0" % name)
            driver.add_constraint("%s.ConDs <= 0" % name)
            driver.add_constraint("%s.ConS0 <= 0" % name)
            driver.add_constraint("%s.ConS1 <= 0" % name)
            driver.add_constraint("%s.SOC[0][0] = %s.SOC[0][-1]" % (name, name))

        # add parameter groups
        cell_param = ["%s.cellInstd" % name for name in names]
        driver.add_parameter(cell_param, low=0, high=1)

        finangles = ["%s.finAngle" % name for name in names]
        driver.add_parameter(finangles, low=0, high=np.pi / 2.)

        antangles = ["%s.antAngle" % name for name in names]
        driver.add_parameter(antangles, low=-np.pi / 4, high=np.pi / 4)

        # add objective
        obj = ''.join(["-%s.Data[0][-1]" % name for name in names])
        driver.add_objective(obj)
