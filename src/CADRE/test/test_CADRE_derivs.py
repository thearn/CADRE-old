
import unittest
import numpy as np

from openmdao.main.api import Assembly, set_as_top, Driver
from openmdao.util.testutil import assert_rel_error
from CADRE.CADRE_assembly import CADRE

import os

class Testcase_CADRE_deriv(unittest.TestCase):

    """ Test run/step/stop aspects of a simple workflow. """

    def test_RK4_issue(self):

        n = 60
        m = 20
        
        solar_raw1 = np.genfromtxt('../data/Solar/Area10.txt')
        solar_raw2 = np.loadtxt("../data/Solar/Area_all.txt")
        comm_rawGdata = np.genfromtxt('../data/Comm/Gain.txt')
        comm_raw = (10 ** (comm_rawGdata / 10.0)).reshape(
            (361, 361), order='F')
        power_raw = np.genfromtxt('../data/Power/curve.dat')
        
        LDs = [5233.5, 5294.5, 5356.5, 5417.5, 5478.5, 5537.5]

        r_e2b_I0s = [np.array([4505.29362, -3402.16069, -3943.74582,
                               4.1923899, -1.56280012,  6.14347427]),
                     np.array(
                         [-1005.46693,  -597.205348, -6772.86532, -0.61047858,
                          -7.54623146,  0.75907455]),
                     np.array(
                         [4401.10539,  2275.95053, -4784.13188, -5.26605537,
                          -1.08194926, -5.37013745]),
                     np.array(
                         [-4969.91222,  4624.84149,  1135.9414,  0.1874654,
                          -1.62801666,  7.4302362]),
                     np.array(
                         [-235.021232,  2195.72976,  6499.79919, -2.55956031,
                          -6.82743519,  2.21628099]),
                     np.array(
                         [-690.314375, -1081.78239, -6762.90367,  7.44316722,
                          1.19745345, -0.96035904])]
        
        top = set_as_top(Assembly())
        top.add('pt', CADRE(n, m, solar_raw1, solar_raw2,
                              comm_raw, power_raw))        

        i = 0
        
        top.pt.set("LD", LDs[i])
        top.pt.set("r_e2b_I0", r_e2b_I0s[i])

        top.pt.run()
        
        inputs = ['BsplineParameters.CP_gamma']
        outputs = ['Comm_DataDownloaded.Data[0][-1]']        
            
        J1 = top.pt.driver.workflow.calc_gradient(inputs, outputs, mode='forward')
        
        top.pt.driver.workflow.config_changed()
        J2 = top.pt.driver.workflow.calc_gradient(inputs, outputs, mode='adjoint')
        
        top.pt.driver.workflow.config_changed()
        Jfd = top.pt.driver.workflow.calc_gradient(inputs, outputs, mode='fd')
        
        print np.max(J1-Jfd)
        print np.max(J2-Jfd)
        print np.max(J1-J2)
        
        self.assertTrue( np.max(J1-J2) < 1.0e-6 )
            
if __name__ == "__main__":

    unittest.main()