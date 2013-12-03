
import unittest
import numpy as np
import pickle
import random
import warnings

from openmdao.main.api import Assembly, set_as_top
from openmdao.util.testutil import assert_rel_error

from CADRE.attitude import Attitude_Angular, Attitude_AngularRates, \
    Attitude_Attitude, Attitude_Roll, Attitude_RotationMtx, \
    Attitude_RotationMtxRates, Attitude_Torque
from CADRE.battery import BatteryConstraints, BatteryPower, BatterySOC
from CADRE.comm import Comm_AntRotation, Comm_AntRotationMtx, Comm_BitRate, \
    Comm_DataDownloaded, Comm_Distance, Comm_EarthsSpin, Comm_EarthsSpinMtx, \
    Comm_GainPattern, Comm_GSposEarth, Comm_GSposECI, Comm_LOS, Comm_VectorAnt, \
    Comm_VectorBody, Comm_VectorECI, Comm_VectorSpherical
from CADRE.orbit import Orbit_Initial, Orbit_Dynamics
from CADRE.parameters import BsplineParameters
from CADRE.power import Power_CellVoltage, Power_SolarPower, Power_Total
from CADRE.reactionwheel import ReactionWheel_Power, \
    ReactionWheel_Torque, ReactionWheel_Dynamics
from CADRE.solar import Solar_ExposedArea
from CADRE.sun import Sun_LOS, Sun_PositionBody, Sun_PositionECI, \
    Sun_PositionSpherical
from CADRE.thermal_temperature import ThermalTemperature


import os

NTIME = 5

# Ignore the numerical warnings from performing the rel error calc.
warnings.simplefilter("ignore")


class Testcase_CADRE(unittest.TestCase):

    """ Test run/step/stop aspects of a simple workflow. """

    def setUp(self):
        """ Called before each test. """
        self.model = set_as_top(Assembly())

    def tearDown(self):
        """ Called after each test. """
        self.model = None

    def setup(self, compname, inputs, state0):

        try:
            self.model.add('comp', eval('%s(NTIME)' % compname))
        except TypeError:
            # At least one comp has no args.
            self.model.add('comp', eval('%s()' % compname))

        self.model.driver.workflow.add('comp')

        #data = pickle.load(open("data1346.pkl", 'rb'))

        for item in inputs + state0:
            val = self.model.comp.get(item)
            #key = "%d:%s" % (NTIME, item)
            if hasattr(val, 'shape'):
                shape1 = val.shape
                #self.model.comp.set(item, data[key])
                self.model.comp.set(item, np.random.random(shape1))
            else:
                #self.model.comp.set(item, data[key][0])
                self.model.comp.set(item, random.random())

    def run_model(self):

        self.model.comp.h = 0.01
        self.model.run()

    def compare_derivatives(self, var_in, var_out, rel_error=False):

        wflow = self.model.driver.workflow
        inputs = ['comp.%s' % v for v in var_in]
        outputs = ['comp.%s' % v for v in var_out]

        # Numeric
        wflow.config_changed()
        Jn = wflow.calc_gradient(inputs=inputs,
                                 outputs=outputs,
                                 mode="fd")
        print Jn
        print '\n'

        # Analytic forward
        wflow.config_changed()
        Jf = wflow.calc_gradient(inputs=inputs,
                                 outputs=outputs,
                                 mode='forward')

        print Jf
        print '\n'

        if rel_error:
            diff = np.nan_to_num(abs(Jf - Jn) / Jn)
        else:
            diff = abs(Jf - Jn)

        #assert_rel_error(self, diff.max(), 0.0, 1e-3)

        # Analytic adjoint
        wflow.config_changed()
        Ja = wflow.calc_gradient(inputs=inputs,
                                 outputs=outputs,
                                 mode='adjoint')

        print Ja

        if rel_error:
            diff = np.nan_to_num(abs(Ja - Jn) / Jn)
        else:
            diff = abs(Ja - Jn)

        assert_rel_error(self, diff.max(), 0.0, 1e-3)

    def test_AAOrbit_Dynamics(self):

        compname = 'Orbit_Dynamics'
        inputs = ['r_e2b_I0']
        outputs = ['r_e2b_I']
        state0 = []

        self.setup(compname, inputs, state0)
        
        shape = self.model.comp.r_e2b_I0.shape
        print shape
        self.model.comp.r_e2b_I0[:3] = np.random.random((3)) * 1e6
        self.model.comp.r_e2b_I0[3:] = np.random.random((3)) * 1e5
        
        self.run_model()
        self.compare_derivatives(inputs+state0, outputs)
        
if __name__ == "__main__":

    unittest.main()
