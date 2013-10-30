
# import unittest
# import sys
# from math import log
# import numpy as np
# from time import time
# import random

# from openmdao.main.api import Assembly, set_as_top
# from openmdao.util.testutil import assert_rel_error

# from CADRE.attitude import Attitude_Angular, Attitude_AngularRates, \
#     Attitude_Attitude, Attitude_Roll, Attitude_RotationMtx, \
#     Attitude_RotationMtxRates, Attitude_Sideslip, Attitude_Torque
# from CADRE.battery import BatteryConstraints, BatteryPower, BatterySOC
# from CADRE.comm import Comm_AntRotation, Comm_AntRotationMtx, Comm_BitRate, \
#     Comm_DataDownloaded, Comm_Distance, Comm_EarthsSpin, Comm_EarthsSpinMtx, \
#     Comm_GainPattern, Comm_GSposEarth, Comm_GSposECI, Comm_LOS, Comm_VectorAnt, \
#     Comm_VectorBody, Comm_VectorECI, Comm_VectorSpherical
# from CADRE.orbit import Orbit_Initial, Orbit_Dynamics
# from CADRE.parameters import BsplineParameters
# from CADRE.power import Power_CellVoltage, Power_SolarPower, Power_Total
# from CADRE.reactionwheel import ReactionWheel_Motor, ReactionWheel_Power, \
#     ReactionWheel_Torque, ReactionWheel_Dynamics
# from CADRE.solar import Solar_ExposedArea
# from CADRE.sun import Sun_LOS, Sun_PositionBody, Sun_PositionECI, \
#     Sun_PositionSpherical
# from CADRE.thermal_temperature import ThermalTemperature


# NTIME = 300
# NTIME2 = NTIME / 2
# NEXEC = 5


# class Testcase_CADRE(unittest.TestCase):

#     """ Test run/step/stop aspects of a simple workflow. """

#     def setUp(self):
#         """ Called before each test. """
#         self.model = set_as_top(Assembly())
#         self.model2 = Assembly()

#     def tearDown(self):
#         """ Called after each test. """
#         self.model = None
#         self.model2 = None

#     def setup(self, compname, inputs, state0):

#         try:
#             self.model.add('comp', eval('%s(NTIME)' % compname))
#         except TypeError:
# At least one comp has no args.
#             self.model.add('comp', eval('%s()' % compname))

#         try:
#             self.model2.add('comp', eval('%s(NTIME2)' % compname))
#         except TypeError:
# At least one comp has no args.
#             self.model2.add('comp', eval('%s()' % compname))

#         self.model.driver.workflow.add('comp')
#         self.model2.driver.workflow.add('comp')

#         for item in inputs + state0:
#             val = self.model.comp.get(item)
#             val2 = self.model2.comp.get(item)
#             if hasattr(val, 'shape'):
#                 shape1 = val.shape
#                 self.model.comp.set(item, np.random.random(shape1))
#                 shape2 = val2.shape
#                 self.model2.comp.set(item, np.random.random(shape2))
#             else:
#                 self.model.comp.set(item, random.random())
#                 self.model2.comp.set(item, random.random())

#     def compare_derivatives(self, var_in, var_out, rel_error=False):

#         inputs = ['comp.%s' % v for v in var_in]
#         outputs = ['comp.%s' % v for v in var_out]
#         self.model.comp.h = 0.01
#         self.model2.comp.h = 0.01

#         print '\n'
#         print type(self.model.comp), "n=%d, averages=%d" % (NTIME, NEXEC)
#         print 30 * '-'

#         tzero = time()
#         for i in range(NEXEC):
#             self.model.run()

#         tt = (time() - tzero) / float(NEXEC)

#         tzero = time()
#         for i in range(NEXEC):
#             self.model2.run()

#         tt2 = (time() - tzero) / float(NEXEC)
#         exp = log(tt / tt2, 2)

#         print "Execution:   ", tt, "(n ** %5.2f)" % exp
#         sys.stdout.flush()

#         self.model.driver.workflow.config_changed()
#         tzero = time()
#         for i in range(NEXEC):
#             self.model.comp.linearize()

#         tt = (time() - tzero) / float(NEXEC)

#         self.model2.driver.workflow.config_changed()
#         tzero = time()
#         for i in range(NEXEC):
#             self.model2.comp.linearize()

#         tt2 = (time() - tzero) / float(NEXEC)
#         exp = log(tt / tt2, 2)

#         print "Linearize:   ", tt, "(n ** %5.2f)" % exp
#         sys.stdout.flush()

#         self.model.driver.workflow.config_changed()
#         tzero = time()
#         for i in range(NEXEC):
#             self.model.driver.workflow.calc_gradient(
#                 inputs=inputs, outputs=outputs)

#         tt = (time() - tzero) / float(NEXEC)

#         self.model2.driver.workflow.config_changed()
#         tzero = time()
#         for i in range(NEXEC):
#             self.model2.driver.workflow.calc_gradient(
#                 inputs=inputs, outputs=outputs)

#         tt2 = (time() - tzero) / float(NEXEC)
#         exp = log(tt / tt2, 2)

#         print "Apply_J  :   ", tt, "(n ** %5.2f)" % exp
#         sys.stdout.flush()

#         self.model.driver.workflow.config_changed()
#         tzero = time()
#         for i in range(NEXEC):
#             self.model.driver.workflow.calc_gradient(
#                 inputs=inputs, outputs=outputs,
#                 mode='adjoint')

#         tt = (time() - tzero) / float(NEXEC)

#         self.model2.driver.workflow.config_changed()
#         tzero = time()
#         for i in range(NEXEC):
#             self.model2.driver.workflow.calc_gradient(
#                 inputs=inputs, outputs=outputs,
#                 mode='adjoint')

#         tt2 = (time() - tzero) / float(NEXEC)
#         exp = log(tt / tt2, 2)

#         print "Apply_JT :   ", tt, "(n ** %5.2f)" % exp
#         sys.stdout.flush()

#     def run_model(self):
#         pass

#     def test_Comm_DataDownloaded(self):

#         compname = 'Comm_DataDownloaded'
#         inputs = ['Dr']
#         outputs = ['Data']
#         state0 = ['Data0']

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_AntRotation(self):

#         compname = 'Comm_AntRotation'
#         inputs = ['antAngle']
#         outputs = ['q_A']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_AntRotationMtx(self):

#         compname = 'Comm_AntRotationMtx'
#         inputs = ['q_A']
#         outputs = ['O_AB']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_BitRate(self):

#         compname = 'Comm_BitRate'
#         inputs = ['P_comm', 'gain', 'GSdist', 'CommLOS']
#         outputs = ['Dr']
#         state0 = []

#         self.setup(compname, inputs, state0)

# These need to be a certain magnitude so it doesn't blow up
#         shape = self.model.comp.P_comm.shape
#         self.model.comp.P_comm = np.ones(shape)
#         shape = self.model.comp.GSdist.shape
#         self.model.comp.GSdist = np.random.random(shape) * 1e3

#         self.run_model()

#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_Distance(self):

#         compname = 'Comm_Distance'
#         inputs = ['r_b2g_A']
#         outputs = ['GSdist']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_EarthsSpin(self):

#         compname = 'Comm_EarthsSpin'
#         inputs = ['t']
#         outputs = ['q_E']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_EarthsSpinMtx(self):

#         compname = 'Comm_EarthsSpinMtx'
#         inputs = ['q_E']
#         outputs = ['O_IE']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_GainPattern(self):

#         compname = 'Comm_GainPattern'
#         inputs = ['azimuthGS', 'elevationGS']
#         outputs = ['gain']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_GSposEarth(self):

#         compname = 'Comm_GSposEarth'
#         inputs = ['lon', 'lat', 'alt']
#         outputs = ['r_e2g_E']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_GSposECI(self):

#         compname = 'Comm_GSposECI'
#         inputs = ['O_IE', 'r_e2g_E']
#         outputs = ['r_e2g_I']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_LOS(self):

#         compname = 'Comm_LOS'
#         inputs = ['r_b2g_I', 'r_e2g_I']
#         outputs = ['CommLOS']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_VectorAnt(self):

#         compname = 'Comm_VectorAnt'
#         inputs = ['r_b2g_B', 'O_AB']
#         outputs = ['r_b2g_A']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_VectorBody(self):

#         compname = 'Comm_VectorBody'
#         inputs = ['r_b2g_I', 'O_BI']
#         outputs = ['r_b2g_B']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_VectorECI(self):

#         compname = 'Comm_VectorECI'
#         inputs = ['r_e2g_I', 'r_e2b_I']
#         outputs = ['r_b2g_I']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Comm_VectorSpherical(self):

#         compname = 'Comm_VectorSpherical'
#         inputs = ['r_b2g_A']
#         outputs = ['azimuthGS', 'elevationGS']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_ThermalTemperature(self):

#         compname = 'ThermalTemperature'
#         inputs = ['exposedArea', 'cellInstd', 'LOS', 'P_comm']
#         outputs = ['temperature']
#         state0 = ['T0']

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Attitude_Angular(self):

#         compname = 'Attitude_Angular'
#         inputs = ['O_BI', 'Odot_BI']
#         outputs = ['w_B']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Attitude_AngularRates(self):

#         compname = 'Attitude_AngularRates'
#         inputs = ['w_B']
#         outputs = ['wdot_B']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Attitude_Attitude(self):

#         compname = 'Attitude_Attitude'
#         inputs = ['r_e2b_I']
#         outputs = ['O_RI']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Attitude_Roll(self):

#         compname = 'Attitude_Roll'
#         inputs = ['Gamma']
#         outputs = ['O_BR']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Attitude_RotationMtx(self):

#         compname = 'Attitude_RotationMtx'
#         inputs = ['O_BR', 'O_RI']
#         outputs = ['O_BI']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Attitude_RotationMtxRates(self):

#         compname = 'Attitude_RotationMtxRates'
#         inputs = ['O_BI']
#         outputs = ['Odot_BI']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Attitude_Sideslip(self):

#         compname = 'Attitude_Sideslip'
#         inputs = ['r_e2b_I', 'O_BI']
#         outputs = ['v_e2b_B']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Attitude_Torque(self):

#         compname = 'Attitude_Torque'
#         inputs = ['w_B', 'wdot_B']
#         outputs = ['T_tot']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Sun_LOS(self):

#         compname = 'Sun_LOS'
#         inputs = ['r_e2b_I', 'r_e2s_I']
#         outputs = ['LOS']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Sun_PositionBody(self):

#         compname = 'Sun_PositionBody'
#         inputs = ['O_BI', 'r_e2s_I']
#         outputs = ['r_e2s_B']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Sun_PositionECI(self):

#         compname = 'Sun_PositionECI'
#         inputs = ['t', 'LD']
#         outputs = ['r_e2s_I']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Sun_PositionSpherical(self):

#         compname = 'Sun_PositionSpherical'
#         inputs = ['r_e2s_B']
#         outputs = ['azimuth', 'elevation']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Solar_ExposedArea(self):

#         compname = 'Solar_ExposedArea'
#         inputs = ['finAngle', 'azimuth', 'elevation']
#         outputs = ['exposedArea']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Power_CellVoltage(self):
# fix
#         compname = 'Power_CellVoltage'
#         inputs = ['LOS', 'temperature', 'exposedArea', 'Isetpt']
#         outputs = ['V_sol']
#         state0 = []

#         self.setup(compname, inputs, state0)

#         shape1 = self.model.comp.temperature.shape
#         shape2 = self.model2.comp.temperature.shape
#         self.model.comp.temperature = np.ones(shape1)
#         self.model2.comp.temperature = np.ones(shape2)

#         shape1 = self.model.comp.temperature.shape
#         shape2 = self.model2.comp.temperature.shape
#         self.model.comp.temperature = np.random.random(shape1) * 40 + 240
#         self.model2.comp.temperature = np.random.random(shape2) * 40 + 240

#         shape1 = self.model.comp.exposedArea.shape
#         shape2 = self.model2.comp.exposedArea.shape
#         self.model.comp.exposedArea = np.random.random(shape1) * 1e-4
#         self.model2.comp.exposedArea = np.random.random(shape2) * 1e-4

#         shape1 = self.model.comp.Isetpt.shape
#         shape2 = self.model2.comp.Isetpt.shape
#         self.model.comp.Isetpt = np.random.random(shape1) * 1e-6
#         self.model2.comp.Isetpt = np.random.random(shape2) * 1e-6

#         self.run_model()
#         self.compare_derivatives(inputs, outputs, rel_error=True)

#     def test_Power_SolarPower(self):

#         compname = 'Power_SolarPower'
#         inputs = ['V_sol', 'Isetpt']
#         outputs = ['P_sol']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Power_Total(self):

#         compname = 'Power_Total'
#         inputs = ['P_sol', 'P_comm', 'P_RW']
#         outputs = ['P_bat']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_ReactionWheel_Motor(self):

#         compname = 'ReactionWheel_Motor'
#         inputs = ['T_RW', 'w_B', 'w_RW']
#         outputs = ['T_m']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_ReactionWheel_Dynamics(self):

#         compname = 'ReactionWheel_Dynamics'
#         inputs = ['w_B', 'T_RW']
#         outputs = ['w_RW']

# keep these at zeros
# state0 = []  # ['w_RW0']

#         self.setup(compname, inputs, state0)

#         shape = self.model.comp.w_B.shape
#         self.model.comp.w_B = np.random.random(shape) * 1e-4
#         shape = self.model.comp.T_RW.shape
#         self.model.comp.T_RW = np.random.random(shape) * 1e-9

#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_ReactionWheel_Power(self):

#         compname = 'ReactionWheel_Power'
#         inputs = ['w_RW', 'T_RW']
#         outputs = ['P_RW']
#         state0 = []

#         self.setup(compname, inputs, state0)

#         self.run_model()
#         self.compare_derivatives(inputs, outputs, rel_error=True)

#     def test_ReactionWheel_Torque(self):

#         compname = 'ReactionWheel_Torque'
#         inputs = ['T_tot']
#         outputs = ['T_RW']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_BatterySOC(self):

#         compname = 'BatterySOC'
#         inputs = ['P_bat', 'temperature']
#         outputs = ['SOC']
#         state0 = ['iSOC']

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_BatteryPower(self):

#         compname = 'BatteryPower'
#         inputs = ['SOC', 'temperature', 'P_bat']
#         outputs = ['I_bat']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_BatteryConstraints(self):

#         compname = 'BatteryConstraints'
#         inputs = ['I_bat', 'SOC']
#         outputs = ['ConCh', 'ConDs', 'ConS0', 'ConS1']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_Orbit_Dynamics(self):
# This comp has no inputs, so no derivs.
#         pass

# compname = 'Orbit_Dynamics'
# inputs = ['']
# outputs = ['r_e2b_I']
# state0 = ['r_e2b_I0']

# self.setup(compname, inputs, state0)
# self.run_model()
# self.compare_derivatives(inputs, outputs)

#     def test_Orbit_Initial(self):

#         compname = 'Orbit_Initial'
#         inputs = ['altPerigee', 'altApogee', 'RAAN', 'Inc', 'argPerigee',
#                   'trueAnomaly']
#         outputs = ['r_e2b_I0']
#         state0 = []

#         self.setup(compname, inputs, state0)
#         self.run_model()
#         self.compare_derivatives(inputs, outputs)

#     def test_bspline_parameters(self):

#         compname = 'BsplineParameters'
#         inputs = ['CP_P_comm', 'CP_gamma', 'CP_Isetpt']
#         outputs = ['P_comm', 'Gamma', 'Isetpt']
#         state0 = []

#         self.model.add('comp', BsplineParameters(NTIME, 5))
#         self.model.driver.workflow.add('comp')
#         self.model2.add('comp', BsplineParameters(NTIME2, 5))
#         self.model2.driver.workflow.add('comp')

#         for item in inputs:
#             val = self.model.comp.get(item)
#             shape1 = val.shape
#             self.model.comp.set(item, np.random.random(shape1))
#             val2 = self.model2.comp.get(item)
#             shape1 = val2.shape
#             self.model2.comp.set(item, np.random.random(shape1))

#         self.run_model()
#         self.compare_derivatives(inputs, outputs)


# if __name__ == "__main__":

#     unittest.main()
