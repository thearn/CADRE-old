''' Battery discipline for CADRE '''

import numpy as np

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import KS
import rk4

# Allow non-standard variable names for scientific calc
# pylint: disable-msg=C0103

#Constants
sigma = 1e-10
eta = 0.99
Cp = 2900.0*0.001*3600.0
IR = 0.9
T0 = 293.0
alpha = np.log(1/1.1**5)


class BatterySOC(rk4.RK4):
    """Computes the time history of the battery state of charge.

    n_time: int
        number of time_steps to take

    time_step: float
        size of each timestep
    """

    def __init__(self, n_times):

        super(BatterySOC, self).__init__()

        # Inputs
        self.add('iSOC',
            Array([0.0], shape=(1, ), dtype=np.float,
                iotype="in", units="unitless", desc="Initial state of charge")
        )

        self.add('P_bat',
            Array(np.zeros((n_times, )), shape=(n_times, ), dtype=np.float,
                units="W", iotype="in", desc="Battery power over time")
        )

        self.add('temperature',
            Array(np.zeros((5, n_times )), shape=(5, n_times ), dtype=np.float,
                units="degK", iotype="in", desc="Battery temperature over time")
        )

        # Outputs
        self.add('SOC',
            Array(np.zeros((1, n_times)), shape=(1, n_times), dtype=np.float,
                iotype="out", units="unitless", desc="Battery state of charge over time")
        )


        self.state_var = "SOC"
        self.init_state_var = "iSOC"
        self.external_vars = ["P_bat", "temperature"]


    def f_dot(self, external, state):
        """Rate of change of SOC"""

        SOC = state[0]
        P = external[0]
        T = external[5]


        voc = 3 + np.expm1(SOC) / (np.e-1)
        #dVoc_dSOC = np.exp(SOC) / (np.e-1)

        V = IR * voc * (2.0 - np.exp(alpha*(T-T0)/T0))
        I = P/V

        soc_dot = -sigma/24*SOC + eta/Cp*I
        return soc_dot

    def df_dy(self, external, state):
        """ State deriative """

        SOC = state[0]
        P = external[0]
        T = external[5]

        voc = 3 + np.expm1(SOC) / (np.e-1)
        dVoc_dSOC = np.exp(SOC) / (np.e-1)

        tmp = 2 - np.exp(alpha*(T-T0)/T0)
        V = IR * voc * tmp
        #I = P/V

        dV_dSOC = IR * dVoc_dSOC * tmp
        dI_dSOC = -P/V**2 * dV_dSOC

        df_dy = -sigma/24 + eta/Cp*dI_dSOC

        return np.array([[df_dy]])

    def df_dx(self, external, state):
        """ Output deriative """

        SOC = state[0]
        P = external[0]
        T = external[1]

        voc = 3 + np.expm1(SOC) / (np.e-1)
        #dVoc_dSOC = np.exp(SOC) / (np.e-1)

        tmp = 2 - np.exp(alpha*(T-T0)/T0)

        V = IR * voc * tmp
        #I = P/V

        dV_dT = - IR * voc * np.exp(alpha*(T-T0)/T0) * alpha/T0
        dI_dT = - P/V**2 * dV_dT
        dI_dP = 1.0/V

        return np.array([[eta/Cp*dI_dP, 0 , 0, 0, 0, eta/Cp*dI_dT]])


class BatteryPower(Component):
    """ Power supplied by the battery"""

    def __init__(self, n=2):
        super(BatteryPower, self).__init__()

        self.n = n

        # Inputs
        self.add('SOC', Array(np.zeros((1, n)), size=(1, n), dtype=np.float,
                              iotype="in", units="unitless", desc="Battery state of charge over time"))

        self.add('temperature', Array(np.zeros((5, n)), size=(n, ),
                                      dtype=np.float, iotype="in",
                                      units="degK", desc="Battery temperature over time"))

        self.add('P_bat', Array(np.zeros((n, )), size=(n, ), dtype=np.float,
                                iotype="in", units="W", desc="Battery power over time"))

        # Outputs
        self.add('I_bat', Array(np.zeros((n, )), size=(n, ), dtype=np.float,
                                iotype="out", units="A", desc="Battery Current over time"))

    def execute(self):
        """ Calculate output. """

        self.exponential = (2.0 - np.exp(alpha*(self.temperature[4, :]-T0)/T0))
        self.voc = 3.0 + np.expm1(self.SOC[0, :]) / (np.e-1)
        self.V = IR * self.voc * self.exponential
        self.I_bat = self.P_bat/self.V

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        #dI_dP
        dV_dvoc = IR * self.exponential
        dV_dT = - IR * self.voc * np.exp(alpha*(self.temperature[4, :] -
                                                T0)/T0) * alpha / T0
        dVoc_dSOC = np.exp(self.SOC[0, :]) / (np.e-1)

        self.dI_dP = 1.0 / self.V
        tmp = -self.P_bat/(self.V**2)
        self.dI_dT = tmp * dV_dT
        self.dI_dSOC = tmp * dV_dvoc * dVoc_dSOC

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'I_bat' in result:
            
            if 'P_bat' in arg:
                result['I_bat'] += self.dI_dP * arg['P_bat']
                
            if 'temperature' in arg:
                result['I_bat'] += self.dI_dT * arg['temperature'][4, :]
                
            if 'SOC' in arg:
                result['I_bat'] += self.dI_dSOC * arg['SOC'][0, :]

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'I_bat' in arg:
            
            if 'P_bat' in result:
                result['P_bat'] += self.dI_dP * arg['I_bat']

            if 'temperature' in result:
                result['temperature'] += np.zeros(self.temperature.shape)
                result['temperature'][4, :] += self.dI_dT * arg['I_bat']

            if 'SOC' in result:
                result['SOC'] += np.zeros(self.SOC.shape)
                result['SOC'][0, :] += self.dI_dSOC * arg['I_bat']



class BatteryConstraints(Component):
    """ Some KS constraints for the battery. I believe this essentially
    replaces a cycle in the graph.
    """

    # Outputs
    ConCh = Float(0.0,
		  iotype="out",
                  units="A",
		  desc="Constraint on charging rate")

    ConDs = Float(0.0,
		  iotype="out",
                  units="A",
		  desc="Constraint on discharging rate")

    ConS0 = Float(0.0,
		  iotype="out",
		  units="unitless",
                  desc="Constraint on minimum state of charge")

    ConS1 = Float(0.0,
		  iotype="out",
		  units="unitless",
                  desc="Constraint on maximum state of charge")

    def __init__(self, n=2):
        super(BatteryConstraints, self).__init__()
        self.n = n

        self.rho = 50
        self.Imin = -10.0
        self.Imax = 5.0
        self.SOC0 = 0.2
        self.SOC1 = 1.0

        # Inputs
        self.add('I_bat', Array(np.zeros((n,)),
		                size=(n,),
				iotype="in",
				units="A",
                                desc="Battery current over time"))

        self.add('SOC', Array(np.zeros((1, n)),
		              size=(1, n),
			      iotype="in",
			      units="unitless",
                              desc="Battery state of charge over time"))

        self.KS_ch = KS.KSfunction()
        self.KS_ds = KS.KSfunction()
        self.KS_s0 = KS.KSfunction()
        self.KS_s1 = KS.KSfunction()

    def execute(self):
        """ Calculate output. """

        self.ConCh = self.KS_ch.compute(self.I_bat - self.Imax, self.rho)
        self.ConDs = self.KS_ds.compute(self.Imin - self.I_bat, self.rho)
        self.ConS0 = self.KS_s0.compute(self.SOC0 - self.SOC, self.rho)
        self.ConS1 = self.KS_s1.compute(self.SOC - self.SOC1, self.rho)

    def linearize(self):
        """ Calculate and save derivatives. (i.e., Jacobian) """

        self.dCh_dg, self.dCh_drho = self.KS_ch.derivatives()
        self.dDs_dg, self.dDs_drho = self.KS_ds.derivatives()
        self.dS0_dg, self.dS0_drho = self.KS_s0.derivatives()
        self.dS1_dg, self.dS1_drho = self.KS_s1.derivatives()

    def apply_deriv(self, arg, result):
        """ Matrix-vector product with the Jacobian. """

        if 'I_bat' in arg:
            if 'ConCh' in result:
                result['ConCh'] += np.dot(self.dCh_dg, arg['I_bat'])
            if 'ConDs' in result:
                result['ConDs'] -= np.dot(self.dDs_dg, arg['I_bat'])
                
        if 'SOC' in arg:
            if 'ConS0' in result:
                result['ConS0'] -= np.dot(self.dS0_dg, arg['SOC'][0, :])
            if 'ConS1' in result:
                result['ConS1'] += np.dot(self.dS1_dg, arg['SOC'][0, :])

    def apply_derivT(self, arg, result):
        """ Matrix-vector product with the transpose of the Jacobian. """

        if 'I_bat' in result:
            if 'ConCh' in arg:
                result['I_bat'] += self.dCh_dg * arg['ConCh']
            if 'ConDs' in arg:
                result['I_bat'] -= self.dDs_dg * arg['ConDs']
                
        if 'SOC' in result:
            if 'ConS0' in arg:
                result['SOC'] -= self.dS0_dg * arg['ConS0']
            if 'ConS1' in arg:
                result['SOC'] += self.dS1_dg * arg['ConS1']




