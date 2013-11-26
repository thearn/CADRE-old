''' Thermal discipline for CADRE '''

import numpy as np

from openmdao.lib.datatypes.api import Array

from CADRE.rk4 import RK4

# Allow non-standard variable names for scientific calc
# pylint: disable-msg=C0103

#constants
m_f = 0.4
m_b = 2.0
cp_f = 0.6e3
cp_b = 2.0e3
A_T = 2.66e-3

alpha_c = 0.9
alpha_r = 0.2
eps_c = 0.87
eps_r = 0.88

q_sol = 1360.0
K = 5.67051e-8

class ThermalTemperature(RK4):
    """ Calculates the temperature distribution on the solar panels."""

    def __init__(self, n_times):
        super(ThermalTemperature, self).__init__()

        # Inputs
        self.add("T0", Array(273.*np.ones((5,)),
		             shape=(5,),
			     dtype=np.float,
			     units="degK",
                             iotype="in",
			     desc="initial temperatures for the 4 fins and body")
        )

        self.add("exposedArea", Array(np.zeros((7, 12, n_times)),
                                      units="m**2",
                                      size=(7, 12 ,n_times),
				      dtype=np.float,
                                      iotype="in",
                                      desc="exposed area for each solar cell")
        )

        self.add("cellInstd", Array(np.ones((7, 12)),
		                    size=(7, 12),
                                    dtype=np.float,
				    iotype="in",
                                    units='unitless',
                                    desc="Cell/Radiator indication",
                                    low=0,
				    high=1)
        )

        self.add("LOS", Array(np.zeros((n_times, )),
		              size=(n_times, ),
                              dtype=np.float,
			      iotype="in",
                              units='unitless',
                              desc="Line of sight to the sun",
			      low=0,
			      high=1)
        )

        self.add("P_comm", Array(np.ones((n_times, )),
		                 size=(n_times, ),
                                 dtype=np.float,
				 iotype="in",
                                 units='W',
                                 desc="Power required by the communication system",
				 low=0,
				 high=1)
        )

        # Outputs
        self.add("temperature", Array(np.zeros((5, n_times)),
                                      shape=(5, n_times), dtype=np.float,
                                      units="degK",
                                      iotype="out", desc="temperature for the 4 fins and body over time",
                                      low=50, high=400)
        )

        self.state_var = "temperature"
        self.init_state_var = "T0"
        self.external_vars = ["exposedArea", "LOS", "P_comm"]
        self.fixed_external_vars = ["cellInstd",]

        # implementation of fixTemps from Thermal_Temperature.f90
        for i in range (0, n_times):
            for k in range (0, 5):
                self.temperature[k, i] = self.T0[k]
                if self.temperature[k, i] < 0:
                    self.temperature[k, i] = 0.

    def f_dot(self, external, state):

        # revised implementation from ThermalTemperature.f90
        exposedArea = external[:84].reshape(7, 12)
        LOS = external[84]
        P_comm = external[85]
        cellInstd = external[86:].reshape(7, 12)

        f = np.zeros((5, ))

        alpha = alpha_c*cellInstd + alpha_r - alpha_r*cellInstd
        eps = eps_c*cellInstd + eps_r - eps_r*cellInstd

        # Panels
        for p in range(0, 12):

            # Body
            if p < 4:
                f_i = 4
                m = m_b
                cp = cp_b

            # Fin
            else:
                f_i = (p+1)%4
                m = m_f
                cp = cp_f
            
            # Cells
            fact1 = q_sol*LOS/(m*cp)
            fact2 = K*A_T*state[f_i]**4/(m*cp)
            f[f_i] += np.sum(alpha[:, p] * exposedArea[:, p]) * fact1
            f[f_i] -= np.sum(eps[:, p]) * fact2

        f[4] += 4.0 * P_comm / m_b / cp_b

        return f


    def df_dy(self, external, state):

        # revised implementation from ThermalTemperature.f90
        cellInstd = external[86:].reshape(7, 12)
        eps = eps_c*cellInstd + eps_r - eps_r*cellInstd

        dfdy = np.zeros((5, 5))

        # Panels
        for p in range(0, 12):

            # Body
            if p < 4:
                f_i = 4
                m = m_b
                cp = cp_b

            # Fin
            else:
                f_i = (p+1)%4
                m = m_f
                cp = cp_f

            # Cells
            fact = 4.0*K*A_T*state[f_i]**3/(m*cp)
            dfdy[f_i, f_i] -= np.sum(eps) * fact

        return dfdy

    def df_dx(self, external, state):

        # revised implementation from ThermalTemperature.f90
        exposedArea = external[:84].reshape(7, 12)
        LOS = external[84]
        cellInstd = external[86:].reshape(7, 12)

        dfdx = np.zeros((5, 170))
        alpha = alpha_c*cellInstd + alpha_r - alpha_r*cellInstd
        dalpha_dw = alpha_c - alpha_r
        deps_dw = eps_c - eps_r

        # Panels
        for p in range(0, 12):

            # Body
            if p < 4:
                f_i = 4
                m = m_b
                cp = cp_b

            # Fin
            else:
                f_i = (p+1)%4
                m = m_f
                cp = cp_f

            # Cells
            fact3 = q_sol/(m*cp)
            fact1 = fact3*LOS
            fact2 = K*A_T*state[f_i]**4/(m*cp)
                
            dfdx[f_i, p:p+84:12] += alpha[:, p] * fact1
            dfdx[f_i, p+86:p+170:12] += dalpha_dw * exposedArea[:, p] * fact1 - \
                                          deps_dw * fact2
            dfdx[f_i, 84] += np.sum(alpha[:, p] * exposedArea[:, p]) * fact3

        dfdx[4, 85] += 4.0 / m_b / cp_b

        return dfdx







