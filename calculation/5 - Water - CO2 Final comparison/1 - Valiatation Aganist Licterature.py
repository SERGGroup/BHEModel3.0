# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_classes import BaseBHE, ReservoirProperties, BHEGeometry
from main_classes.constant import CALCULATION_DIR
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import numpy as np
import os


# %%------------   IMPORT VALIDATION DATA AND INITIALIZE CALCULATION    ---------------------------------------------> #

#   Validation data from:
#       "Techno-economic analysis of Advanced Geothermal Systems (AGS)" (in "0 - Resources/0 - Validation/res" folder)
#
#   Data extracted from plot using:
#       WebPlotDigitizer (https://automeris.io/WebPlotDigitizer)
#

# Input Parameters
# General Params
t_in = 10                   # [°C]
time = 30                   # [year]
s_year = 365 * 24 * 3600    # [s / year]
flow_rates_ovr = np.linspace(0, 200, 50)

# Geometric Params
bhe_geom = BHEGeometry()
bhe_geom.depth = 3500        # [m]
bhe_geom.l_horiz = 5000      # [m]
bhe_geom.d_well = 0.5        # [m]
bhe_geom.n_wells = 4

# Reservoir Params
res_prop = ReservoirProperties()
res_prop.grad = 0.035       # [°C/m]
res_prop.k_rock = 2.1       # [W/(m K)]
res_prop.c_rock = 1000      # [J/(kg K)]
res_prop.rho_rock = 2650    # [kg/m^3]

base_res_folder = os.path.join(CALCULATION_DIR, '5 - Water - CO2 Final comparison', '0 - Resources', '0 - Validation')
dataset_file = os.path.join(base_res_folder, 'res', 'digitalized_dataset.csv')

data = np.genfromtxt(dataset_file, delimiter=',', skip_header=2, filling_values=np.nan)
temperatures = data[:, :2]
powers = data[:, 2:]

# For internal reason in WebPlotDigitizer the power data has been extracted using the same axes as the temperature
# hence it has to be converted back (instead of MW as in the graph, power has been converted in kW)
powers[:, 1] = powers[:, 1] / 120 * 400


# %%------------   EVALUATE CO2                                         ---------------------------------------------> #
bhe_in = ThermodynamicPoint(["CarbonDioxide"], [1], unit_system="MASS BASE SI")
bhe_in.set_variable("T", t_in + 273.15)
bhe_in.set_variable("Q", 0)

CO2_bhe = BaseBHE(bhe_in, reservoir_properties=res_prop, geometry=bhe_geom)
CO2_bhe.set_HX_condition()
CO2_dict = CO2_bhe.evaluate_HXG(times=time * s_year, m_dot=flow_rates_ovr)
CO2_dict = CO2_bhe.evaluate_surface_condition(CO2_dict)
CO2_dict = CO2_bhe.evaluate_direct_expansion(CO2_dict)


# %%------------   EVALUATE WATER                                       ---------------------------------------------> #
bhe_in = ThermodynamicPoint(["Water"], [1], unit_system="MASS BASE SI")
bhe_in.set_variable("T", t_in + 273.15)
bhe_in.set_variable("P", 1e5)

water_bhe = BaseBHE(bhe_in, reservoir_properties=res_prop, geometry=bhe_geom)
water_bhe.set_HX_condition()
water_dict = water_bhe.evaluate_HXG(times=time * s_year, m_dot=flow_rates_ovr)
water_dict = water_bhe.evaluate_surface_condition(water_dict)


# %%------------   PLOT TEMPERATURE                                     ---------------------------------------------> #
plt.plot(water_dict["m_dot"], (water_dict["T_out"] - 273.15), label="water")
plt.plot(CO2_dict["m_dot"], (CO2_dict["T_out"] - 273.15), label="CO2")
plt.scatter(temperatures[:, 0], temperatures[:, 1], c="black", marker="+", label="Malek et Al.")
plt.xlim([0, 100])
plt.legend()
plt.show()


# %%------------   PLOT POWER                                           ---------------------------------------------> #
eta_CO2 = 0.8
eta_water = 0.15
# plt.plot(water_dict["m_dot"]/1.15, water_dict["dex"] * water_dict["m_dot"] * eta_water / 1e3, label="water")
plt.plot(CO2_dict["m_dot"] / 2, CO2_dict["dh_turb"] * CO2_dict["m_dot"] * eta_CO2 / 1e3, label="CO2")
plt.scatter(powers[:, 0], powers[:, 1], c="black", marker="+", label="Malek et Al.")
plt.xlim([0, 100])
plt.legend()
plt.show()
