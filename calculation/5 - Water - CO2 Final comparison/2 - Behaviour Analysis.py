# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_classes import BaseBHE, ReservoirProperties, BHEGeometry
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import numpy as np


# %%------------   IMPORT VALIDATION DATA AND INITIALIZE CALCULATION    ---------------------------------------------> #
# Input Parameters
# General Params
t_in = 10                   # [°C]
time = 20                   # [year]
s_year = 365 * 24 * 3600    # [s / year]
flow_rates_ovr = np.linspace(0, 200, 101)[1:]

# Geometric Params
bhe_geom = BHEGeometry()
bhe_geom.depth = 1800        # [m]
bhe_geom.l_horiz = 5000      # [m]
bhe_geom.d_well = 0.5        # [m]
bhe_geom.n_wells = 4

# Reservoir Params
res_prop = ReservoirProperties()
res_prop.grad = 0.035       # [°C/m]
res_prop.k_rock = 2.1       # [W/(m K)]
res_prop.c_rock = 1000      # [J/(kg K)]
res_prop.rho_rock = 2650    # [kg/m^3]


# %%------------   EVALUATE CO2                                         ---------------------------------------------> #
bhe_in = ThermodynamicPoint(["CarbonDioxide"], [1], unit_system="MASS BASE SI")
bhe_in.set_variable("T", t_in + 273.15)
bhe_in.set_variable("Q", 0)
p_in_co2 = bhe_in.get_variable("P")

CO2_bhe = BaseBHE(bhe_in, reservoir_properties=res_prop, geometry=bhe_geom)
CO2_bhe.set_HX_condition()
CO2_dict = CO2_bhe.evaluate_HXG(times=time * s_year, m_dot=flow_rates_ovr)
CO2_dict = CO2_bhe.evaluate_surface_condition(CO2_dict)
CO2_dict = CO2_bhe.evaluate_direct_expansion(CO2_dict)


# %%------------   EVALUATE WATER                                       ---------------------------------------------> #
bhe_in = ThermodynamicPoint(["Water"], [1], unit_system="MASS BASE SI")
bhe_in.set_variable("T", t_in + 273.15)
bhe_in.set_variable("P", 1e5)
p_in_water = bhe_in.get_variable("P")

water_bhe = BaseBHE(bhe_in, reservoir_properties=res_prop, geometry=bhe_geom)
water_bhe.set_HX_condition()
water_dict = water_bhe.evaluate_HXG(times=time * s_year, m_dot=flow_rates_ovr)
water_dict = water_bhe.evaluate_surface_condition(water_dict)


# %%------------   PLOT TEMPERATURE                                     ---------------------------------------------> #
plt.plot(water_dict["m_dot"], (water_dict["T_out"] - 273.15), label="water")
plt.plot(CO2_dict["m_dot"], (CO2_dict["T_out"] - 273.15), label="CO2")
plt.legend()
plt.show()


# %%------------   PLOT Pressure                                        ---------------------------------------------> #
plt.plot(water_dict["m_dot"], (water_dict["p_out"]-p_in_water)/1e5, label="water")
plt.plot(CO2_dict["m_dot"], (CO2_dict["p_out"]-p_in_co2)/1e5, label="CO2")
plt.legend()
plt.show()


# %%------------   PLOT DH perc                                         ---------------------------------------------> #
plt.plot(water_dict["m_dot"], water_dict["dh_perc"], label="Water")
plt.plot(CO2_dict["m_dot"], CO2_dict["dh_perc"], label="CO2")
plt.legend()
plt.show()

# %%------------   PLOT rho                                             ---------------------------------------------> #
plt.plot(water_dict["m_dot"], water_dict["drho"], label="Water")
plt.plot(CO2_dict["m_dot"], CO2_dict["drho"], label="CO2")
plt.legend()
plt.show()


# %%------------   PLOT EXERGY                                          ---------------------------------------------> #
plt.plot(water_dict["m_dot"], water_dict["dex"] * water_dict["m_dot"] / 1e3, label="Water")
plt.plot(CO2_dict["m_dot"], CO2_dict["dex"] * CO2_dict["m_dot"] / 1e3, label="CO2")
plt.plot(CO2_dict["m_dot"], CO2_dict["dh_turb"] * CO2_dict["m_dot"] / 1e3, label="CO2-direct")
plt.legend()
plt.show()


# %%------------   PLOT POWER                                           ---------------------------------------------> #
plt.plot(water_dict["m_dot"], water_dict["dh"] * water_dict["m_dot"] / 1e3, label="Water")
plt.plot(CO2_dict["m_dot"], CO2_dict["dh"] * CO2_dict["m_dot"] / 1e3, label="CO2")
plt.legend()
plt.show()

