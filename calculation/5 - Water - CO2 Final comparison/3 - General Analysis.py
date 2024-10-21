# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from main_classes import BaseBHE, ReservoirProperties, BHEGeometry
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np


# %%------------   IMPORT VALIDATION DATA AND INITIALIZE CALCULATION    ---------------------------------------------> #
# Input Parameters
# General Params
t_in = 10                   # [°C]
time = 20                   # [year]
s_year = 365 * 24 * 3600    # [s / year]
flow_rates_ovr = np.linspace(0, 300, 401)[1:]
depths = np.linspace(2000, 5000, 50)    # [m]
grads = np.linspace(10, 50, 50) / 1e3   # [°C/m]
depths, grads = np.meshgrid(depths, grads, indexing='ij')

# Geometric Fixed Params
bhe_geom = BHEGeometry()
bhe_geom.l_horiz = 5000      # [m]
bhe_geom.d_well = 0.5        # [m]
bhe_geom.n_wells = 4

# Reservoir Fixed Params
res_prop = ReservoirProperties()
res_prop.k_rock = 2.1       # [W/(m K)]
res_prop.c_rock = 1000      # [J/(kg K)]
res_prop.rho_rock = 2650    # [kg/m^3]


# %%------------   INTI BHEs                                            ---------------------------------------------> #
bhe_in_co2 = ThermodynamicPoint(["CarbonDioxide"], [1], unit_system="MASS BASE SI")
bhe_in_co2.set_variable("T", t_in + 273.15)
bhe_in_co2.set_variable("Q", 0)

bhe_in_water = ThermodynamicPoint(["Water"], [1], unit_system="MASS BASE SI")
bhe_in_water.set_variable("T", t_in + 273.15)
bhe_in_water.set_variable("P", 1e5)

CO2_bhe = BaseBHE(bhe_in_co2, reservoir_properties=res_prop, geometry=bhe_geom)
water_bhe = BaseBHE(bhe_in_water, reservoir_properties=res_prop, geometry=bhe_geom)


# %%------------   CALCULATE RESULTS                                    ---------------------------------------------> #
pbar = tqdm(desc="Calculating Points", total=grads.shape[0]*grads.shape[1])

res_shape = np.append(grads.shape, 2)
dexs = np.empty(res_shape)
dhs = np.empty(res_shape)
dps = np.empty(res_shape)
dts = np.empty(res_shape)
m_dot_opts = np.empty(res_shape)

for i, sub_grads in enumerate(grads):

    for j, grad in enumerate(sub_grads):

        CO2_bhe.geom.depth = depths[i, j]
        water_bhe.geom.depth = depths[i, j]

        CO2_bhe.res_prop.grad = grad
        water_bhe.res_prop.grad = grad

        CO2_bhe.set_HX_condition()
        CO2_dict = CO2_bhe.evaluate_HXG(times=time * s_year, m_dot=flow_rates_ovr)

        k_max = np.argmax(np.nan_to_num(CO2_dict["dex"] * CO2_dict["m_dot"], nan=-np.inf))
        CO2_dict = CO2_bhe.evaluate_surface_condition(CO2_dict, indices=[k_max])

        dexs[i, j, 0] = CO2_dict["dex"][k_max]
        dhs[i, j, 0] = CO2_dict["dh"][k_max]
        dps[i, j, 0] = CO2_dict["dp_overall"][k_max]
        dts[i, j, 0] = CO2_dict["dt_overall"][k_max]
        m_dot_opts[i, j, 0] = CO2_dict["m_dot"][k_max]

        water_bhe.set_HX_condition()
        water_dict = water_bhe.evaluate_HXG(times=time * s_year, m_dot=flow_rates_ovr)

        k_max = np.argmax(np.nan_to_num(water_dict["dex"] * water_dict["m_dot"], nan=-np.inf))
        water_dict = water_bhe.evaluate_surface_condition(water_dict, indices=[k_max])

        dexs[i, j, 1] = water_dict["dex"][k_max]
        dhs[i, j, 1] = water_dict["dh"][k_max]
        dps[i, j, 1] = water_dict["dp_overall"][k_max]
        dts[i, j, 1] = water_dict["dt_overall"][k_max]
        m_dot_opts[i, j, 1] = water_dict["m_dot"][k_max]

        pbar.update(1)

pbar.close()


# %%------------   PLOT OPTIMAL FLOW RATE                         ---------------------------------------------> #
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

for i, fluid in enumerate(["Carbon Dioxide", "Water"]):

    cont = axs[i].contourf(grads * 1e3, depths / 1e3, m_dot_opts[:, :, i])
    axs[i].set_title("{}".format(fluid))
    axs[i].set_xlabel("Gradient (°C/km)")
    axs[i].set_ylabel("Depth (km)")
    fig.colorbar(cont, ax=axs[i])

fig.suptitle("Optimal Flow Rate (kg/s)", fontsize=14)
plt.tight_layout(pad=1)
plt.show()


# %%------------   PLOT EXERGY OUTPUT                                   ---------------------------------------------> #
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

for i, fluid in enumerate(["Carbon Dioxide", "Water"]):
    cont = axs[i].contourf(grads * 1e3, depths / 1e3, dexs[:, :, i] * m_dot_opts[:, :, i] / 1e3)
    axs[i].set_title("{}".format(fluid))
    axs[i].set_xlabel("Gradient (°C/km)")
    axs[i].set_ylabel("Depth (km)")
    fig.colorbar(cont, ax=axs[i])

fig.suptitle("Optimal Exergy Output (kW)", fontsize=14)
plt.tight_layout(pad=1)
plt.show()


# %%------------   PLOT ENERGY OUTPUT                                   ---------------------------------------------> #
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

for i, fluid in enumerate(["Carbon Dioxide", "Water"]):
    cont = axs[i].contourf(grads * 1e3, depths / 1e3, dhs[:, :, i] * m_dot_opts[:, :, i] / 1e3)
    axs[i].set_title("{}".format(fluid))
    axs[i].set_xlabel("Gradient (°C/km)")
    axs[i].set_ylabel("Depth (km)")
    fig.colorbar(cont, ax=axs[i])

fig.suptitle("Optimal Energy Output (kW)", fontsize=14)
plt.tight_layout(pad=1)
plt.show()


# %%------------   PLOT COMPARISON                                      ---------------------------------------------> #
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

delta_hs = dhs[:, :, 0] * m_dot_opts[:, :, 0] - dhs[:, :, 1] * m_dot_opts[:, :, 1]
cont1 = ax1.contourf(grads * 1e3, depths / 1e3, delta_hs / 1e3)
ax1.set_title("Energy Output Difference (kW)")
fig.colorbar(cont1, ax=ax1)
ax1.set_xlabel("Gradient (°C/km)")
ax1.set_ylabel("Depth (km)")

delta_exs = dexs[:, :, 0] * m_dot_opts[:, :, 0] - dexs[:, :, 1] * m_dot_opts[:, :, 1]
cont2 = ax2.contourf(grads * 1e3, depths / 1e3, delta_exs / 1e3)
ax2.set_title("Exergy Output Difference (kW)")
fig.colorbar(cont2, ax=ax2)
ax2.set_xlabel("Gradient (°C/km)")
ax2.set_ylabel("Depth (km)")

fig.suptitle("Comparison Between CO2 and Water (kW)", fontsize=14)
plt.tight_layout(pad=1)
plt.show()
