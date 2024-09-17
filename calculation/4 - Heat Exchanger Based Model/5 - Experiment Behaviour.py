# %%-------------------------------------   IMPORT MODULES                          ---------------------------------> #
from REFPROPConnector import ThermodynamicPoint, DiagramPlotter, DiagramPlotterOptions
from main_classes import BaseBHE
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import g


# %%-------------------------------------   INPUT DATA                              ---------------------------------> #
w_dt_max = 5e5      # [W]
t_in = 20           # [°C]
depth = 2.5         # [m]

dTs = np.linspace(0, 40, 40)

# %%-------------------------------------   EVALUATE                                ---------------------------------> #
fluids = ["Ethane", "Carbon Dioxide"]
p_ins = [5.5, 10]
ovr_results = list()
rel_results = list()
i = 0
for fluid in fluids:

    p_in = p_ins[i]
    i += 1

    input_point = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
    input_point.set_variable("T", t_in + 273.15)
    input_point.set_variable("P", p_in * 1e6)

    down_point = input_point.duplicate()
    down_point.set_variable("H", input_point.get_variable("H") + g * depth)
    down_point.set_variable("S", input_point.get_variable("S"))

    tmp_points = [down_point.duplicate(), down_point.duplicate()]

    result = list()
    for dT in dTs:

        tmp_points[0].set_variable("T", input_point.get_variable("T") + dT)
        tmp_points[0].set_variable("P", down_point.get_variable("P"))

        tmp_points[1].set_variable("H", tmp_points[0].get_variable("H") - g * depth)
        tmp_points[1].set_variable("S", tmp_points[0].get_variable("S"))

        result.append(tmp_points[1].get_variable("P") - input_point.get_variable("P"))

    result = np.array(result) / 1e3
    ovr_results.append(result)
    rel_results.append(result/np.max(result))

ovr_results = np.array(ovr_results)
rel_results = np.array(rel_results)

# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
for j in [0, 1]:

    plt.plot(

        dTs, rel_results[j, :],
        label="{}".format(fluids[j])

    )

plt.legend()
plt.xlabel("$\Delta T$ (°C)")
plt.ylabel("% max flow rate")
plt.title("{} - {} Comparison".format(fluids[0], fluids[1]))
plt.show()