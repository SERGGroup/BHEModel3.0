# %%------------   IMPORT MODULES                         -----------------------------------------------------------> #
from main_classes.support.error_band import draw_error_band
from main_classes.constant import CALCULATION_DIR
from REFPROPConnector import ThermodynamicPoint
from main_classes import BaseBHE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os, math


# %%------------   IMPORT VALIDATION DATA                 -----------------------------------------------------------> #

#   Validation data from:
#       "20221107 HOCLOOP - CO2 case well" (in res folder)
#
#   Data extracted from plot using:
#       WebPlotDigitizer (https://automeris.io/WebPlotDigitizer)

RES_FOLDER = os.path.join(CALCULATION_DIR, "4 - Heat Exchanger Based Model", "0 - res")
VALIDATION_FILE = os.path.join(RES_FOLDER, "validation_datasets", "t_over_time_validation_datasets.csv")
VALIDATION_ARR = np.array(pd.read_csv(VALIDATION_FILE).drop(0), dtype=float)
VALIDATION_DICT = {}
keys = [10, 5, 1]

for i in range(int(VALIDATION_ARR.shape[1] / 2)):

    j = 0

    for j in range(VALIDATION_ARR.shape[0]):

        if math.isnan(VALIDATION_ARR[j, 2 * i]):

            j -= 1
            break

    VALIDATION_DICT.update({

        keys[i]: {

            "t [years]": VALIDATION_ARR[:j+1, 2 * i],
            "T [°C]": VALIDATION_ARR[:j+1, 2 * i + 1]

        }

    })


# %%------------   INPUT DATA DEFINITION                  -----------------------------------------------------------> #
t_in = 30           # [C]
depth = 4500        # [m]
l_overall = 10000   # [m]
r_cur = 500         # [m]

# !!! NOTE !!!
#
#   In the document provided by Ola (in the res folder) t_rock = 172°C while from the graphs used for the
#   validation (at pg. 2) such temperature appears to be around 160°C. This discrepancy, make me doubts about
#   the consistency of the other data.
#
#   For this reason, a correction factor has been introduced for the flow rate to overlap the calculated curve
#   and the validation points. The flow rate correction has been introduced only because it was easier to
#   implement, it simply reflects the fact that other parameters (such as rock conductivity or density) may
#   be different.
#
#   This has to be clarified with Ola.

t_rock = 160        # [C]
k_rock = 2.44       # [W/(m K)]
c_rock = 1.085      # [kJ/(kg K)]
rho_rock = 2542     # [kg/m^3]

l_horiz = l_overall - depth
gradient = (t_rock - t_in) / depth


# %%------------   INIT WELL CALCULATION                  -----------------------------------------------------------> #
bhe_in = ThermodynamicPoint(["Water"], [1], unit_system="MASS BASE SI")
bhe_in.set_variable("T", t_in + 273.15)
bhe_in.set_variable("P", 2e5)

base_bhe = BaseBHE(bhe_in)
base_bhe.geom.depth = depth

base_bhe.res_prop.grad = gradient
base_bhe.res_prop.c_rock = c_rock * 1e3
base_bhe.res_prop.rho_rock = rho_rock
base_bhe.res_prop.k_rock = k_rock

base_bhe.set_HX_condition()


# %%-------------------------------------   ANALYZE SYSTEM                      -------------------------------------> #
RESULTS_DICT = {}
n_points = 30
shape = 2.5
time_points = np.power((np.array(range(n_points)) + 1) / n_points, shape) * 10

for key in VALIDATION_DICT.keys():

    time_list = list()
    t_out_list = list()
    m_dot = key / 1.55

    dh_percs, output_points, drho_down = base_bhe.evaluate_HXG(time_points * 3.154e4, m_dot=key)

    i = 0
    for output_point in output_points:

        time_list.append(time_points[i])
        t_out_list.append(output_point.get_variable("T")-273.15)

        i += 1

    RESULTS_DICT.update({

        key: {

            "t [years]": time_list,
            "T [°C]": t_out_list

        }

    })

# %%------------   PLOT RESULTS                           -----------------------------------------------------------> #
fig = plt.figure(dpi=150)
fig.set_size_inches(10, 8)
ax = fig.add_subplot(1, 1, 1)

lines = list()
labels = list()

for key in RESULTS_DICT.keys():

    x_calc = np.array(RESULTS_DICT[key]["t [years]"])
    y_calc = np.array(RESULTS_DICT[key]["T [°C]"])

    x_val = np.array(VALIDATION_DICT[key]["t [years]"])
    y_val = np.array(VALIDATION_DICT[key]["T [°C]"])

    line = draw_error_band(ax, x_calc, y_calc, 1.5, err_fixed=True, alpha=0.4, zorder=5)
    ax.scatter(x_val, y_val, s=80, marker="x", zorder=10, linewidths=1.5, color=line.get_color())

    lines.append(line)
    labels.append("$m_{{dot}}$ = {}[kg/s]".format(key))

ax.set_xlabel(r'time [years]', fontsize='large', loc='center')
ax.set_ylabel(r'$T_{out}$ [°C]', fontsize='large', loc='center')
plt.legend(lines, labels)
plt.show()


# %%------------   SAVE IMAGE                             -----------------------------------------------------------> #
output_directory = os.path.join(RES_FOLDER, "output")

if not os.path.isdir(output_directory):

    os.mkdir(output_directory)

fig.savefig(os.path.join(output_directory, "2 - Model Validation.png"))