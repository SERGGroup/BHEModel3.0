# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import BaseBHE, economicEvaluator
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import numpy as np


# %%-------------------------------------   INIT CALCULATION                    -------------------------------------> #
t_in = 10           # [°C]
depth = 1000         # [m]
gradient = 0.05     # [°C/m]

bhe_in = ThermodynamicPoint(["CarbonDioxide"], [1], unit_system="MASS BASE SI")
bhe_in.set_variable("T", t_in + 273.15)
bhe_in.set_variable("Q", 0)

base_bhe = BaseBHE(bhe_in)
base_bhe.geom.depth = depth
base_bhe.res_prop.grad = gradient

base_bhe.set_HX_condition()

param = base_bhe.integrator.params
h_rel = param[0][:-1]
integral = param[2][:-1]
plt.scatter(h_rel, integral)

base_bhe.integrator.solve_analytically = True
int_2 = base_bhe.integrator.evaluate_integral(h_rel)
plt.plot(h_rel, int_2)
plt.show()


# %%-------------------------------------   ANALYZE SYSTEM                      -------------------------------------> #
times = np.logspace(-4, 1, 100)
base_bhe.integrator.solve_analytically = False

for Pe in [0, 0.1, 1, 10]:

    base_bhe.res_prop.pe = Pe

    dh_percs, output_points = base_bhe.evaluate_HXG(times * 3.154e+7, m_dot=10)
    plt.plot(times, dh_percs*base_bhe.integrator.dh_max*10/1e6, label="Pe={}".format(Pe))

plt.legend()
plt.xscale("log")
plt.xlabel("Time (year)")
plt.ylabel("power (MW)")
plt.show()
