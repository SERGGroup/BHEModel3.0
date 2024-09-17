# %%-------------------------------------   IMPORT MODULES                          ---------------------------------> #
from REFPROPConnector import ThermodynamicPoint, DiagramPlotter, DiagramPlotterOptions
from main_classes import BaseBHE
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


# %%-------------------------------------   INIT CALCULATION                        ---------------------------------> #
w_dt_max = 5e5      # [W]
t_in = 0           # [°C]
depth = 500         # [m]
l_horiz = 3500      # [m]
main_fluid = "Ethane"
other_fluid = "Carbon Dioxide"
times = np.array([1 / 365, 10])

mix_conc = np.array([0, 1])
gradients = np.linspace(25, 75, 20) / 1e3

mix_conc, gradients = np.meshgrid(mix_conc, gradients)


# %%-------------------------------------   INIT BHE                                ---------------------------------> #
bhe_calculators = list()
p_ins = np.zeros(mix_conc.shape)
pbar = tqdm(desc="Calculating Points", total=mix_conc.shape[0]*mix_conc.shape[1])
for i in range(mix_conc.shape[0]):

    for j in range(mix_conc.shape[1]):

        conc = mix_conc[i, j]
        if conc < 0.5:

            fluid = main_fluid

        else:

            fluid = other_fluid

        bhe_in_mix = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
        bhe_in_mix.set_variable("T", t_in + 273.15)
        bhe_in_mix.set_variable("Q", 0)
        p_ins[i, j] = bhe_in_mix.get_variable("P")

        mix_bhe = BaseBHE(bhe_in_mix)
        mix_bhe.geom.depth = depth
        mix_bhe.geom.l_horiz = l_horiz
        mix_bhe.res_prop.grad = gradients[i, j]
        mix_bhe.set_HX_condition()
        bhe_calculators.append(mix_bhe)
        pbar.update(1)


pbar.close()


# %%-------------------------------------   RETRIEVE RESULTS                        ---------------------------------> #
pbar = tqdm(desc="Calculating Points", total=len(times)*mix_conc.shape[0]*mix_conc.shape[1])
results = np.empty((4, len(times) + 1, gradients.shape[0], gradients.shape[1]))
results[:] = np.nan

for i in range(mix_conc.shape[0]):

    for j in range(mix_conc.shape[1]):

        base_bhe = bhe_calculators[j + i * mix_conc.shape[1]]
        m_dot_curr = w_dt_max / base_bhe.integrator.dh_max
        dh_percs, output_points, drhos = base_bhe.evaluate_HXG(times * 3.154e+7, m_dot=m_dot_curr)

        results[0, 0, i, j] = base_bhe.integrator.dh_max
        results[3, 0, i, j] = base_bhe.ideal_exergy_efficiency

        results[0, 1:, i, j] = dh_percs * base_bhe.integrator.dh_max

        point_ideal = base_bhe.ideal_points[-1]
        results[1, 0, i, j] = point_ideal.get_variable("T") - (t_in + 273.15)
        results[2, 0, i, j] = point_ideal.get_variable("P") - p_ins[i, j]

        n = 1
        for point in output_points:

            results[1, n, i, j] = point.get_variable("T") - (t_in + 273.15)
            results[2, n, i, j] = point.get_variable("P") - p_ins[i, j]
            pbar.update(1)
            n += 1

pbar.close()


# %%-------------------------------------   PLOT IDEAL PRESSURE INCREASE        -------------------------------------> #
time_i = 0
dPs = results[2, time_i, :, :] + p_ins[:, :]
dPs[np.where(dPs < 0)] = np.nan
dPs = dPs - p_ins
betas = dPs / p_ins + 1
betas_corr = (dPs - 2e6) / p_ins + 1
dPs_MPA = dPs / 1e6
# dPs_MPA[3:4, 1] = np.nan
#betas[3, 1] = np.nan
#betas[4, 1] = np.nan
fluids = [main_fluid, other_fluid]

for j in [0, 1]:

    plt.plot(

        gradients[:, j] * depth + t_in, betas[:, j],
        label="{}".format(fluids[j])

    )

plt.legend()
plt.xlabel("Temperature @ {}m (°C)".format(depth))
plt.ylabel("Beta (-)")
plt.title("{} - {} Comparison".format(main_fluid, other_fluid))
plt.show()


# %%-------------------------------------   PLOT POWER AND EXERGY OUTPUT        -------------------------------------> #
time_i = 0
dh_max = results[3, time_i, :, :] #/ 1e3
dPs = results[2, time_i, :, :] + p_ins[:, :]
dh_max[np.where(dPs < 0)] = np.nan

for j in [0, 1]:

    plt.plot(

        gradients[:, j] * depth + t_in, dh_max[:, j],
        label="{}".format(fluids[j])

    )

plt.legend()
plt.xlabel("Temperature @ {}m (°C)".format(depth))
plt.ylabel("rel power output (kJ/kg)")
plt.title("{} - {} Comparison".format(main_fluid, other_fluid))
plt.show()