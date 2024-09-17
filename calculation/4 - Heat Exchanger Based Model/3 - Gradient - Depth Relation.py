# %%-------------------------------------   IMPORT MODULES                          ---------------------------------> #
from main_classes import BaseBHE, economicEvaluator
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


# %%-------------------------------------   INIT CALCULATION                        ---------------------------------> #
m_dot = 1
t_in = 10           # [°C]
depth = 3000         # [m]
l_horiz = 3500      # [m]
main_fluid = "Carbon Dioxide"
other_fluid = "Propane"
times = np.array([1 / 365, 10])

mix_conc = np.array([0, 0.1, 0.25, 0.5, 0.75, 1])
gradients = np.linspace(25, 75, 20) / 1e3

mix_conc, gradients = np.meshgrid(mix_conc, gradients)


# %%-------------------------------------   INIT BHE                                ---------------------------------> #
bhe_calculators = list()
p_ins = np.zeros(mix_conc.shape)
pbar = tqdm(desc="Calculating Points", total=mix_conc.shape[0]*mix_conc.shape[1])
for i in range(mix_conc.shape[0]):

    for j in range(mix_conc.shape[1]):

        conc = mix_conc[i, j]
        bhe_in_mix = ThermodynamicPoint([main_fluid, other_fluid], [1 - conc, conc], unit_system="MASS BASE SI")
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
results = np.empty((3, len(times), gradients.shape[0], gradients.shape[1]))
results[:] = np.nan

for i in range(mix_conc.shape[0]):

    for j in range(mix_conc.shape[1]):

        base_bhe = bhe_calculators[j + i * mix_conc.shape[1]]
        dh_percs, output_points, drhos = base_bhe.evaluate_HXG(times * 3.154e+7, m_dot=m_dot)

        results[0, :, i, j] = dh_percs

        n = 0
        for point in output_points:

            results[1, n, i, j] = point.get_variable("T") - (t_in + 273.15)
            results[2, n, i, j] = point.get_variable("P") - p_ins[i, j]
            pbar.update(1)
            n += 1

pbar.close()


# %%-------------------------------------   PLOT PRESSURE INCREASE              -------------------------------------> #
dPs = results[2, 0, :, :] + p_ins[:, :]
dPs[np.where(dPs < 0)] = np.nan
dPs = dPs - p_ins
betas = dPs / p_ins + 1
dPs_MPA = dPs / 1e6

for j in [0, 1, 2, 3, 4]:

    plt.plot(

        gradients[:, j]*1e3, betas[:, j],
        label="{}={}$\%_{{mass}}$".format(

            other_fluid,
            int(mix_conc[0, j]*100)

        )

    )

plt.legend()
plt.xlabel("Gradient (°C/km)")
plt.ylabel("Beta (-)")
plt.title("$CO_2$ - {} Mixture".format(other_fluid))
plt.show()


# %%-------------------------------------   PLOT INTEGRAL PROFILE               -------------------------------------> #
base_bhes = [

    bhe_calculators[0 + 10 * mix_conc.shape[1]],
    bhe_calculators[2 + 10 * mix_conc.shape[1]],
    bhe_calculators[3 + 10 * mix_conc.shape[1]]

]

labels = [

    "Pure $CO_2$",
    "{}={}$\%_{{mass}}$".format(

        other_fluid,
        int(mix_conc[0, 2]*100)

    ),
    "Pure Ethane"

]

i = 0
for base_bhe in base_bhes:

    param = base_bhe.integrator.params
    h_rel = param[0][:-1]
    integral = param[2][:-1]

    plt.plot(h_rel, 1/integral, label=labels[i])
    i += 1


d = base_bhe.geom.d_well
UAs = base_bhe.res_prop.evaluate_rel_resistance(times=times * 3.154e+6, d=d) * np.pi * d

colors = ["lightgray", "gray"]
linestyles = ["--", ":"]

i = 0
for length in [350, 1000]:

    j = 0
    for UA in UAs:

        DT_mean = np.ones(integral.shape) * UA / length * (m_dot * base_bhe.integrator.dh_max)
        plt.plot(h_rel, DT_mean, linestyles[j], color=colors[i])
        j += 1

    i += 1

plt.yscale("log")
plt.xlim((0.31, 1))
plt.legend()

plt.xlabel("% of Heat Exchange")
plt.ylabel("$DT_{mean}$ (°C)")
plt.title("$CO_2$ - {} Mixture".format(other_fluid))
plt.show()
