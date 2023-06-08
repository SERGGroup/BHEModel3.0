# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from itertools import cycle
import scipy.constants
from tqdm import tqdm
import numpy as np

# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
n_points = 200
fluids = ["methane"]
# fluids = ["water", "carbonDioxide"]
# fluids = ["water", "carbonDioxide", "methane", "ammonia", "n-pentane"]
R = scipy.constants.R


# %%-------------------------------------   CALCULATE REDLICH-KONG              -------------------------------------> #
b1 = np.power(2, 1/3) - 1
v_mods = np.logspace(-2, 2, n_points) + b1
t_mods = np.array([0.9, 1, 1.1, 1.5, 2, 10])

p_rels = np.empty((n_points, len(t_mods), len(fluids) + 1))
p_mods = np.empty((n_points, len(t_mods), len(fluids) + 1))
cvs = np.empty((n_points, len(t_mods), len(fluids) + 1))
dcps = np.empty((n_points, len(t_mods), len(fluids) + 1))
dhdps = np.empty((n_points, len(t_mods), len(fluids) + 1))

p_rels[:, :, :] = np.nan
p_mods[:, :, :] = np.nan

xv, yv = np.meshgrid(v_mods, t_mods, indexing='ij')

alpha = (xv - b1)
beta = b1 * np.power(yv, 1/2) * xv
gamma = xv + b1
p_rels[:, :, 0] = (3 * yv / alpha - 1 /(beta * gamma))

# delta = 1 / ((xv - b1)**2 * (yv ** -2))
# corr = (np.power(delta/10000 + 1, 1/2) + 1) / 2
# p_mods[:, :, 0] = p_rels[:, :, 0] / corr

dp_dt = 3 / alpha + 1 / (2 * beta * gamma * yv)
dp_dv = - 3 * yv / alpha ** 2 + (2 * xv + b1) / (beta * gamma ** 2 * xv)
dcps[:, :, 0] = - yv * dp_dt ** 2 / dp_dv
dhdps[:, :, 0] = 1 / 3 * (xv + yv * dp_dt / dp_dv)
cvs[:, :, 0] = 1
dcps[dcps < 0] = np.inf


# %%-------------------------------------   CALCULATE REFPROP                   -------------------------------------> #
pbar = tqdm(desc="Calculating Points", total=n_points * len(t_mods) * len(fluids))

for i in range(len(fluids)):

    fluid = fluids[i]
    tp = ThermodynamicPoint([fluid], [1], unit_system="MOLAR BASE SI")
    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    tp.set_variable("T", t_crit)
    tp.set_variable("P", p_crit)
    v_crit = 1 / tp.get_variable("rho")

    for j in range(n_points):

        for k in range(len(t_mods)):

            t_real = yv[j, k] * t_crit
            rho_real = 1 / (xv[j, k] * v_crit)

            try:

                tp.set_variable("T", t_real)
                tp.set_variable("rho", rho_real)

                p_real = tp.get_variable("P")
                p_rel = p_real / p_crit
                cv = tp.get_variable("CV")
                cp = tp.get_variable("CP")
                dh_real = tp.get_derivative("H", "P", "T")

            except Exception:

                p_rel = -1
                cv = -1
                cp = -1
                dh_real = -1

            if p_rel > 0:
                p_rels[j, k, i + 1] = p_rel
                p_mods[j, k, i + 1] = p_rel

            if cv > 0:
                cvs[j, k, i + 1] = cv / (3/2 * R)

            if cp > 0:
                dcps[j, k, i + 1] = (cp - cv) * t_crit / (v_crit * p_crit) * 0.8

            if dh_real > 0:
                dhdps[j, k, i + 1] = dh_real / (R * t_crit)

            pbar.update(1)

pbar.close()

# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
h_fig = 10
n_sub_plots = 2
fig, axis = plt.subplots(1, n_sub_plots)
fig.set_size_inches(n_sub_plots * h_fig, h_fig)
plot_every = 15

colors = pl.cm.viridis(np.linspace(0,1,len(t_mods)+1))

for i in range(len(t_mods)):

    k = len(t_mods) - i - 1
    axis[0].plot(xv[:, k] - b1, p_rels[:, k, 0], color=colors[k])
    axis[1].plot(xv[:, k] - b1, dhdps[:, k, 0], color=colors[k])

    axis[0].set_ylabel("$p_{rel}$ [-]", size="xx-large")
    axis[1].set_ylabel("$c_{p_{rel}}$ - $c_{v_{rel}}$ [-]", size="xx-large")

    markers = cycle(["^", "+", "x", ".", "o", "s"])

    for j in range(len(fluids)):

        marker = next(markers)
        sub_x = xv[:, i] - b1
        sub_y = p_rels[:, i, j + 1]
        axis[0].scatter(sub_x[::plot_every], sub_y[::plot_every], color=colors[i], marker=marker, s=10)

        sub_y = dhdps[:, i, j + 1] * 5E6
        axis[1].scatter(sub_x[::plot_every], sub_y[::plot_every], color=colors[i], marker=marker, s=10)

    for ax in axis:

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("$v_{rel}$ [-]", size="xx-large")

plt.ylim((0.001, 100))
plt.show()

# %%-------------------------------------   PLOT CVS                            -------------------------------------> #
h_fig = 10
n_sub_plots = 1
fig, ax = plt.subplots(1, n_sub_plots)
fig.set_size_inches(n_sub_plots * h_fig, h_fig)
plot_every = 1

colors = pl.cm.viridis(np.linspace(0,1,len(t_mods)+1))

for i in range(len(t_mods)):

    k = len(t_mods) - i - 1
    ax.plot(yv[:, k], cvs[:, k, 0], color=colors[k])

    markers = cycle(["^", "+", "x", ".", "o", "s"])

    for j in range(len(fluids)):

        marker = next(markers)
        sub_x = yv[:, i]
        sub_y = cvs[:, i, j + 1]
        ax.scatter(sub_x[::plot_every], sub_y[::plot_every], color=colors[i], marker=marker, s=10)

ax.set_ylabel("$c_{v_{rel}}$ [-]", size="xx-large")
ax.set_xscale("log")
#ax.set_yscale("log")
ax.set_xlabel("$T_{rel}$ [-]", size="xx-large")

plt.ylim((0, 10))
plt.show()