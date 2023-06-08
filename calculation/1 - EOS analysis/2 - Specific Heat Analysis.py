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
#fluids = ["ammonia"]
fluids = ["water", "carbonDioxide"]
R = scipy.constants.R
f = 3

# %%-------------------------------------   CALCULATE VAN DER WALLS             -------------------------------------> #
v_mods = np.logspace(np.log10(0.4), np.log10(9), n_points)
t_mods = np.array([0.9, 1, 1.1, 1.5, 2, 10])

p_rels = np.empty((n_points, len(t_mods), len(fluids) + 1))
cvs = np.empty((n_points, len(t_mods), len(fluids) + 1))
cps = np.empty((n_points, len(t_mods), len(fluids) + 1))

p_rels[:, :, :] = np.nan
cvs[:, :, :] = np.nan
cps[:, :, :] = np.nan

xv, yv = np.meshgrid(v_mods, t_mods, indexing='ij')

alpha = yv * xv ** 3
beta = (3 * xv - 1) ** 2
cvs[:, :, 0] = np.ones(xv.shape) * f / 2 * R
cps[:, :, 0] = cvs[:, :, 0] + 4 * R * alpha / (4 * alpha - beta)
p_rels[:, :, 0] = 8 * yv / (3 * xv - 1) - 3 / (xv ** 2)
cps[cps<0] = np.inf


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

            except Exception:

                p_rel = -1
                cv = -1
                cp = -1

            if p_rel > 0:
                p_rels[j, k, i + 1] = p_rel

            if cv > 0:
                cvs[j, k, i + 1] = cv

            if cp > 0:
                cps[j, k, i + 1] = cp

            pbar.update(1)

pbar.close()

# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
h_fig = 10
n_sub_plots = 2
fig, axis = plt.subplots(1, n_sub_plots)
fig.set_size_inches(n_sub_plots * h_fig, h_fig)
plot_every = 15

colors = pl.cm.viridis(np.linspace(0,1,len(t_mods)+1))
markers = cycle(["^", "+", "x", ".", "o", "s"])

for i in range(len(t_mods)):

    k = len(t_mods) - i - 1
    axis[0].plot(xv[:, k], p_rels[:, k, 0], color=colors[k])
    axis[1].plot(xv[:, k], cps[:, k, 0] - cvs[:, k, 0], color=colors[k])

    markers = cycle(["^", "+", "x", ".", "o", "s"])
    for j in range(len(fluids)):

        marker = next(markers)
        sub_x = xv[:, i]
        sub_y = p_rels[:, i, j + 1]
        axis[0].scatter(sub_x[::plot_every], sub_y[::plot_every], color=colors[i], marker=marker, s=10)

        sub_y = cps[:, i, j + 1] - cvs[:, i, j + 1]
        axis[1].scatter(sub_x[::plot_every], sub_y[::plot_every], color=colors[i], marker=marker, s=10)

    for ax in axis:

        ax.set_xscale("log")
        ax.set_yscale("log")

plt.ylim((1, 10000))
plt.show()
