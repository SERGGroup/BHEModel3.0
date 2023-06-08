# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from itertools import cycle
import scipy.constants
from tqdm import tqdm
import numpy as np

# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
n_points = 250
fluids = [

    "water", "CarbonDioxide", "CarbonMonoxide", "Nitrogen", "Oxygen",
    "methane", "Ethane", "Pentane", "Isobutane", "Propane",
    "Ammonia", "Nitrogen", "Helium", "Hydrogen", "Neon",
    "R134a", "MM", "R1234ze(E)"

]
R = scipy.constants.R

b1 = np.power(2, 1/3) - 1
vs = np.logspace(-2, 5, n_points) + b1
ts = np.array([0.5, 0.7, 0.8, 0.9, 1, 1.1, 1.5, 2, 5, 10, 50, 200])
v_rels, t_rels = np.meshgrid(vs, ts, indexing='ij')

p_rels = np.empty((n_points, len(ts), len(fluids)))
cvs = np.empty((n_points, len(ts), len(fluids)))
dcps = np.empty((n_points, len(ts), len(fluids)))
z_rels = np.empty((n_points, len(ts), len(fluids)))


# %%-------------------------------------   CALCULATE REFPROP                   -------------------------------------> #
pbar = tqdm(desc="Calculating Points", total=n_points * len(ts) * len(fluids))

for i in range(len(fluids)):

    fluid = fluids[i]
    tp = ThermodynamicPoint([fluid], [1], unit_system="MOLAR BASE SI")
    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    tp.set_variable("T", t_crit)
    tp.set_variable("P", p_crit)
    v_crit = 1 / tp.get_variable("rho")

    t_reals = t_rels * t_crit
    rho_reals = 1 / (v_rels * v_crit)

    for j in range(n_points):

        for k in range(len(ts)):

            try:

                tp.set_variable("T", t_reals[j, k])
                tp.set_variable("rho", rho_reals[j, k])

                p_real = tp.get_variable("P")
                p_rel = p_real / p_crit
                cv = tp.get_variable("CV")
                cp = tp.get_variable("CP")

            except Exception as a:

                p_rel = -1
                cv = -1
                cp = -1

            if p_rel > 0:
                p_rels[j, k, i] = p_rel

            if cv > 0:
                cvs[j, k, i] = cv / (3/2 * R)

            if cp > 0:
                dcps[j, k, i] = (cp - cv) * t_crit / (v_crit * p_crit) * 0.8

            pbar.update(1)

    z_rels[:,:,i] = p_rels[:,:,i] * v_rels / t_rels

pbar.close()


# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #
h_fig = 10
n_sub_plots = 2
fig, axis = plt.subplots(1, n_sub_plots)
fig.set_size_inches(n_sub_plots * h_fig, h_fig)
plot_every = 2

colors = pl.cm.viridis(np.linspace(0,1,len(ts)+1))

for i in range(len(ts)):

    markers = cycle(["^", "+", "x", ".", "o", "s"])

    for j in range(len(fluids)):

        marker = next(markers)
        sub_x = v_rels[:, i] #- b1
        sub_y = p_rels[:, i, j]
        axis[0].plot(sub_x[::plot_every], sub_y[::plot_every], color=colors[i])#, marker=marker, s=10)

        sub_x = p_rels[:, i, j] #- b1
        sub_y = z_rels[:, i, j]
        axis[1].plot(sub_x[::plot_every], sub_y[::plot_every], color=colors[i])#, marker=marker, s=10)

    for ax in axis:

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("$v_{rel}$ [-]", size="xx-large")

axis[0].set_ylim((0.0001, 100000))
axis[1].set_ylim((0.0001, 100))
plt.show()


