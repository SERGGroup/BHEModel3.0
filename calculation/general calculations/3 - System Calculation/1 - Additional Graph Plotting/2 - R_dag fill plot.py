# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import VdWEOS, RKEOS, SRKEOS, PREOS, convert_cp_coeff
from main_classes.constant import CALCULATION_DIR
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
from itertools import cycle
from tqdm import tqdm
import numpy as np
import os


# %%-------------------------------------   CALCULATIONS OPTIONS                -------------------------------------> #
n_points = 1000
t_rels = np.concatenate(

    (np.logspace(0, 1, 150),
     np.logspace(np.log10(0.2), 0, 100)[:-1]),
    axis=None

)

t_rels.sort()
v_rels = np.logspace(-3, 3, n_points)
p_rels = np.logspace(-2, 4, n_points)

t_tv_mesh, v_tv_mesh = np.meshgrid(t_rels, v_rels, indexing='ij')
t_tp_mesh, p_tp_mesh = np.meshgrid(t_rels, p_rels, indexing='ij')


# %%-------------------------------------   INIT FLUID                          -------------------------------------> #
fluid = "Methane"

tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
tp.set_variable("T", tp.RPHandler.TC)
tp.set_variable("P", tp.RPHandler.PC)
cp0_crit = tp.get_variable("CP0")

t_crit = tp.RPHandler.TC
p_crit = tp.RPHandler.PC

acntr = 0.011
m_mol = 0.01604
base_coeff = [0.81558724]
cp_coeffs = convert_cp_coeff(base_coeff, cp0_crit, t_crit)
fluid_curr = PREOS(

    t_crit=t_crit, p_crit=p_crit,
    cp_ideal=cp_coeffs, m_molar=m_mol,
    acntr=acntr

)

crit_state = fluid_curr.get_state(t=t_crit, p=p_crit)


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
k = 0
res_tp_mesh = np.empty(p_tp_mesh.shape)
res_tv_mesh = np.empty(v_tv_mesh.shape)
res_tp_mesh[:] = np.nan
res_tv_mesh[:] = np.nan
p_sat_arr = np.zeros((len(t_rels), 1))
v_sat_arr = np.zeros((len(t_rels), 2))

pbar = tqdm(desc="Calculating Points", total=n_points * len(t_rels) * 2)

for i in range(len(t_rels)):

    t_curr = crit_state.t * t_rels[i]
    sat_states = fluid_curr.get_sat_state(t=t_curr)
    p_sat_arr[i, 0] = sat_states[0].p
    v_sat_arr[i, 0] = sat_states[0].v
    v_sat_arr[i, 1] = sat_states[1].v

    for j in range(len(v_rels)):

        v_curr = v_rels[j] * (crit_state.v - fluid_curr.b) + fluid_curr.b
        curr_state = fluid_curr.get_state(v=v_curr, t=t_curr, sat_states=sat_states)

        if not curr_state.bifase:
            res_tv_mesh[i, j] = curr_state.r / curr_state.cp

        pbar.update(1)

    for k in range(len(p_rels)):

        p_curr = p_rels[k] * crit_state.p
        curr_state = fluid_curr.get_state(p=p_curr, t=t_curr, sat_states=sat_states)

        if not curr_state.bifase:
            res_tp_mesh[i, k] = curr_state.r / curr_state.cp

        pbar.update(1)


p_sat_arr = p_sat_arr / fluid_curr.p_crit
v_sat_arr = (v_sat_arr - fluid_curr.b) / (fluid_curr.v_crit - fluid_curr.b)

pbar.close()

# %%-------------------------------------   IDENTIFY MAX VALUES                 -------------------------------------> #
v_max = [list(), list()]
p_max = [list(), list()]

for i in range(len(t_rels)):

    if t_rels[i] > 1:

        max_r_v_values = np.max(res_tv_mesh[i, :])
        j = np.where(res_tv_mesh[i, :] == max_r_v_values)
        v_max[0].append(t_tv_mesh[i, j])
        v_max[1].append(v_tv_mesh[i, j])

        max_r_t_values = np.max(res_tp_mesh[i, :])
        k = np.where(res_tp_mesh[i, :] == max_r_t_values)
        p_max[0].append(t_tp_mesh[i, k])
        p_max[1].append(p_tp_mesh[i, k])

v_max = np.array(v_max)[:, :, 0, 0]
p_max = np.array(p_max)[:, :, 0, 0]

# %%-------------------------------------   PLOT CONTOUR                        -------------------------------------> #
y_labels = ["$R^{\\dagger}\\ /\\ cp$ [-]"]
fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

n_level = 21
n_ticks_level = 11
n_isoline_level = 6
isoline_width = 0.75
isoline_alpha = 0.6
plot_contour_lines = False

levels = np.linspace(0, 1, n_level)
ticks_levels = np.linspace(0,1, n_ticks_level)
isoline_levels = np.linspace(0,1, n_isoline_level)

axs[0].contourf(p_tp_mesh, t_tp_mesh, res_tp_mesh, levels=levels, zorder=1)
cs = axs[1].contourf(v_tv_mesh, t_tv_mesh, res_tv_mesh, levels=levels, zorder=1)

cbar = fig.colorbar(cs)
cbar.set_ticks(ticks_levels)
cbar.ax.set_ylabel("$R^{\\dagger}/cp$ [-]")

if plot_contour_lines:

    cl_p = axs[0].contour(

        p_tp_mesh, t_tp_mesh, res_tp_mesh,
        levels=isoline_levels, colors="black",
        alpha=isoline_alpha, linewidths=isoline_width,
        zorder=1

    )
    cl_v = axs[1].contour(

        v_tv_mesh, t_tp_mesh, res_tv_mesh,
        levels=isoline_levels, colors="black",
        alpha=isoline_alpha, linewidths=isoline_width,
        zorder=1

    )
    labels_p = axs[0].clabel(cl_p, inline=0.1, fontsize=7, zorder=1)
    labels_v = axs[1].clabel(cl_v, inline=0.1, fontsize=7, zorder=1)

    for l in labels_p + labels_v:
        l.set_rotation(0)

axs[0].plot(p_sat_arr[:, 0], t_rels, color="black", linewidth=2, zorder=2)
axs[1].plot(v_sat_arr[:, 0], t_rels, color="black", linewidth=2, zorder=2)
axs[1].plot(v_sat_arr[:, 1], t_rels, color="black", linewidth=2, zorder=2)

plots_each_max = 20
tmp_list = np.linspace(0, len(v_max[0]) - 1, len(v_max[0]))
mark_indices = np.where(np.mod(tmp_list, plots_each_max) == 0)
axs[0].plot(p_max[1], p_max[0], "-", color="tab:orange", zorder=3)
axs[1].plot(v_max[1], v_max[0], "-", color="tab:orange", zorder=3)
axs[0].scatter(

    p_max[1, mark_indices], p_max[0, mark_indices],
    s=20, edgecolors="tab:orange",
    facecolors='white', zorder=100

)
axs[1].scatter(

    v_max[1, mark_indices], v_max[0, mark_indices],
    s=20, edgecolors="tab:orange",
    facecolors='white', zorder=100

)

for ax in axs:

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("$T_{rel}$ [-]")

axs[0].set_xlabel("$p_{rel}$ [-]")
axs[1].set_xlabel("$(v - b) / (v_{crit} - b)$ [-]")
axs[1].set_xlim((np.min(v_rels), np.max(v_rels)))

plt.tight_layout()
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
filename = "r_dag_behaviour_contour.png"
filepath = os.path.join(CALCULATION_DIR, "3 - System Calculation", "output", filename)
fig.savefig(filepath)
