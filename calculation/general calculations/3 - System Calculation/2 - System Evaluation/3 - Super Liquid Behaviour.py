# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import RKEOS, evaluate_surface, evaluate_system
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
import scipy.constants
from tqdm import tqdm
import numpy as np


# %%-------------------------------------   INIT ARRAYS                         -------------------------------------> #
n_grad = 150
n_depth = 3
n_t_rel = 5

grad_nd_list = np.logspace(0, 2, n_grad)
dz_nd_list = np.logspace(-2, -1, n_depth)
t_rel_list = np.logspace(np.log10(0.5), np.log10(20), n_t_rel)

grad_nd, dz_nd = np.meshgrid(grad_nd_list, dz_nd_list, indexing="ij")

carnot_factor = 1 - 1 / (1 + grad_nd * dz_nd)

spc_work_gas = (grad_nd - 1) * dz_nd
spc_ex_gas = spc_work_gas - np.log((1 + grad_nd * dz_nd)/(1 + dz_nd))
ex_eta_gas = spc_ex_gas / (spc_work_gas * carnot_factor)

spc_work_liq = grad_nd * dz_nd
spc_ex_liq = spc_work_liq - np.log(1 + spc_work_liq)
ex_eta_liq = spc_ex_liq / (spc_work_liq * carnot_factor)

if np.isnan(ex_eta_gas[0, 0]):

    ex_eta_gas[0, :] = np.ones(ex_eta_gas[0, :].shape)

base_shape = np.array(grad_nd.shape)
res_shape = np.append([len(t_rel_list)], base_shape)

r_cp_in_arr = np.empty(len(t_rel_list))
r_cp_in_arr[:] = np.nan

grad_rocks = np.empty(res_shape)
depth = np.empty(res_shape)
t_rocks = np.empty(res_shape)
w_dot_nds = np.empty(res_shape)
ex_dot_nds = np.empty(res_shape)
eta_exs = np.empty(res_shape)
w_dex_mins = np.empty(res_shape)
w_dex_maxs = np.empty(res_shape)

grad_rocks[:] = np.nan
depth[:] = np.nan
t_rocks[:] = np.nan
w_dot_nds[:] = np.nan
ex_dot_nds[:] = np.nan
eta_exs[:] = np.nan
w_dex_mins[:] = np.nan
w_dex_maxs[:] = np.nan


# %%-------------------------------------   CALCULATE                           -------------------------------------> #
p_rel = 10**2
tp = ThermodynamicPoint(["Methane"], [1], unit_system="MASS BASE SI")
fluid = RKEOS(

    p_crit=tp.RPHandler.PC, t_crit=tp.RPHandler.TC,
    cp_ideal=845.85, m_molar=0.01604,
    acntr=0.011

)

bhe_list = list()
pbar = tqdm(desc="Calculating Points", total=len(t_rel_list) * len(grad_nd_list) * len(dz_nd_list))
for i in range(len(t_rel_list)):

    t_in = t_rel_list[i] * tp.RPHandler.TC
    p_in = p_rel * tp.RPHandler.PC
    in_state = fluid.get_state(t=t_in, p=p_in)

    r_cp_in_arr[i] = 1 - in_state.r / in_state.cp
    grad_rocks[i, :, :] = grad_nd * scipy.constants.g / in_state.cp
    depth[i, :, :] = dz_nd * in_state.cp * t_in / scipy.constants.g
    t_rocks[i, :, :] = t_in + depth[i, :, :] * grad_rocks[i, :, :]

    bhe_sub_list = list()
    for j in range(len(grad_nd_list)):

        bhe_subsub_list = list()
        for k in range(len(dz_nd_list)):

            try:

                states = evaluate_system(fluid, in_state, depth[i, j, k], t_rocks[i, j, k])
                surface_states = evaluate_surface(fluid, states, evaluate_with_flash=False)
                bhe_subsub_list.append(states)

            except:

                pass

            else:

                w_dot = states[3].h - states[0].h
                ex_dot = w_dot - states[0].t * (states[3].s - states[0].s)
                ex_dot_nd = ex_dot / (states[0].cp * states[0].t)

                if not (w_dot < 0 or ex_dot < 0):

                    w_dot_nds[i, j, k] = w_dot / (states[0].cp * states[0].t)
                    ex_dot_nds[i, j, k] = ex_dot / (states[0].cp * states[0].t)
                    eta_exs[i, j, k] = ex_dot / (w_dot * carnot_factor[j, k])

                    w_dex_mins[i, j, k] = (surface_states[1].h - states[0].h) / w_dot
                    w_dex_maxs[i, j, k] = (states[3].h - surface_states[1].h) / w_dot

            pbar.update(1)

        bhe_sub_list.append(bhe_subsub_list)

    bhe_list.append(bhe_sub_list)

pbar.close()


# %%-------------------------------------   PLOT INITIAL COMPARISON             -------------------------------------> #
fig, base_axs = plt.subplots(2, 2, dpi=300)
fig.set_size_inches(10, 5)
gs = base_axs[0, 1].get_gridspec()
cmap = plt.get_cmap('viridis')
norm = plt.Normalize(0, 1)

for ax in base_axs[:, 1]:

    ax.remove()

axbig = fig.add_subplot(gs[:, 1])
axs = [base_axs[0, 0], base_axs[1, 0], axbig]
x_values = grad_nd - 1

t_rel_label = "$T_{{rel}} = 10^{{ {:0.1f} }}$"
cmp_labels = ["Liquid", t_rel_label, "Ideal Gas"]
cmp_y_values = [

    [spc_work_liq, spc_ex_liq, ex_eta_liq],
    [w_dot_nds, ex_dot_nds, eta_exs],
    [spc_work_gas, spc_ex_gas, ex_eta_gas],

]

k = 2
lines = [list(), list(), list()]

for m in range(len(cmp_y_values[0])):

    for n in range(len(cmp_labels)):

        if cmp_labels[n] == t_rel_label:

            for i in range(len(t_rel_list)):

                lines[m].append(

                    axs[m].plot(

                        x_values[:, k], cmp_y_values[n][m][i, :, k], "-",
                        label=t_rel_label.format(np.log10(t_rel_list[i])),
                        color=cmap(norm((i + 1)/(len(t_rel_list) + 1)))

                    )[0]

                )

        else:

            lines[m].append(

                axs[m].plot(

                    x_values[:, k], cmp_y_values[n][m][:, k], "-",
                    label=cmp_labels[n], color=cmap(norm(n)),
                    linewidth=2

                )[0]

            )


y_names = ["${\\dot{w}}^{\\#}$ [-]", "${\\dot{e}_x}^{\\#}$ [-]", "${\\eta}_{ex}$ [-]"]
axs[0].get_xaxis().set_visible(False)

for k in range(len(axs)):

    axs[k].set_xscale("log")
    axs[k].set_ylabel(y_names[k])
    axs[k].set_xlabel("${\\nabla T_{rocks}}^{\\#} - 1$ [-]")

    if not y_names[k] == y_names[-1]:
        axs[k].set_yscale("log")

    if not y_names[k] == y_names[1]:
        axs[k].legend(handles=lines[k], fontsize="8")

axs[-1].set_ylim((

    0.3,
    1

))

plt.tight_layout(pad=2)
plt.subplots_adjust(hspace=0)
plt.show()


# %%-------------------------------------   EVALUATE COMPRESSIBILITY EFFECT     -------------------------------------> #
t_in = 5 * tp.RPHandler.TC
p_in = 0.1 * tp.RPHandler.PC
in_state = fluid.get_state(t=t_in, p=p_in)
ig_comp = 1 - in_state.r / in_state.cp

comp_rel = (ig_comp - r_cp_in_arr) / ig_comp

plt.plot(t_rel_list, eta_exs[:, 0, 0])
plt.xscale("log")
plt.show()
