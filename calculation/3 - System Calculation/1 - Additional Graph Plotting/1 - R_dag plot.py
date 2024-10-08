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
n_points = 20000
t_mods = np.concatenate((np.logspace(0, 1, 15), np.logspace(-1, 0, 10)[:-1]), axis=None)
v_rels = np.logspace(-3, 4, n_points)

use_base_cp = True
plot_v_rel = True


# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
fluids = ["Methane"]
m_mols = [0.01604]
acntrs = [0.011]

if use_base_cp:

    cps = list()
    # base_coeff = [0.81558724, 0.17805891, 0.01111623]
    base_coeff = [0.81558724]

    for fluid in fluids:

        tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
        tp.set_variable("T", tp.RPHandler.TC)
        tp.set_variable("P", tp.RPHandler.PC)
        cp0_crit = tp.get_variable("CP0")
        t_crit = tp.RPHandler.TC

        coeff = convert_cp_coeff(base_coeff, cp0_crit, t_crit)
        cps.append(coeff)

else:

    cps = [

        [1.20012469e+03, 3.24812968e+00, 7.48129676e-04, - 7.04488778e-07]

    ]

eos_classes = [PREOS]
eos_names = ["PR eos"]
alphas = [1]


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
k = 0
y_labels = ["$R^{\\dagger}\\ /\\ cp$ [-]"]
fig, axs = plt.subplots(len(y_labels), 2 * len(fluids), figsize=(6 * 2 * len(fluids), 5 * len(y_labels)), dpi=300)
pbar = tqdm(desc="Calculating Points", total=n_points * len(t_mods) * len(fluids))

if len(np.shape(axs)) == 1:

    axs = [axs]

elif len(np.shape(axs)) == 0:

    axs = [[axs]]

for fluid in fluids:

    tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
    styles = cycle(["--", "-.", "-", ":"])

    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    eos_fluids = list()
    crit_states = list()

    for eos_class in eos_classes:

        fluid_curr = eos_class(

            t_crit=t_crit, p_crit=p_crit,
            cp_ideal=cps[k], m_molar=m_mols[k],
            acntr=acntrs[k]

        )
        crit_states.append(fluid_curr.get_state(t=t_crit, p=p_crit))
        eos_fluids.append(fluid_curr)

    lines = list()
    max_r = list()
    max_xi = list()
    max_v = list()
    max_xi_ovr = list()
    max_p = list()

    for t_mod in t_mods:

        t_curr = t_crit * t_mod
        res_rels = np.zeros((len(eos_classes), n_points))
        xi_rels = np.zeros((len(eos_classes), n_points))
        p_rels = np.zeros((len(eos_classes), n_points))

        for i in range(len(v_rels)):

            j = 0
            v_curr = v_rels[i] * (eos_fluids[0].v_crit - eos_fluids[0].b) + eos_fluids[0].b
            for curr_fluid in eos_fluids:

                curr_state = curr_fluid.get_state(v=v_curr, t=t_curr)
                xi_rels[j, i] = curr_state.xi
                p_rels[j, i] = curr_state.p / crit_states[j].p

                if curr_state.bifase:
                    res_rels[j, i] = np.nan

                else:
                    res_rels[j, i] = curr_state.r / curr_state.cp

                j += 1

            pbar.update(1)

        style = next(styles)

        if t_mod < 1:

            eos_colors = ["tab:green"]

        else:

            max_value = np.max(res_rels[0, :])
            max_i = np.where(res_rels[0, :] == max_value)[0][0]

            max_r.append(max_value)

            max_v.append(v_rels[max_i])
            max_p.append(p_rels[0, max_i])
            max_xi.append(xi_rels[0, max_i])
            max_xi_ovr.append(np.max(xi_rels))

            if t_mod == 1:

                eos_colors = ["black"]
                style = "-"

            else:

                eos_colors = ["tab:blue"]

        lines = list()
        for j in range(len(eos_fluids)):

            nan_mask = np.isfinite(res_rels[j, :])
            lines.append(axs[0][2 * k].plot(

                p_rels[j, nan_mask], res_rels[j, nan_mask],
                linestyle=style, color=eos_colors[j],
                label=eos_names[j], alpha=alphas[j]

            )[0])

            if plot_v_rel:
                plt_x = v_rels / max_v[0]

            else:
                plt_x = xi_rels[j, :] / max_xi[0]

            axs[0][2 * k + 1].plot(

                plt_x[nan_mask], res_rels[j, nan_mask],
                linestyle=style, color=eos_colors[j],
                label=eos_names[j], alpha=alphas[j]

            )

    labs = [l.get_label() for l in lines]

    for j in range(len(y_labels)):

        t_mods = np.array(t_mods)
        max_r = np.array(max_r)
        max_p = np.array(max_p)
        max_xi = np.array(max_xi) / max_xi[0]
        max_v = np.array(max_v) / max_v[0]
        accpt_indxs = np.where((np.array(max_xi_ovr) > max_xi))[0]

        max_xi[0] = 1
        axs[j][2 * k].set_xlim((10**-2, 10**5))
        ax_ins = axs[j][2 * k].inset_axes(

            [0.13, 0.55, 0.37, 0.42],

        )

        ax_ins.plot(t_mods[accpt_indxs], max_r[accpt_indxs], "o-", color="tab:orange", mfc='white')
        axs[j][2 * k].plot(max_p[accpt_indxs], max_r[accpt_indxs], "o-", color="tab:orange", mfc='white', alpha=0.7)

        if plot_v_rel:
            axs[j][2 * k + 1].invert_xaxis()
            axs[j][2 * k + 1].set_xlabel("$(v - b) / (v_{crit} - b)$ [-]")
            axs[j][2 * k + 1].plot(

                max_v[accpt_indxs], max_r[accpt_indxs], "o-",
                color="tab:orange", mfc='white', alpha=0.7

            )

        else:
            axs[j][2 * k + 1].set_xlabel("$\\xi_{rel}$ [-]")
            axs[j][2 * k + 1].plot(

                max_xi[accpt_indxs], max_r[accpt_indxs],
                "o-", color="tab:orange", mfc='white', alpha=0.7

            )

        ax_ins.set_xlabel("$T_{rel}$ [-]")
        axs[j][2 * k].set_xlabel("$p_{rel}$ [-]")

        for n in range(2):
            axs[j][2 * k + n].invert_xaxis()
            axs[j][2 * k + n].set_ylabel(y_labels[j])
            axs[j][2 * k + n].set_ylim((0, 1))
            axs[j][2 * k + n].set_xscale("log")
            axs[j][2 * k + n].grid(visible=True, which="major")
            axs[j][2 * k + n].grid(visible=True, which="minor", linewidth=0.75, alpha=0.25)

        ax_ins.grid(visible=True, which="major")
        ax_ins.grid(visible=True, which="minor", linewidth=0.75, alpha=0.25)
        ax_ins.set_ylabel("R / cp$_{max}$ [-]")

    k += 1

pbar.close()
plt.tight_layout()
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
if use_base_cp:
    filename = "r_dag_behaviour_fixed_cp.png"
else:
    filename = "r_dag_behaviour.png"

filepath = os.path.join(CALCULATION_DIR, "3 - System Calculation", "output", filename)
fig.savefig(filepath)
