# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import VdWEOS, RKEOS, SRKEOS, PREOS
from main_classes.constant import CALCULATION_DIR
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
from itertools import cycle
from tqdm import tqdm
import numpy as np
import os


# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
n_points = 200
t_mods = [0.5, 0.7, 0.9, 1, 2]
p_rels = np.logspace(-4, 1, n_points)

fluids = ["Water", "Carbon Dioxide", "Methane"]
m_mols = [0.01801528, 0.04401, 0.01604]
acntrs = [0.344, 0.239, 0.011]
cps = [1864.84159, 845.846, 2230.12]

eos_classes = [VdWEOS, RKEOS, SRKEOS, PREOS]
eos_names = ["VdW EoS", "RK EoS", "SRK EoS", "PR EoS"]
eos_colors = ["tab:blue", "tab:orange", "tab:green", "tab:gray"]
alphas = [0.35, 1, 1, 1]


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
k = 0
x_labels = ["h [kJ/kg]", "s [kJ/(kg K)]"]
fig, axs = plt.subplots(len(x_labels), len(fluids), figsize=(6 * len(fluids), 5 * len(x_labels)), dpi=300)
pbar = tqdm(desc="Calculating Points", total=n_points * len(t_mods) * len(fluids))

for fluid in fluids:

    tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
    styles = cycle(["--", "-.", "-", ":"])

    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    tp.set_variable("T", t_crit)
    tp.set_variable("P", p_crit)
    h_crit_rp = tp.get_variable("h")
    s_crit_rp = tp.get_variable("s")

    fluids = list()
    crit_states = list()

    for eos_class in eos_classes:

        fluid_curr = eos_class(t_crit=t_crit, p_crit=p_crit, cp_ideal=cps[k], m_molar=m_mols[k], acntr=acntrs[k])
        crit_states.append(fluid_curr.get_state(t=t_crit, p=p_crit))
        fluids.append(fluid_curr)

    for t_mod in t_mods:

        h_rels = np.zeros((5, n_points))
        s_rels = np.zeros((5, n_points))

        t_curr = fluids[0].t_crit * t_mod

        for i in range(len(p_rels)):

            p_curr = p_rels[i] * fluids[0].p_crit

            tp.set_variable("P", p_curr)
            tp.set_variable("T", t_curr)
            h_res = tp.get_variable("h")
            s_res = tp.get_variable("s")
            rho_res = tp.get_variable("rho")

            if rho_res < 0:

                h_res = np.nan
                s_res = np.nan

            h_rels[0, i] = h_res - h_crit_rp
            s_rels[0, i] = s_res - s_crit_rp

            j = 0
            for curr_fluid in fluids:

                curr_state = curr_fluid.get_state(p=p_curr, t=t_curr)
                h_rels[j + 1, i] = curr_state.h - crit_states[j].h
                s_rels[j + 1, i] = curr_state.s - crit_states[j].s
                j += 1

            pbar.update(1)

        style = next(styles)

        h_rels = h_rels / 1e3
        s_rels = s_rels / 1e3

        lines = list()
        for j in range(len(fluids)):

            dh = h_rels[1 + j, 0] - h_rels[0, 0]
            ds = s_rels[1 + j, 0] - s_rels[0, 0]

            axs[1, k].plot(s_rels[1 + j, :] - ds, p_rels, linestyle=style, color=eos_colors[j], alpha=alphas[j])
            line_eos = axs[0, k].plot(

                h_rels[1 + j, :] - dh, p_rels,
                linestyle=style, color=eos_colors[j],
                label=eos_names[j], alpha=alphas[j]

            )

            lines.append(line_eos[0])

        line_rp = axs[0, k].plot(h_rels[0, :], p_rels, label="REFPROP", linestyle=style, color="black")
        axs[1, k].plot(s_rels[0, :], p_rels, linestyle=style, color="black")
        lines.append(line_rp[0])

    labs = [l.get_label() for l in lines]

    for j in range(len(x_labels)):

        axs[j, k].set_yscale("log")
        axs[j, k].grid(visible=True, which="major")
        axs[j, k].grid(visible=True, which="minor", linewidth=0.75, alpha=0.25)
        axs[j, k].set_xlabel(x_labels[j])
        axs[j, k].set_ylabel("$p_{rel}$ [-]")
        axs[j, k].set_title(fluid)
        axs[j, k].legend(lines, labs)

    k += 1


pbar.close()
plt.tight_layout()
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
filepath = os.path.join(CALCULATION_DIR, "1 - EOS analysis", "output", "8 - EoS Comparison - Enthalpy and Entropy.png")
fig.savefig(filepath)
