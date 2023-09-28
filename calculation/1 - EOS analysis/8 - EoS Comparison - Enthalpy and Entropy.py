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
n_points = 200
t_mods = [0.5, 0.7, 0.9, 1, 2]
p_rels = np.logspace(-4, 1, n_points)

use_base_cp = False
simplify_reading = False


# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
# fluids = ["Water", "Carbon Dioxide", "Methane"]
# m_mols = [0.01801528, 0.04401, 0.01604]
# acntrs = [0.344, 0.239, 0.011]
# cps = [
#
#     [1.78959195e+03, 1.06576195e-01, 5.88389412e-04, -1.99830366e-07],
#     [4.49897751e+02, 1.66780277e+00, -1.27243808e-03, 3.90820268e-07],
#     [1.20012469e+03, 3.24812968e+00, 7.48129676e-04, - 7.04488778e-07]
#
# ]

fluids = ["Carbon Dioxide", "Methane"]
m_mols = [0.04401, 0.01604]
acntrs = [0.239, 0.011]
cps = [

    [4.49897751e+02, 1.66780277e+00, -1.27243808e-03, 3.90820268e-07],
    [1.20012469e+03, 3.24812968e+00, 7.48129676e-04, - 7.04488778e-07]

]

if use_base_cp:

    cps = list()
    # base_coeff = [0.81558724, 0.17805891, 0.01111623]
    base_coeff = [0.79655508, 0.21270104]

    for fluid in fluids:

        tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
        tp.set_variable("T", tp.RPHandler.TC)
        tp.set_variable("P", tp.RPHandler.PC)
        cp0_crit = tp.get_variable("CP0")
        t_crit = tp.RPHandler.TC

        coeff = convert_cp_coeff(base_coeff, cp0_crit, t_crit)
        cps.append(coeff)

eos_classes = [VdWEOS, RKEOS, SRKEOS, PREOS]
eos_names = ["VdW eos", "RK eos", "SRK eos", "PR eos"]
eos_colors = ["tab:blue", "tab:orange", "tab:green", "tab:gray"]
alphas = [0.35, 0.5, 1, 1]


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
    t_ref = t_crit*2
    p_ref = p_crit*0.01

    eos_fluids = list()
    crit_states = list()
    ref_states = list()

    for eos_class in eos_classes:

        fluid_curr = eos_class(

            t_crit=t_crit, p_crit=p_crit,
            cp_ideal=cps[k], m_molar=m_mols[k],
            acntr=acntrs[k]

        )
        crit_states.append(fluid_curr.get_state(t=t_crit, p=p_crit))
        ref_states.append(fluid_curr.get_state(t=t_ref, p=p_ref))
        eos_fluids.append(fluid_curr)

    tp.set_variable("T", t_ref)
    tp.set_variable("P", p_ref)
    h_ref_rp = tp.get_variable("h")
    s_ref_rp = tp.get_variable("s")

    for t_mod in t_mods:

        h_rels = np.zeros((5, n_points))
        s_rels = np.zeros((5, n_points))

        t_curr = t_crit * t_mod

        for i in range(len(p_rels)):

            p_curr = p_rels[i] * eos_fluids[0].p_crit

            tp.set_variable("P", p_curr)
            tp.set_variable("T", t_curr)
            h_res = tp.get_variable("h")
            s_res = tp.get_variable("s")
            rho_res = tp.get_variable("rho")

            if rho_res < 0:

                h_res = np.nan
                s_res = np.nan

            h_rels[0, i] = h_res - h_ref_rp
            s_rels[0, i] = s_res - s_ref_rp

            j = 0
            for curr_fluid in eos_fluids:

                curr_state = curr_fluid.get_state(p=p_curr, t=t_curr)
                h_rels[j + 1, i] = curr_state.h - ref_states[j].h
                s_rels[j + 1, i] = curr_state.s - ref_states[j].s
                j += 1

            pbar.update(1)

        style = next(styles)

        h_rels = h_rels / 1e3
        s_rels = s_rels / 1e3

        lines = list()
        for j in range(len(eos_fluids)):

            if simplify_reading:
                dh = h_rels[1 + j, 0] - h_rels[0, 0]
                ds = s_rels[1 + j, 0] - s_rels[0, 0]

            else:
                dh = 0.
                ds = 0.

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
if use_base_cp:

    cp_text = "simple cp"

else:

    cp_text = "complex cp"

filepath = os.path.join(CALCULATION_DIR, "1 - EOS analysis", "output", "8 - HS comparison {}.png".format(cp_text))
fig.savefig(filepath)
