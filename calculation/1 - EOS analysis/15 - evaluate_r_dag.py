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
n_points = 500
t_mods = [0.9, 1, 1.1, 1.2, 1.5,  2, 3, 10]
p_rels = np.logspace(-3, 2, n_points)

use_base_cp = False


# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
fluids = ["Carbon Dioxide", "Methane"]  # ["Water", ["Carbon Dioxide", "Methane"]
m_mols = [0.04401, 0.01604]     # [0.01801528, [0.04401, 0.01604]
acntrs = [0.239, 0.011]     # [0.344, [0.239, 0.011]

if use_base_cp:

    cps = list()
    base_coeff = [0.81558724, 0.17805891, 0.01111623]
    #base_coeff = [0.79655508, 0.21270104]

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

        [1.78959195e+03, 1.06576195e-01, 5.88389412e-04, -1.99830366e-07],
        [4.49897751e+02, 1.66780277e+00, -1.27243808e-03, 3.90820268e-07],
        [1.20012469e+03, 3.24812968e+00, 7.48129676e-04, - 7.04488778e-07]

    ]

eos_classes = [VdWEOS]  # [VdWEOS, RKEOS, SRKEOS, PREOS]
eos_names = ["PR eos"]  # ["VdW eos", "RK eos", "SRK eos", "PR eos"]
eos_colors = ["tab:green"]  # ["tab:blue", "tab:orange", "tab:green", "tab:gray"]
alphas = [1]    # [0.35, 0.5, 1, 1]


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
k = 0
y_labels = ["R / cp [-]"]
fig, axs = plt.subplots(len(y_labels), len(fluids), figsize=(6 * len(fluids), 5 * len(y_labels)), dpi=300)
pbar = tqdm(desc="Calculating Points", total=n_points * len(t_mods) * len(fluids))

if len(np.shape(axs)) == 1:

    axs = [axs]

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
    for t_mod in t_mods:

        t_curr = t_crit * t_mod
        res_rels = np.zeros((len(eos_classes) + 1, n_points))

        for i in range(len(p_rels)):

            p_curr = p_rels[i] * eos_fluids[0].p_crit

            tp.set_variable("P", p_curr)
            tp.set_variable("T", t_curr)
            rho_res = tp.get_variable("rho")

            dpdrho = tp.get_derivative("P", "rho", "T")
            dpdt = tp.get_derivative("P", "T", "rho")
            dpdv = - rho_res ** 2 * dpdrho
            gamma = tp.get_variable("gamma")

            # rp_res = (- t_curr * (dpdt ** 2) / dpdv) / tp.get_variable("cp")
            # rp_res = (tp.get_variable("cp") - tp.get_variable("cv")) / tp.get_variable("cp")
            rp_res = -dpdv

            if rho_res < 0:

                rp_res = np.nan

            res_rels[0, i] = rp_res

            j = 0
            for curr_fluid in eos_fluids:

                curr_state = curr_fluid.get_state(p=p_curr, t=t_curr)
                # res_rels[j + 1, i] = curr_state.r / curr_state.cp
                res_rels[j + 1, i] = -curr_state.dpdv
                j += 1

            pbar.update(1)

        style = next(styles)

        label_curr = "$t_{{rel}}$ = {}[-]".format(t_mod)
        for j in range(len(eos_fluids)):

            lines.append(axs[0][k].plot(

                p_rels, res_rels[1 + j, :],
                linestyle=style, color=eos_colors[j],
                label=label_curr, alpha=alphas[j]

            )[0])

        # lines.append(axs[0][k].plot(p_rels, res_rels[0, :], label=label_curr, linestyle=style, color="black")[0])
        axs[0][k].plot(p_rels, res_rels[0, :], label=label_curr, linestyle=style, color="black")

    labs = [l.get_label() for l in lines]

    for j in range(len(y_labels)):

        axs[j][k].set_title(fluid)
        axs[j][k].set_xscale("log")
        axs[j][k].invert_xaxis()
        axs[j][k].legend(lines, labs)
        axs[j][k].set_ylabel(y_labels[j])
        axs[j][k].set_xlabel("$p_{rel}$ [-]")
        # axs[j][k].set_ylim((0, 1))
        axs[j][k].set_yscale("log")
        axs[j][k].grid(visible=True, which="major")
        axs[j][k].grid(visible=True, which="minor", linewidth=0.75, alpha=0.25)

    k += 1

pbar.close()
plt.tight_layout()
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
if use_base_cp:

    cp_text = "simple cp"

else:

    cp_text = "complex cp"

filepath = os.path.join(CALCULATION_DIR, "1 - EOS analysis", "output", "15 - R_dag_evaluation {}.png".format(cp_text))
fig.savefig(filepath)
