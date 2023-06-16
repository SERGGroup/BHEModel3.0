# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.constant import CALCULATION_DIR
from REFPROPConnector import ThermodynamicPoint
from main_classes import VdWEOS, RKEOS, SRKEOS
import matplotlib.pyplot as plt
from itertools import cycle
from tqdm import tqdm
import numpy as np
import os


# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
n_points = 200
t_mods = [0.5, 0.7, 0.9, 1.5, 2]
p_rels = np.logspace(-4, 1, n_points)

fluids = ["Water", "Carbon Dioxide", "Methane"]
m_mols = [0.01801528, 0.04401, 0.01604]
acntrs = [0.344, 0.239, 0.011]
cps = [1864.84159, 845.846, 2230.12]
eos_class = RKEOS


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
k = 0
x_labels = ["R [J/(kg K)]", "cp [J/(kg K)]"]
fig, axs = plt.subplots(len(x_labels), len(fluids), figsize=(6 * len(fluids), 5 * len(x_labels)), dpi=300)
pbar = tqdm(desc="Calculating Points", total=n_points * len(t_mods) * len(fluids))

for fluid in fluids:

    tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
    styles = cycle(["--", "-.", "-", ":"])

    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    tp.set_variable("T", t_crit)
    tp.set_variable("P", p_crit)

    eos_fluid = eos_class(

        t_crit=t_crit, p_crit=p_crit,
        cp_ideal=tp.get_variable("CP0"),
        m_molar=m_mols[k], acntr=acntrs[k]

    )
    crit_state = eos_fluid.get_state(t=t_crit, p=p_crit)

    for t_mod in t_mods:

        r_res = np.zeros((2, n_points))
        cp_res = np.zeros((2, n_points))
        dpdts = np.zeros((2, n_points))
        dpdvs = np.zeros((2, n_points))

        t_curr = eos_fluid.t_crit * t_mod

        for i in range(len(p_rels)):

            p_curr = p_rels[i] * eos_fluid.p_crit

            tp.set_variable("P", p_curr)
            tp.set_variable("T", t_curr)
            rho_res = tp.get_variable("rho")

            if rho_res < 0:

                r_dag = np.nan
                cp = np.nan
                dpdt = np.nan
                dpdv = np.nan

            else:

                dpdt = tp.get_derivative("P", "T", "rho")
                dpdrho = tp.get_derivative("P", "rho", "T")
                dpdv = - rho_res ** 2 * dpdrho

                r_dag = (- t_curr * dpdt ** 2 / dpdv)
                cp = (tp.get_variable("CP") - tp.get_variable("CP0"))

            r_res[0, i] = r_dag
            cp_res[0, i] = cp
            dpdts[0, i] = dpdt
            dpdvs[0, i] = dpdv

            curr_state = eos_fluid.get_state(p=p_curr, t=t_curr)
            r_res[1, i] = curr_state.r
            cp_res[1, i] = curr_state.cp - curr_state.fluid_solver.cp_ideal
            dpdts[1, i] = curr_state.dpdt
            dpdvs[1, i] = curr_state.dpdv

            pbar.update(1)

        style = next(styles)

        line_rp = axs[0, k].plot(r_res[0,:], p_rels, label="REFPROP", linestyle=style, color="black")
        line_eos = axs[0, k].plot(r_res[1,:], p_rels, label="RK EoS", linestyle=style, color="tab:orange")

        axs[1, k].plot(cp_res[0, :], p_rels, label="REFPROP", linestyle=style, color="black")
        axs[1, k].plot(cp_res[1, :], p_rels, label="RK EoS", linestyle=style, color="tab:orange")

    lines = [line_rp[0], line_eos[0]]
    labs = [l.get_label() for l in lines]

    for j in range(len(x_labels)):

        axs[j, k].set_yscale("log")
        axs[j, k].grid(visible=True, which="major")
        axs[j, k].grid(visible=True, which="minor", linewidth=0.75, alpha=0.25)
        axs[j, k].set_xlabel(x_labels[j])
        axs[j, k].set_ylabel("$p_{rel}$ [-]")
        axs[j, k].set_title(fluid)
        axs[j, k].legend(lines, labs)

    #axs[0, k].set_xlim([0, 1e7])
    #axs[1, k].set_xlim([0, 1e7])

    k += 1


pbar.close()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
plt.tight_layout()
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
filepath = os.path.join(CALCULATION_DIR, "1 - EOS analysis", "output", "8 - EoS Comparison - Other.png")
fig.savefig(filepath)
