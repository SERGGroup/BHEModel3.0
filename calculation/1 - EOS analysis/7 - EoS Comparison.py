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

# fluids = ["Water", "Carbon Dioxide", "Methane"]
# m_mols = [0.01801528, 0.04401, 0.01604]
# acntrs = [0.344, 0.239, 0.011]
# cps = [1864.84159, 845.846, 2230.12]
fluids = ["Carbon Dioxide", "Methane"]
m_mols = [0.04401, 0.01604]
acntrs = [0.239, 0.011]
cps = [845.846, 2230.12]

eos_names = ["VdW eos", "RK eos", "SRK eos", "PR eos"]
eos_colors = ["tab:blue", "tab:orange", "tab:green", "tab:gray"]
alphas = [0.35, 1, 1, 1]


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
k = 0
fig, axs = plt.subplots(1, len(fluids), figsize=(6 * len(fluids), 5), dpi=300)
pbar = tqdm(desc="Calculating Points", total=n_points * len(t_mods) * len(fluids))

for fluid in fluids:

    tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
    styles = cycle(["--", "-.", "-", ":"])
    ax = axs[k]

    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    tp.set_variable("T", t_crit)
    tp.set_variable("P", p_crit)
    v_crit_rp = 1 / tp.get_variable("rho")

    vdw_fluid = VdWEOS(t_crit=t_crit, p_crit=p_crit, cp_ideal=cps[k], m_molar=m_mols[k], acntr=acntrs[k])
    srk_fluid = SRKEOS(t_crit=t_crit, p_crit=p_crit, cp_ideal=cps[k], m_molar=m_mols[k], acntr=acntrs[k])
    rk_fluid = RKEOS(t_crit=t_crit, p_crit=p_crit, cp_ideal=cps[k], m_molar=m_mols[k], acntr=acntrs[k])
    pr_fluid = PREOS(t_crit=t_crit, p_crit=p_crit, cp_ideal=cps[k], m_molar=m_mols[k], acntr=acntrs[k])
    k += 1

    for t_mod in t_mods:

        v_rels = np.zeros((5, n_points))
        t_curr = vdw_fluid.t_crit * t_mod

        for i in range(len(p_rels)):

            p_curr = p_rels[i] * vdw_fluid.p_crit

            tp.set_variable("P", p_curr)
            tp.set_variable("T", t_curr)
            v_real = 1 / tp.get_variable("rho")

            if v_real < 0:

                v_real = np.nan

            v_rels[0, i] = vdw_fluid.v(t=t_curr, p=p_curr) / v_crit_rp
            v_rels[1, i] = rk_fluid.v(t=t_curr, p=p_curr) / v_crit_rp
            v_rels[2, i] = srk_fluid.v(t=t_curr, p=p_curr) / v_crit_rp
            v_rels[3, i] = pr_fluid.v(t=t_curr, p=p_curr) / v_crit_rp
            v_rels[4, i] = v_real / v_crit_rp

            pbar.update(1)

        style = next(styles)
        lines = list()

        for j in range(len(eos_names)):

            line_eos = ax.plot(

                v_rels[j, :], p_rels,
                label=eos_names[j], linestyle=style,
                color=eos_colors[j], alpha=alphas[j]

            )

            lines.append(line_eos[0])

        line_rp = ax.plot(v_rels[4, :], p_rels, label="REFPROP", linestyle=style, color="black")
        lines.append(line_rp[0])

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(visible=True, which="major")
    ax.grid(visible=True, which="minor", linewidth=0.75, alpha=0.25)
    ax.set_xlabel("$v_{rel}$ [-]")
    ax.set_ylabel("$p_{rel}$ [-]")
    ax.set_title(fluid)

    labs = [l.get_label() for l in lines]
    ax.legend(lines, labs)

pbar.close()

plt.tight_layout()
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
fig.savefig(os.path.join(CALCULATION_DIR, "1 - EOS analysis", "output", "7 - eos Comparison.png"))