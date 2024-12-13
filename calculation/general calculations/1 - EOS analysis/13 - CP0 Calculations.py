# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.constant import CALCULATION_DIR
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
from itertools import cycle
from tqdm import tqdm
import numpy as np
import os


# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
n_points = 200
t_mods = np.logspace(-0.5, 0.5, n_points)

fluids = ["Water", "Carbon Dioxide", "Methane"]
eos_colors = ["tab:blue", "tab:green", "tab:orange"]


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
k = 0
lines = list()
fig, axs = plt.subplots(1, 1, figsize=(6, 5), dpi=300)
pbar = tqdm(desc="Calculating Points", total=len(t_mods) * len(fluids))

CP0_rels_ovrl = list()
t_rel_ovrl = list()
CP0_crits = list()

for fluid in fluids:

    tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")
    styles = cycle(["--", "-.", "-", ":"])

    t_crit = tp.RPHandler.TC
    p_crit = tp.RPHandler.PC

    tp.set_variable("T", t_crit)
    tp.set_variable("P", p_crit)
    CP0_crits.append(tp.get_variable("CP0"))

    CP0_rels = np.empty(t_mods.shape)
    CP0_rels[:] = np.nan

    for i in range(len(t_mods)):

        t_curr = t_crit * t_mods[i]
        tp.set_variable("P", p_crit)
        tp.set_variable("T", t_curr)
        CP0 = tp.get_variable("CP0")

        if CP0 > 0:

            CP0_rels[i] = tp.get_variable("CP0") / CP0_crits[-1]

            CP0_rels_ovrl.append(CP0_rels[i])
            t_rel_ovrl.append(t_mods[i])

        pbar.update(1)

    style = next(styles)

    line_eos = axs.plot(t_mods, CP0_rels, linestyle="-", color=eos_colors[k], label=fluid)
    lines.append(line_eos[0])
    k += 1

coef = np.polyfit(t_rel_ovrl, CP0_rels_ovrl, 1)
reg = np.poly1d(coef)

line_eos = axs.plot(t_mods, reg(t_mods), linestyle="--", color="tab:gray", label="Linear", alpha=0.75)
lines.append(line_eos[0])

labs = [l.get_label() for l in lines]

# axs.set_yscale("log")
# axs.set_xscale("log")
axs.grid(visible=True, which="major")
axs.grid(visible=True, which="minor", linewidth=0.75, alpha=0.25)
axs.set_ylabel("$CP0_{rel}$ [-]")
axs.set_xlabel("$t_{rel}$ [-]")
axs.legend(lines, labs)

pbar.close()
plt.tight_layout()
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
filepath = os.path.join(CALCULATION_DIR, "1 - EOS analysis", "output", "13 - CP0 Calculations.png")
fig.savefig(filepath)

