# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes.constant import CALCULATION_DIR
from main_classes.cubic_eos import CubicEOS
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import os


# %%-------------------------------------   T_SAT CALCULATIONS                  -------------------------------------> #
fluid = CubicEOS(

    p_crit=1e6,
    t_crit=1e5,
    cp_ideal=850,
    m_molar=0.044

)

p_rels = np.logspace(0, -30, num=1500)
t_rels = np.zeros(p_rels.shape)

for i in range(len(p_rels)):

    t_sat, zl, zv= fluid.t_sat(p_rels[i] * fluid.p_crit)
    t_rels[i] = t_sat / fluid.t_crit

# %%-------------------------------------   PLOT RESULTS                        -------------------------------------> #

reg = np.poly1d(np.polyfit(np.log(p_rels), np.log(t_rels), 10))
def f_reg(x_in):

    x = np.log(x_in)
    y = 0.

    for c in reg.coeffs:

        y = y * x + c

    return np.exp(y)

fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

for ax in axs:

    ax.plot(p_rels, t_rels, label="$T_{sat}$")
    ax.plot(p_rels, f_reg(p_rels), label="Reg. Results")
    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("$p_{rel}$ [-]")
    ax.set_ylabel("$T_{rel}$ [-]")
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))
    ax.yaxis.set_minor_formatter(ticker.StrMethodFormatter("{x:.2f}"))
    ax.grid(visible=True, which="both")
    ax.legend()

cu_lims = {

    "x": (0.5, 1.1),
    "y": (0.9,  1.01),

}

axs[-1].set_xlim(cu_lims["x"])
axs[-1].set_ylim(cu_lims["y"])
axs[-1].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))
axs[-1].xaxis.set_minor_formatter(ticker.StrMethodFormatter("{x:.2f}"))

axs[0].add_patch(

    Rectangle(

        (min(cu_lims["x"]),min(cu_lims["y"])),
        max(cu_lims["x"]) - min(cu_lims["x"]),
        max(cu_lims["y"]) - min(cu_lims["y"]),
        edgecolor='black',
        linestyle='dashed',
        fill=False,
        lw=1.,
        zorder=100

    )

)
plt.tight_layout()
plt.show()

# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
fig.savefig(os.path.join(CALCULATION_DIR, "EOS analysis", "output", "6 - Saturation Polinomial.png"))