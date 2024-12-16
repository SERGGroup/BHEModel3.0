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
t_rels = np.logspace(np.log10(0.3), np.log10(4), n_points)
t_sats = np.logspace(np.log10(min(t_rels)), 0, n_points)
p_rels = [

    10 ** -6, 10 ** -5, 10 ** -4,
    10 ** -3, 10 ** -2, 10 ** -1,
    10 ** 0,  10 ** 1,  10 ** 2

]


# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
fluid = "Carbon Dioxide"
m_mol = 0.04401
acntr = 0.239

reduced_cp_coeff = False
cps = [4.49897751e+02, 1.66780277e+00, -1.27243808e-03, 3.90820268e-07]

eos_class = PREOS
tp = ThermodynamicPoint([fluid], [1], unit_system="MASS BASE SI")

t_crit = tp.RPHandler.TC
p_crit = tp.RPHandler.PC
t_ref = t_crit*2
p_ref = p_crit*0.01

fluid_curr = eos_class(

    t_crit=t_crit, p_crit=p_crit,
    cp_ideal=cps, m_molar=m_mol,
    acntr=acntr

)

crit_state = fluid_curr.get_state(t=t_crit, p=p_crit)
ref_state = fluid_curr.get_state(t=t_ref, p=p_ref)


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
alpha = 1
eos_name = "PR eos"
eos_color = "tab:green"
styles = cycle(["--", "-.", "-", ":"])
x_label = "$s [kJ/(kg K)]$"
y_label = "$T_{rel}$ [-]"

fig, ax = plt.subplots(1, 1, figsize=(6, 5), dpi=300)
pbar = tqdm(desc="Calculating Points", total=(1 + n_points) * len(p_rels))
lines = list()

for p_rel in p_rels:

    h_rels = np.zeros((1, n_points))
    s_rels = np.zeros((1, n_points))
    p_curr = p_rel * fluid_curr.p_crit

    for i in range(len(t_rels)):

        t_curr = t_rels[i] * fluid_curr.t_crit
        curr_state = fluid_curr.get_state(p=p_curr, t=t_curr)
        h_rels[0, i] = curr_state.h - ref_state.h
        s_rels[0, i] = curr_state.s - ref_state.s

        pbar.update(1)

    style = next(styles)

    h_rels = h_rels / 1e3
    s_rels = s_rels / 1e3
    new_line = ax.plot(

        s_rels[0][:], t_rels,
        linestyle=style, color=eos_color,
        alpha=alpha, label="$p_{{rel}}=10^{{{:.0f}}}$".format(np.log10(p_rel))

    )[0]
    lines.append(new_line)

s_sat_liq = np.zeros(len(t_sats))
s_sat_vap = np.zeros(len(t_sats))
for i in range(len(t_sats)):

    t_curr = t_sats[i] * fluid_curr.t_crit
    liq_state, vap_state = fluid_curr.get_sat_state(t=t_curr)
    s_sat_liq[i] = liq_state.s - ref_state.s
    s_sat_vap[i] = vap_state.s - ref_state.s
    pbar.update(1)

s_sat_liq = s_sat_liq / 1e3
s_sat_vap = s_sat_vap / 1e3

ax.plot(s_sat_liq, t_sats, color="black")
ax.plot(s_sat_vap, t_sats, color="black")

pbar.close()

labs = [l.get_label() for l in lines]

ax.set_yscale("log")
ax.grid(visible=True, which="major")
ax.grid(visible=True, which="minor", linewidth=0.75, alpha=0.25)
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
ax.legend(lines, labs, fontsize=8)
plt.tight_layout()
plt.show()


# %%-------------------------------------   SAVE PLOT                           -------------------------------------> #
filepath = os.path.join(CALCULATION_DIR, "1 - EOS analysis", "output", "14 - TS Diagram.png")
fig.savefig(filepath)
