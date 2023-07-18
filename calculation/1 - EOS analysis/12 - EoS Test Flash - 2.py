# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import RKEOS, SRKEOS, PREOS
from main_classes import calculate_flash
import matplotlib.pyplot as plt
import numpy as np


# %%-------------------------------------   INIT CALCULATIONS                   -------------------------------------> #
tr_in = 0.7
p_ratio = 0.9
n_points = 100
eos_class = RKEOS


# %%-------------------------------------   PLOT FLASH ITERATION                -------------------------------------> #
fluid = eos_class()
in_state = fluid.get_sat_state(t=tr_in * fluid.t_crit, which="liq")

p_out = p_ratio * in_state.p
out_liq, out_vap = fluid.get_sat_state(p=p_out, which="both")

v_rels = np.linspace(0, 1, n_points + 2)
v_rels = v_rels[1:-1]
entropies = np.full((n_points, 2), np.nan)

for i in range(n_points):

    v_rel = v_rels[i]

    for j in [0, 1]:

        try:

            if j == 0:
                v = v_rel * out_liq.v

            else:
                v = 1 / (v_rel / out_vap.v)

            out_state_curr = fluid.get_state(p=p_out, v=v)

        except Exception as exception:

            print(exception)

        else:

            entropies[i, j] = out_state_curr.h


fig, ax = plt.subplots(1, 1, figsize=(10, 8))

ax.plot(v_rels, entropies[:, 0], label='liquid')
ax.plot(v_rels, entropies[:, 1], label='vapour')
ax.legend()

plt.tight_layout()
plt.show()
