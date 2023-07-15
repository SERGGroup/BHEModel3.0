# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import VdWEOS, RKEOS, SRKEOS, PREOS
from main_classes import iteration_flash


# %%-------------------------------------   CALCULATE                           -------------------------------------> #
tr_in = 0.5
pr_in = 1.2
error = list()

for eos_class in [VdWEOS, RKEOS, SRKEOS, PREOS]:

    fluid = eos_class()
    in_state = fluid.get_sat_state(t=tr_in * fluid.t_crit, which="liq")
    out_state = iteration_flash(fluid, in_state.s, search_h=False, p=pr_in * in_state.p)
    error.append(out_state.s - in_state.s)

print(error)