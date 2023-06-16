# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import VdWEOS, RKEOS, SRKEOS, PREOS
from main_classes import iteration_flash


# %%-------------------------------------   CALCULATE                           -------------------------------------> #
tr_in = 0.8

error = list()
for eos_class in [VdWEOS, RKEOS, SRKEOS, PREOS]:

    fluid = eos_class()
    in_state = fluid.get_sat_state(t=tr_in * fluid.t_crit, which="liq")
    out_state = iteration_flash(fluid, in_state.s, search_h=False, p=0.7*in_state.p)
    error.append(out_state.s - in_state.s)

print(error)