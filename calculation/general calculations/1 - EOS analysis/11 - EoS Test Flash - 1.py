# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import VdWEOS, RKEOS, SRKEOS, PREOS
from main_classes import calculate_flash


# %%-------------------------------------   CALCULATE LIQUID                    -------------------------------------> #
tr_in = 0.8

error = list()
for eos_class in [RKEOS, SRKEOS, PREOS]:

    fluid = eos_class()
    in_state = fluid.get_sat_state(t=tr_in * fluid.t_crit, which="liq")
    out_state = calculate_flash(fluid, "SP", in_state.s, 2*in_state.p)
    error.append(out_state.s - in_state.s)

print(error)

# %%-------------------------------------   CALCULATE SATURATION                -------------------------------------> #
tr_in = 0.8

error = list()
for eos_class in [RKEOS, SRKEOS, PREOS]:

    fluid = eos_class()
    in_state = fluid.get_sat_state(t=tr_in * fluid.t_crit, which="liq")
    out_state = calculate_flash(fluid, "SP", in_state.s, in_state.p)
    error.append(out_state.s - in_state.s)

print(error)