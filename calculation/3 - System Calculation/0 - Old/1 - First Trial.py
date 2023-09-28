# %%-------------------------   IMPORT MODULES                      -------------------------------------------------> #
from main_classes.geothermal_system.simple_geothermal_system import evaluate_system, evaluate_surface
from main_classes.eos.subclasses import SRKEOS, RKEOS


# %%-------------------------   INITIALIZE SYSTEM                   -------------------------------------------------> #
fluid = RKEOS()

depth = 3500.
grad_rel = 400
grad = grad_rel * fluid.t_crit / depth
t_rel_in = 0.85

t_inj = t_rel_in * fluid.t_crit
t_rock = t_inj + grad * depth / 1000

liq_state, vap_state = fluid.get_sat_state(t=t_inj)
input_state = fluid.get_state(t=liq_state.t, v=liq_state.v * 0.99)

sys_states = evaluate_system(fluid, input_state, depth, t_rock)
res = evaluate_surface(fluid, sys_states)
