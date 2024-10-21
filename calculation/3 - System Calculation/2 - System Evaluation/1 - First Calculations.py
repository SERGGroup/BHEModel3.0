# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from main_classes import RKEOS, evaluate_system, evaluate_surface


# %%-------------------------------------   CALCULATIONS                        -------------------------------------> #
fluid = RKEOS()

depth = 345.50
grad_rocks = 0.011
v_in = 0.01 * (fluid.v_crit - fluid.b) + fluid.b
in_state = fluid.get_state(t=0.9*fluid.t_crit, v=v_in)

t_in = in_state.t
t_rocks = t_in + depth * grad_rocks
states = evaluate_system(fluid, in_state, depth, t_rocks)
surface_states = evaluate_surface(fluid, states, evaluate_with_flash=False)

w_dot = (states[3].h - states[0].h)
ex_dot = w_dot - states[0].t * (states[3].s - states[0].s)
cf = 1 - t_in / t_rocks

w_dot_nd = w_dot / (states[0].cp * states[0].t)
ex_dot = ex_dot / (states[0].cp * states[0].t)
eta_ex = ex_dot / (w_dot * cf)

w_dex_ratios = [

    (surface_states[0].h - states[0].h) / w_dot,
    (states[3].h - surface_states[1].h) / w_dot

]
