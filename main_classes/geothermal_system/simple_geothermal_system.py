from main_classes.support.simple_integrator import SimpleIntegrator
from main_classes.eos import CubicEOS, FluidState, calculate_flash
from scipy.integrate import RK45
import scipy.constants


def calculate_expansion(fluid_eos: CubicEOS, state_in: FluidState, p_out):

    state = fluid_eos.get_state(t=state_in.t, v=state_in.v)

    def rk_overall_der(z, y):

        v_curr = y[0]

        state.update_state(t=fluid_eos.t(v=v_curr, p=z), v=v_curr)
        dv = (1 - state.r / state.cp) / state.dpdv

        return [dv]

    if p_out < fluid_eos.p_crit:

        liq_state, vap_state = fluid_eos.get_sat_state(p=p_out)
        x_curr = (state_in.s - liq_state.s) / (vap_state.s - liq_state.s)

        if 0 < x_curr < 1:

            v_curr = x_curr * (vap_state.v - liq_state.v) + liq_state.v
            state.update_state(p=p_out, v=v_curr)
            return state

    integrator = RK45(rk_overall_der, state_in.p, [state_in.v], p_out)
    output = integrator.y
    state.update_state(p=p_out, v=output[0])

    return state


def calculate_vertical(fluid_eos: CubicEOS, state_in: FluidState, res_depth: float, downward=True):

    state = fluid_eos.get_state(t=state_in.t, v=state_in.v)

    def rk_overall_der(z, y):

        p_curr = y[0]
        v_curr = y[1]

        state.update_state(t=fluid_eos.t(v=v_curr, p=p_curr), v=v_curr)

        dp = - scipy.constants.g / v_curr
        dv = (1 - state.r / state.cp) / state.dpdv * dp

        return [dp, dv]

    def rk_overall_bifase(z, y):

        p_curr = y[0]
        h_curr = state_in.h - z * scipy.constants.g

        liq_state, vap_state = fluid_eos.get_sat_state(p=p_curr)
        x_curr = (h_curr - liq_state.h) / (vap_state.h - liq_state.h)
        v_curr = x_curr * (vap_state.v - liq_state.v) + liq_state.v

        state.update_state(p=p_curr, v=v_curr)
        dp = - scipy.constants.g / v_curr

        return [dp]

    if downward:

        res_depth = -abs(res_depth)

    else:

        res_depth = abs(res_depth)

    i = 0.
    prev_z = 0.
    while True:

        if state.bifase:

            integrator = RK45(rk_overall_bifase, prev_z, [state.p], res_depth)

        else:

            integrator = RK45(rk_overall_der, prev_z, [state.p, state.v], res_depth)

        prev_bifase_condition = state.bifase
        while integrator.status == 'running':

            if not state.bifase == prev_bifase_condition:

                break

            else:

                prev_z = integrator.t
                integrator.step()

        if not integrator.status == 'running':

            output = integrator.y
            break

        else:

            if i > 25:
                output = integrator.y
                break

            i += 1

    return fluid_eos.get_state(p=output[0], v=output[1])


def evaluate_system(fluid_eos: CubicEOS, in_state: FluidState, depth_res: float, t_res: float):

    res_in = calculate_vertical(fluid_eos, in_state, depth_res, downward=True)
    res_out = fluid_eos.get_state(p=res_in.p, t=t_res)
    sys_out = calculate_vertical(fluid_eos, res_out, depth_res, downward=False)

    return [in_state, res_in, res_out, sys_out]


def evaluate_surface(fluid_eos: CubicEOS, geo_system_points, evaluate_with_flash=False):

    geo_in = geo_system_points[0]
    geo_out = geo_system_points[-1]

    if evaluate_with_flash:

        # Expansion Work Evaluation
        exp_out = calculate_flash(fluid_eos, "PS", geo_in.p, geo_out.s)
        cool_out = calculate_flash(fluid_eos, "PS", geo_out.p, geo_in.s)

    else:

        # Expansion Work Evaluation
        exp_out = calculate_expansion(fluid_eos, geo_out, geo_in.p)
        cool_out = calculate_expansion(fluid_eos, geo_in, geo_out.p)

    return [cool_out, exp_out]
