from cubic_eos import CubicEOS, FluidState
from scipy.integrate import RK45
import scipy.constants


def calculate_vertical(fluid_eos: CubicEOS, state_in: FluidState, res_depth: float, downward=True):
    def rk_overall_der(z, y):

        t_curr = y[0]
        v_curr = y[1]

        state = fluid_eos.get_state(t=t_curr, v=v_curr)

        dp = - scipy.constants.g / v_curr
        dt = state.r / state.cp / state.dpdt * dp
        dv = (1 - state.r / state.cp) / state.dpdv * dp

        return [dt, dv]

    if downward:

        res_depth = -abs(res_depth)

    else:

        res_depth = abs(res_depth)

    integrator = RK45(rk_overall_der, 0, [state_in.t, state_in.v], res_depth, rtol=1e-06, atol=1e-07)

    while integrator.status == 'running':

        integrator.step()

    output = integrator.y

    return fluid_eos.get_state(t=output[0], v=output[1])


def evaluate_system(fluid_eos: CubicEOS, in_state: FluidState, depth_res: float, t_res: float):

    res_in = calculate_vertical(fluid_eos, in_state, depth_res, downward=True)
    res_out = fluid_eos.get_state(p=res_in.p, t=t_res)
    sys_out = calculate_vertical(fluid_eos, res_out, depth_res, downward=False)

    return [in_state, res_in, res_out, sys_out]
