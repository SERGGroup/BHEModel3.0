from main_classes.cubic_eos import CubicEOS
from scipy.optimize import root_scalar
import numpy as np


def flash(fluid: CubicEOS, input_code, input_1, input_2):

    pass


def direct_flash(fluid: CubicEOS, t=None, p=None, v=None):

    return fluid.get_state(t=t, p=p, v=v)


def base_quality_flash(fluid: CubicEOS, q, t=None, p=None):

    liq_state, vap_state = fluid.get_sat_state(p=p, t=t)
    v = (vap_state.v - liq_state.v) * q + liq_state.v
    return fluid.get_state(p=p, t=t, v=v)


def iteration_flash(fluid: CubicEOS, target_value, search_h=True, t=None, p=None):

    # Check Saturation
    liq_state, vap_state = fluid.get_sat_state(p=p, t=t)

    if search_h:

        liq_value = liq_state.h
        vap_value = vap_state.h

        def sub_err_function(v):

            state = fluid.get_state(p=p, t=t, v=v)
            return target_value - state.h

    else:

        liq_value = liq_state.s
        vap_value = vap_state.s

        def sub_err_function(v):

            state = fluid.get_state(p=p, t=t, v=v)
            return target_value - state.s

    if liq_value < target_value < vap_value:

        q = (target_value - liq_value) / (vap_value - liq_value)
        v = (vap_state.v - liq_state.v) * q + liq_state.v
        return fluid.get_state(p=p, t=t, v=v)

    if target_value < liq_value:

        interval = [fluid.b, liq_state.v]
        err_function = sub_err_function

    else:

        interval = [0, 1 / vap_state.v]
        err_function = lambda rho: sub_err_function(1 / rho)

    sol = __bisect(err_function, interval)

    if target_value < liq_value:

        v = sol

    else:

        v = 1 / sol

    return fluid.get_state(p=p, t=t, v=v)


def __bisect(err_function, interval):

    n = 0
    curr_value = np.mean(interval)

    while n < 25:

        curr_value = np.mean(interval)
        error = err_function(curr_value)
        n += 1

        if error < 0:

            interval[1] = curr_value

        else:

            interval[0] = curr_value

    return curr_value
