from main_classes.eos.cubic_eos import CubicEOS
import numpy as np

ACCEPTED_CODES = ["PT", "PV", "VT", "QT", "QP", "HT", "HP", "SP", "ST"]
VARIABLE_CODES = ["P", "T", "H", "S", "Q", "V"]


def calculate_flash(fluid: CubicEOS, input_code, input_1, input_2):

    __check_input_code(input_code)
    extr_values = __extract_values(input_code, input_1, input_2)

    if extr_values["H"] is not None:

        return iteration_flash(fluid, extr_values["H"], t=extr_values["T"], p=extr_values["P"], search_h=True)

    elif extr_values["S"] is not None:

        return iteration_flash(fluid, extr_values["S"], t=extr_values["T"], p=extr_values["P"], search_h=False)

    elif extr_values["Q"] is not None:

        return base_quality_flash(fluid, q=extr_values["Q"], t=extr_values["T"], p=extr_values["P"])

    else:

        return direct_flash(fluid, t=extr_values["T"], p=extr_values["P"], v=extr_values["V"])


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
            return state.h - target_value

    else:

        liq_value = liq_state.s
        vap_value = vap_state.s

        def sub_err_function(v):

            state = fluid.get_state(p=p, t=t, v=v)
            return state.s - target_value

    if liq_value < target_value < vap_value:

        q = (target_value - liq_value) / (vap_value - liq_value)
        v = (vap_state.v - liq_state.v) * q + liq_state.v
        return fluid.get_state(p=p, t=t, v=v)

    if target_value < liq_value:

        interval = [0, liq_state.v]
        err_function = lambda v: sub_err_function(v=v)

    else:

        interval = [0, 1/vap_state.v]
        err_function = lambda rho: sub_err_function(v=1/rho)

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


def __extract_values(input_code, input_1, input_2):

    inputs = [input_1, input_2]
    return_dict = dict()

    for var in VARIABLE_CODES:
        return_dict.update({var: None})

    for i in range(2):

        for var in return_dict.keys():

            if var in input_code[i]:
                return_dict.update({var: inputs[i]})
                break

    return return_dict


def __check_input_code(input_code):

    reversed_code = input_code[::-1]

    if not len(input_code) == 2:
        __raise_error(input_code)

    if (input_code not in ACCEPTED_CODES) and (reversed_code not in ACCEPTED_CODES):
        __raise_error(input_code)


def __raise_error(input_code):

    raise Exception(

        "\n\n'{INPUT_CODE}' is not an acceptable input code\nAccepted codes are: {CODES_LIST}".format(

            INPUT_CODE=input_code,
            CODES_LIST=__accepted_code_str()

        )

    )


def __accepted_code_str():
    accepted_str = "\n\n" + "\033[1m"

    for var in VARIABLE_CODES:
        accepted_str += "\t" + var

    accepted_str += "\033[0m" + "\n"

    for var_1 in VARIABLE_CODES:

        accepted_str += "\033[1m" + var_1 + "\033[0m"

        for var_2 in VARIABLE_CODES:

            if (var_1 + var_2) in ACCEPTED_CODES or (var_2 + var_1) in ACCEPTED_CODES:

                accepted_str += "\t" + str(u'\u2714')

            else:

                accepted_str += "\t" + str(u'\u2718')

        accepted_str += "\n"

    accepted_str += "\n" + "For example, 'PT' and 'HT' are acceptable but 'HS' or 'PP' are not!"
    return accepted_str
