from scipy.optimize import root_scalar
import numpy as np

COEFFICIENTS = {

    "RK": [

        1.49803361e-16,  5.77595975e-14,  9.70857882e-12,
        9.35085484e-10,  5.71092803e-08,  2.31578922e-06,
        6.37056724e-05, 1.20490603e-03,  1.61578963e-02,
        1.79115746e-01, -4.21340902e-04

    ]

}

def t_sat_rel(p_rel, cs):

    if p_rel > 0:

        x = np.log(p_rel)
        y = 0.

        for c in cs:

            y = y * x + c

        return np.exp(y)

    else:

        return -np.inf

def p_sat_rel(t_rel, cs):

    def zero_function(p_rel_guess):

        if p_rel_guess > 0:

            values = t_sat_rel(p_rel_guess, cs) - t_rel

        else:

            values = -np.inf

        return values

    sol = root_scalar(zero_function, bracket=[0, 1], fprime2=False, rtol=1E-3)
    print(sol)

def __get_der_cs(cs):

    ncs = len(cs)
    cs_1 = np.zeros(ncs - 1)

    for i in range(len(cs_1)):
        cs_1[i] = (ncs - 1 - i) * cs[i]

    return cs_1

