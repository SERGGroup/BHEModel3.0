# %%-------------------------------------   IMPORT MODULES                      -------------------------------------> #
from REFPROPConnector import ThermodynamicPoint
import matplotlib.pyplot as plt
from itertools import cycle
import scipy.constants
import pandas as pd
import numpy as np


# %%-------------------------------------   FUNCTION DEFINITION                 -------------------------------------> #

class RKFluid:

    def __init__(self, p_crit, t_crit, cp_ideal, m_molar):

        self.r_spc = scipy.constants.R / m_molar
        self.p_crit = p_crit
        self.t_crit = t_crit
        self.v_crit = self.r_spc * self.t_crit / (3 * self.p_crit)

        self.cp_ideal = cp_ideal
        self.m_molar = m_molar

        alpha = self.r_spc * self.t_crit
        beta = np.power(2, 1/3) - 1

        self.a = alpha ** 2 * np.sqrt(self.t_crit) / self.p_crit / (9 * beta)
        self.b = beta * alpha / (3  * self.p_crit)

    def p(self, t, v, check_sat=True):

        if check_sat and t < self.t_crit:

            p_sat, z_l, z_v = self.p_sat(t)

            v_liq = z_l * self.r_spc * t / p_sat
            v_vap = z_v * self.r_spc * t / p_sat

            if v_liq < v < v_vap:

                return p_sat

        p = self.r_spc * t / (v - self.b) - self.a / (np.sqrt(t) * v * (v + self.b))
        return p

    def p_sat(self, t):

        if t < self.t_crit:

            p_rels = [0, 1]
            n = 0

            while n < 25:

                p_rel = np.mean(p_rels)
                p_curr = p_rel * self.p_crit

                err, z_l, z_v = self.error_p(t, p_curr)

                if err < 0:

                    p_rels[1] = p_rel

                else:

                    p_rels[0] = p_rel

                n += 1

            return p_curr, z_l, z_v

        return np.nan, np.nan, np.nan

    def f(self, t, v):

        p = self.p(t, v)
        z = p * v / (self.r_spc * t)

        return np.exp(self.f_z(t, z))

    def f_z(self, t, z):

        a = z - 1
        b = np.log(z - self.b)
        c = self.a / (self.b * np.sqrt(t)) * np.log(z/(z + self.b))

        return a - b + c

    def error_p(self, t, p):

        A = self.a * p / (np.sqrt(t)*self.r_spc ** 2 * t ** 2)
        B = self.b * p / (self.r_spc * t)
        res = np.roots([1, -1, A - B - B**2, -A*B])

        z_l = np.min(np.real(res))
        z_v = np.max(np.real(res))

        f_l = self.f_z(t,np.min(np.real(res)))
        f_v = self.f_z(t,np.max(np.real(res)))

        f_l_1 = z_l - 1 - np.log(z_l - B) - A / B * np.log((z_l + B) / z_l)
        f_v_1 = z_v - 1 - np.log(z_v - B) - A / B * np.log((z_v + B) / z_v)

        return f_l_1 - f_v_1, z_l, z_v

    def error(self, t, v_gas):

        v_g = v_gas / self.b
        alpha = self.a / (np.sqrt(t) * (self.b * self.r_spc * t))
        beta = 1 / (v_g - 1) - alpha * (1 / (v_g * (v_g + 1)))
        p = self.p(t, v_gas)
        beta_1 = p * self.b / (self.r_spc * t)

        v_l, v_m, v_g = self.__v_l(alpha, beta, v_g)

        if v_l > 0:

            z_l = beta * v_l
            f_l = beta * v_l - 1 - np.log(beta*(v_l - 1)) - alpha * self.r_spc * t * np.log((v_l + 1)/v_l)
            f_g = beta * v_g - 1 - np.log(beta*(v_g - 1)) - alpha * self.r_spc * t * np.log((v_g + 1)/v_g)

            p_l = self.p(t, v_l * self.b)
            p_g = self.p(t, v_g * self.b)

            return (f_l - f_g)/f_g, np.nan

        else:

            return np.nan, np.nan

    def __v_l(self, alpha, beta, v_g):

        a = beta
        b = beta * v_g - 1
        c = b * v_g + alpha - beta - 1

        b = b/2
        delta = b ** 2 - a * c

        if delta >= 0:

            delta_sqr = np.sqrt(delta)
            x1 = (delta_sqr - b) / a
            x2 = -(delta_sqr + b) / a

            arr = np.array([x1, x2, v_g])
            arr.sort()

            return arr

        else:

            return - np.inf, - np.inf, v_g


# %%-------------------------------------   CALCULATE POINTS                    -------------------------------------> #
rp_co2 = ThermodynamicPoint(["Carbon Dioxide"], [1], unit_system="MASS BASE SI")
t_crit = rp_co2.RPHandler.TC
p_crit = rp_co2.RPHandler.PC

rp_co2.set_variable("T", t_crit)
rp_co2.set_variable("P", p_crit)
v_crit = 1 / rp_co2.get_variable("rho")
cp_ideal = rp_co2.get_variable("CP0")
m_molar = 0.04401

rk_co2 = RKFluid(p_crit, t_crit, cp_ideal, m_molar)

n_points = 50
v_rels = np.logspace(-1, 4, n_points)
styles = cycle(["--", "-.", "-", ":"])

for t_mod in [0.5, 0.7, 0.8, 0.9, 1, 2]:

    p_res = np.zeros((2, n_points))
    f_res = np.zeros((1, n_points))

    for i in range(len(v_rels)):

        v_curr = v_rels[i] * (rk_co2.v_crit - rk_co2.b) + rk_co2.b
        t_curr = t_crit * t_mod
        rp_co2.set_variable("t", t_curr)
        rp_co2.set_variable("rho", 1/v_curr)

        p_res[0, i] = rp_co2.get_variable("P")
        p_res[1, i] = rk_co2.p(t_curr, v_curr)

        f_res[0, i] = rk_co2.f(t_curr, v_curr)

    style = next(styles)
    plt.plot(v_rels, p_res[0,:], label="real", linestyle=style, color="tab:blue")
    plt.plot(v_rels, p_res[1,:], label="RK", linestyle=style, color="tab:orange")

plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()


# %%-------------------------------------   PLOT DENSITY RESULTS                -------------------------------------> #
n_points = 10000

rho_rels = np.logspace(-6, 10, n_points)
v_rels = 1 / rho_rels
t_curr = rk_co2.t_crit * 0.7

values = np.zeros((4, n_points))

for i in range(len(v_rels)):

    v_liq = v_rels[i] * (rk_co2.v_crit - rk_co2.b) + rk_co2.b
    values[0, i] = rk_co2.p(t_curr, v_liq) / rk_co2.p_crit
    values[1, i] = v_liq

    v_l, v_l2 = rk_co2.error(t_curr, v_liq)
    values[2, i] = v_l
    values[3, i] = v_l2

df = pd.DataFrame(values[0:3, :].T)
df = df.dropna(axis=0)
x_values = np.log10(df[0])

plt.scatter(df[0], df[2])
plt.grid(visible=True, which="major", axis="both")
plt.xscale("log")
plt.show()


# %%-------------------------------------   PLOT PRESSURE RESULTS               -------------------------------------> #
n_points = 100

p_rels = np.logspace(0, -2, n_points)
t_curr = rk_co2.t_crit * 0.9

values = np.zeros((4, n_points))

for i in range(len(p_rels)):

    p_curr = p_rels[i] * rk_co2.p_crit
    values[0, i] = p_curr / rk_co2.p_crit

    error = rk_co2.error_p(t_curr, p_curr)
    values[2, i] = error

df = pd.DataFrame(values[0:3, :].T)
df = df.dropna(axis=0)
x_values = np.log10(df[0])

plt.scatter(df[0], df[2])
plt.grid(visible=True, which="major", axis="both")
plt.xscale("log")
plt.show()
