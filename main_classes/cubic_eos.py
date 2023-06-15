from main_classes.support.p_sat_calculation import p_sat_rel, t_sat_rel
from scipy.optimize import fsolve, root_scalar
from abc import ABC, abstractmethod
import scipy.constants
import numpy as np


def get_real_res(res):

    res_real = np.isreal(res)
    if np.sum(res_real) == 1:

        res_1 = np.real(np.sum(res, where=res_real))
        res_2 = res_1

    else:

        res_1 = np.min(np.real(res))
        res_2 = np.max(np.real(res))

    return res_1, res_2


class CubicEOS(ABC):

    def __init__(self, p_crit=7380000., t_crit=304.1, cp_ideal=845.85, m_molar=0.04401, acntr=0.239):

        """
        :param p_crit: critical pressure [Pa]
        :param t_crit: critical temperature [K]
        :param cp_ideal: ideal gas specific heat capacity [J/kg K]
        :param m_molar: molar mass [kg/mol]
        :param acntr: acentricity factor (omega) - if required by the EoS [-]
        """

        self.r_spc = scipy.constants.R / m_molar

        self.p_crit = p_crit
        self.t_crit = t_crit
        self.cp_ideal = cp_ideal
        self.m_molar = m_molar
        self.acntr = acntr

        self.__init_coefficients()

    def __init_coefficients(self):

        self.a_vdw = 27 / 64 * (self.r_spc * self.t_crit) ** 2 / self.p_crit
        self.b_vdw = self.r_spc * self.t_crit / (8 * self.p_crit)

        # Default values to be updated in init_coefficient()
        self.a_0 = self.a_vdw
        self.b = self.b_vdw
        self.z_crit = 1/3
        self.r_1 = 0
        self.r_2 = -1

        self.init_coefficient()

        self.v_crit = self.z_crit * self.r_spc * self.t_crit / self.p_crit

        self.r12 = self.r_1 * self.r_2
        self.r1r2 = self.r_1 + self.r_2

    def get_state(self, p=None, t=None, v=None):

        if (t is not None) and (v is not None):

            return FluidState(self, t, v)

        elif (t is not None) and (p is not None):

            return FluidState(self, t, self.v(t, p))

        elif (p is not None) and (v is not None):

            return FluidState(self, self.t(v, p), v)

        return FluidState(self)

    def get_sat_state(self, p=None, t=None, which="liq"):

        if (t is None) and (p is None):

            return FluidState(self)

        if t is not None:

            t_sat = t
            p_sat, z_l, z_v = self.p_sat(t)

        else:

            p_sat = p
            t_sat, z_l, z_v = self.t_sat(p)

        if which.lower() == "both":

            v_l = (z_l * self.r_spc * t_sat) / p_sat
            v_v = (z_v * self.r_spc * t_sat) / p_sat

            return FluidState(self, t_sat, v_l), FluidState(self, t_sat, v_v)

        else:

            if which.lower() == "liq" or "liquid" or "l":

                v_sat = (z_l * self.r_spc * t_sat) / p_sat

            else:

                v_sat = (z_v * self.r_spc * t_sat) / p_sat

            return FluidState(self, t_sat, v_sat)

    def p(self, t, v):

        if t < self.t_crit:

            p_sat, z_l, z_v = self.p_sat(t, get_approximate=True)

            if p_sat is None:

                p_sat, z_l, z_v = self.p_sat(t,get_approximate=False)
                v_liq = z_l * self.r_spc * t / p_sat
                v_vap = z_v * self.r_spc * t / p_sat

                if v_liq < v < v_vap:

                    return p_sat

            else:

                v_liq = z_l * self.r_spc * t / p_sat
                v_vap = z_v * self.r_spc * t / p_sat

                if v_liq * 0.98 < v < v_vap * 1.02:

                    p_sat, z_l, z_v = self.p_sat(t, get_approximate=False)
                    v_liq = z_l * self.r_spc * t / p_sat
                    v_vap = z_v * self.r_spc * t / p_sat

                    if v_liq < v < v_vap:

                        return p_sat

        z_l, z_v = self.z(t=t, v=v)
        return z_l * self.r_spc * t / v

    def t(self, v, p):

        if p < self.p_crit:

            t_sat, z_l, z_v = self.t_sat(p, get_approximate=True)

            if t_sat is None:

                t_sat, z_l, z_v = self.t_sat(p, get_approximate=False)
                v_liq = z_l * self.r_spc * t_sat / p
                v_vap = z_v * self.r_spc * t_sat / p

                if v_liq < v < v_vap:
                    return t_sat

            else:

                v_liq = z_l * self.r_spc * t_sat / p
                v_vap = z_v * self.r_spc * t_sat / p

                if v_liq * 0.98 < v < v_vap * 1.02:

                    t_sat, z_l, z_v = self.t_sat(p, get_approximate=False)
                    v_liq = z_l * self.r_spc * t_sat / p
                    v_vap = z_v * self.r_spc * t_sat / p

                    if v_liq < v < v_vap:
                        return t_sat

        z_l, z_v = self.z(v=v, p=p)
        return p * v / (self.r_spc * z_l)

    def v(self, t, p):

        z_l, z_v = self.z(t=t, p=p)

        f_l = self.fug(t, p, z_l)
        f_v = self.fug(t, p, z_v)

        if f_v > f_l:

            z = z_l

        else:

            z = z_v

        return z * self.r_spc * t / p

    def z(self, p=None, t=None, v=None):

        if (t is not None) and (v is not None):

            p = self.r_spc * t / (v - self.b) - self.a(t) / ((v - self.b * self.r_1) * (v - self.b * self.r_2))
            z_l = p * v / (self.r_spc * t)
            z_v = z_l

        elif (t is not None) and (p is not None):

            alpha = self.a(t) / (self.b * self.r_spc * t)
            beta = self.b * p / (self.r_spc * t)

            A_1 = beta * (self.r1r2 + 1) + 1
            A_2 = beta * (beta * self.r12 + alpha + self.r1r2 * (beta + 1))
            A_3 = beta ** 2 * (self.r12 * (beta + 1) + alpha)

            res = np.roots([1, -A_1, A_2, -A_3])
            z_l, z_v = get_real_res(res)

        elif (p is not None) and (v is not None):

            t = self.iterate_t(p, v)
            z_l = p * v / (self.r_spc * t)
            z_v = z_l

        else:

            z_l = np.nan
            z_v = np.nan

        return z_l, z_v

    def p_sat(self, t, get_approximate=False):

        return self.__get_sat(t, get_approximate=get_approximate, get_p=True)

    def t_sat(self, p, get_approximate=False):

        return self.__get_sat(p, get_approximate=get_approximate, get_p=False)

    def __get_sat(self, y, get_approximate=False, get_p=True):

        if get_approximate:

            # return self.__approx_sat(y, get_p=get_p)
            return self.__iterate_sat(y, get_p=get_p)

        else:

            return self.__iterate_sat(y, get_p=get_p)

    def __iterate_sat(self, y, get_p):

        if get_p:

            x_crit = self.p_crit
            y_crit = self.t_crit

        else:

            x_crit = self.t_crit
            y_crit = self.p_crit

        if y < y_crit:

            x_rels = [0, 1]
            x_curr, z_l, z_v = np.zeros(3)

            n = 0
            while n < 25:

                x_rel = np.mean(x_rels)
                x_curr = x_rel * x_crit

                if get_p:

                    err, z_l, z_v = self.__error_sat(t=y, p=x_curr)

                else:

                    err, z_l, z_v = self.__error_sat(t=x_curr, p=y)
                    err = -err

                if err < 0:

                    x_rels[1] = x_rel

                else:

                    x_rels[0] = x_rel

                n += 1

            return x_curr, z_l, z_v

        elif y == y_crit:

            return x_crit, self.z_crit, self.z_crit

        return np.nan, np.nan, np.nan

    def __approx_sat(self, y, get_p):

        cs = self.sat_coefficients

        if self.sat_coefficients is not None:

            if get_p:

                x = p_sat_rel(y / self.t_crit, cs) * self.p_crit
                z_l, z_v = self.z(t=y, p=x)

            else:

                x = t_sat_rel(y / self.p_crit, cs) * self.t_crit
                z_l, z_v = self.z(t=x, p=y)

            return x, z_l, z_v

        else:

            return None, None, None

    def __error_sat(self, t, p):

        z_l, z_v = self.z(t=t, p=p)

        if not z_l == z_v:

            f_l = self.fug(t, p, z_l)
            f_v = self.fug(t, p, z_v)

            return f_l - f_v, z_l, z_v

        else:

            v = z_l * self.r_spc * t / p
            error = -(self.a(t) / (v - self.b) - self.a(self.t_crit) / (self.v_crit - self.b))

            return error, z_l, z_v

    def ddp_ddv(self, t, p, z):

        v = z * self.r_spc * t / p
        alpha = self.a(t) / (self.b * self.r_spc * t)
        gamma = v / self.b

        a_0 = 2 * self.r_spc * t / self.b ** 3
        a_1 = 1 / (gamma - 1) ** 3
        a_2 = (self.r_1 ** 2 + self.r_2 ** 2 + 3 * self.r12) - gamma * self.r1r2 + 3 * gamma ** 2
        a_3 = (gamma - self.r_1) ** 3 * (gamma - self.r_2) ** 3

        return a_0 * (a_1 - a_2 / a_3 * alpha)

    def fug(self, t, p, z):

        v = z * self.r_spc * t / p
        alpha = self.a(t) / (self.b * self.r_spc * t)
        beta = self.b * p / (self.r_spc * t)

        f_v = (np.log(v - self.b * self.r_1) - np.log(v - self.b * self.r_2)) / (self.r_1 - self.r_2)
        return z - 1 + alpha * f_v - np.log(z - beta)

    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #
    #                                   FUNCTION TO BE UPDATED DEPENDING ON THE EOS
    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #

    @abstractmethod
    def a(self, t):

        return self.a_0 / np.sqrt(t)

    @abstractmethod
    def da(self, t):

        return -self.a_0 / (2 * np.power(t, 3/2))

    @abstractmethod
    def dda(self, t):

        return 3 * self.a_0 / (4 * np.power(t, 5/2))

    @abstractmethod
    def init_coefficient(self):

        alpha = self.r_spc * self.t_crit
        beta = np.power(2, 1 / 3) - 1

        self.a_0 = alpha ** 2 * np.sqrt(self.t_crit) / self.p_crit / (9 * beta)
        self.b = beta * alpha / (3 * self.p_crit)

        self.z_crit = 1/3
        self.r_1 = 0
        self.r_2 = -1

    @property
    def sat_coefficients(self):

        return None

    def iterate_t(self, p, v):

        def root_funct(t_rel):

            t = t_rel * self.t_crit
            return p - self.r_spc * t / (v - self.b) + self.a(t) / (v * (v + self.b))

        A = self.r_spc / (v - self.b_vdw)
        B = self.a_vdw / (v * (v + self.b_vdw))
        t_0 = (p + B) / A

        sol = root_scalar(root_funct, x0=t_0 / self.t_crit, x1=1)
        return sol.root * self.t_crit

    def iterate_coefficients(self, x_0):

        self.r12 = self.r_1 * self.r_2
        self.r1r2 = self.r_1 + self.r_2

        beta_1 = self.p_crit / (self.r_spc * self.t_crit)
        alpha_1 = self.a(self.t_crit) / (self.r_spc * self.t_crit) / self.a_0
        a_3 = self.r_1 ** 2 + self.r_2 ** 2 + 3 * self.r_1 * self.r_2

        def f_iter(x):

            if len(x) > 2:

                self.v_crit, self.a_0, self.b = x

            else:

                self.a_0, self.b = x

            alpha = alpha_1 * self.a_0 / self.b
            beta = beta_1 * self.b
            eta = self.v_crit / self.b

            a1 = (eta - 1)
            a2 = (eta - self.r_1) * (eta - self.r_2)

            p = 1 / a1 - alpha / a2 - beta
            dpdv = (2 * eta + self.r1r2) * alpha / a2 ** 2 - 1 / a1 ** 2
            ddpddv = 1 / a1 ** 3 - (3 * eta ** 2 - eta * self.r1r2 + a_3) * alpha / a2 ** 3

            if len(x) > 2:

                return [p, dpdv, ddpddv]

            else:

                return [dpdv, ddpddv]

        if len(x_0) > 2:
            self.v_crit, self.a_0, self.b = fsolve(f_iter, x_0)

        else:
            self.a_0, self.b = fsolve(f_iter, x_0)


class FluidState:

    def __init__(self, fluid_solver: CubicEOS, t=None, v=None):

        self.fluid_solver = fluid_solver

        self.__t = t
        self.__v = v

        self.__init_parameters()

    def __init_parameters(self):

        # State Variables
        self.__p = None
        self.__s = None
        self.__h = None
        self.__a = None
        self.__g = None

        # Additional State Variables
        self.__r = None
        self.__cp = None

        # Derivatives
        self.__dpdt = None
        self.__dpdv = None
        self.__ddpddt = None
        self.__ddpddv = None

        # Integrals
        self.__int_t_ddpddt = None
        self.__int_t_dpdt = None
        self.__int_rv = None
        self.__int_rtv = None

        self.__init_saturation_condition()

        # Recurring Parameters
        self.__u = None
        self.__g = None
        self.__alpha = None
        self.__beta = None
        self.__eta = None

    def __init_saturation_condition(self):

        self.__x = None
        self.__liquid_state = None
        self.__vapour_state = None

        if self.__t < self.fluid_solver.t_crit:

            p_sat, z_l, z_v = self.fluid_solver.p_sat(self.__t)
            v_liq = z_l * self.fluid_solver.r_spc * self.__t / p_sat
            v_vap = z_v * self.fluid_solver.r_spc * self.__t / p_sat

            if v_liq < self.__v < v_vap:

                self.__p = p_sat
                self.__x = (self.__v - v_liq) / (v_vap - v_liq)
                self.__liquid_state = FluidState(self.fluid_solver, t=self.__t, v=v_liq)
                self.__vapour_state = FluidState(self.fluid_solver, t=self.__t, v=v_vap)

    def update_state(self, t=None, v=None):

        self.__t = t
        self.__v = v

        self.__init_parameters()

    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #
    #                                                 PROPERTIES
    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #

    @property
    def t(self):

        return self.__t

    @property
    def v(self):

        return self.__v

    @property
    def p(self):

        if self.__p is None:

            self.__p = self.fluid_solver.p(self.__t, self.__v)

        return self.__p

    @property
    def h(self):

        if self.__h is None:

            if self.bifase:

                self.__h = self.__liquid_state.h * (1 - self.__x) + self.__vapour_state.h * self.__x

            else:

                h_ideal = self.fluid_solver.cp_ideal * self.__t
                h_dep = self.int_t_dpdt + self.p * self.__v - self.fluid_solver.r_spc * self.__t
                self.__h = h_ideal + h_dep

        return self.__h

    @property
    def s(self):

        if self.__s is None:

            if self.bifase:

                self.__s = self.__liquid_state.s * (1 - self.__x) + self.__vapour_state.s * self.__x

            else:

                s_ideal = self.fluid_solver.cp_ideal * np.log(self.__t) - self.fluid_solver.r_spc * np.log(self.p)
                s_dep = - self.int_rv
                self.__s = s_ideal + s_dep

        return self.__s

    @property
    def r(self):

        if self.__r is None:

            self.__r = - self.__t * self.dpdt ** 2 / self.dpdv

        return self.__r

    @property
    def cp(self):

        if self.bifase:

            return np.inf

        if self.__cp is None:

            self.__cp = self.fluid_solver.cp_ideal + self.r - self.fluid_solver.r_spc + self.int_t_ddpddt

        return self.__cp

    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #
    #                                                 DERIVATIVES
    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #

    @property
    def dpdt(self):

        if self.__dpdt is None:

            t = self.__t
            b = self.fluid_solver.b
            r = self.fluid_solver.r_spc
            da = self.fluid_solver.da(t)

            eta = self.eta
            eta_r1 = eta - self.fluid_solver.r_1
            eta_r2 = eta - self.fluid_solver.r_2

            self.__dpdt = (r / (eta - 1) - da / (b * eta_r1 * eta_r2)) / b

        return self.__dpdt

    @property
    def dpdv(self):

        if self.__dpdv is None:

            eta_r1 = self.eta - self.fluid_solver.r_1
            eta_r2 = self.eta - self.fluid_solver.r_2
            r1r2 = self.fluid_solver.r1r2

            a0 = self.p / self.beta
            a1 = 1 / ((self.eta - 1) ** 2)
            a2 = (2 * self.eta + r1r2) / ((eta_r1 * eta_r2) ** 2)

            self.__dpdv = a0 * (a2 * self.alpha - a1)

        return self.__dpdv

    @property
    def ddpddt(self):

        if self.__ddpddt is None:

            t = self.__t
            dda = self.fluid_solver.dda(t)
            b = self.fluid_solver.b

            eta = self.eta
            eta_r1 = eta - self.fluid_solver.r_1
            eta_r2 = eta - self.fluid_solver.r_2

            self.__ddpddt = - dda / (b ** 2 * eta_r1 * eta_r2)

        return self.__ddpddt

    @property
    def ddpddv(self):

        if self.__ddpddv is None:

            t = self.__t
            b = self.fluid_solver.b
            r = self.fluid_solver.r_spc

            eta = self.eta
            eta_r1 = eta - self.fluid_solver.r_1
            eta_r2 = eta - self.fluid_solver.r_2
            r1r2 = self.fluid_solver.r1r2
            r1r2_sqr = self.fluid_solver.r_1 ** 2 + self.fluid_solver.r_2 ** 2 + 3 * self.fluid_solver.r12

            a_0 = 2 * r * t / b ** 3
            a_1 = 1 / (eta - 1) ** 3
            a_2 = (3 * eta ** 2 - eta * r1r2 + r1r2_sqr) / (eta_r1 * eta_r2) ** 3

            self.__ddpddv = a_0 * (a_1 - a_2 * self.alpha)

        return self.__ddpddv

    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #
    #                                                 INTEGRALS
    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #

    @property
    def int_t_ddpddt(self):

        if self.__int_t_ddpddt is None:

            self.__int_t_ddpddt = - self.__t * self.fluid_solver.dda(self.__t) * self.u

        return self.__int_t_ddpddt

    @property
    def int_t_dpdt(self):

        if self.__int_t_dpdt is None:

            t = self.__t
            unl_a = self.fluid_solver.a(t) - t * self.fluid_solver.da(t)
            self.__int_t_dpdt = unl_a * self.u

        return self.__int_t_dpdt

    @property
    def int_rv(self):

        if self.__int_rv is None:

            self.__int_rv = self.g + self.fluid_solver.da(self.__t) * self.u

        return self.__int_rv

    @property
    def int_rtv(self):

        if self.__int_rtv is None:

            t = self.__t
            self.__int_rtv = t * self.g + self.fluid_solver.a(t) * self.u

        return self.__int_rtv

    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #
    #                                                RECURRING PARAMETERS
    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #

    @property
    def u(self):

        if self.__u is None:

            eta = self.eta
            eta_r1 = eta - self.fluid_solver.r_1
            eta_r2 = eta - self.fluid_solver.r_2

            try:

                if eta_r1 == eta_r2:

                    self.__u = 1 / (- self.fluid_solver.b * eta_r1)

                else:

                    self.__u = np.log(eta_r1 / eta_r2) / (self.fluid_solver.b * (eta_r2 - eta_r1))

            except:

                self.__u = np.inf

        return self.__u

    @property
    def g(self):

        if self.__g is None:

            self.__g = self.fluid_solver.r_spc * np.log(self.__v / (self.__v - self.fluid_solver.b))

        return self.__g

    @property
    def alpha(self):

        if self.__alpha is None:

            self.__alpha = self.fluid_solver.a(self.__t) / (self.fluid_solver.b * self.fluid_solver.r_spc * self.__t)

        return self.__alpha

    @property
    def beta(self):

        if self.__beta is None:
            self.__beta = self.fluid_solver.b * self.p / (self.fluid_solver.r_spc * self.__t)

        return self.__beta

    @property
    def eta(self):

        if self.__eta is None:
            self.__eta = self.__v / self.fluid_solver.b

        return self.__eta

    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #
    #                                             BIFASE PROPERTIES
    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #

    @property
    def bifase(self):

        return (self.__liquid_state is not None) and (self.__vapour_state is not None) and (self.__x is not None)

    @property
    def liquid_phase(self):

        if self.bifase:

            return self.__liquid_state

        else:

            return self

    @property
    def vapour_phase(self):

        if self.bifase:

            return self.__vapour_state

        else:

            return self
