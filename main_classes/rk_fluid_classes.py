from scipy.integrate import RK45
import scipy.constants
import numpy as np


class RKFluid:

    def __init__(self, p_crit, t_crit, cp_ideal, m_molar):

        self.r_spc = scipy.constants.R / m_molar
        self.p_crit = p_crit
        self.t_crit = t_crit
        self.v_crit = self.r_spc * self.t_crit / (3 * self.p_crit)

        self.cp_ideal = cp_ideal
        self.m_molar = m_molar

        alpha = self.r_spc * self.t_crit
        beta = np.power(2, 1 / 3) - 1

        self.a = alpha ** 2 * np.sqrt(self.t_crit) / self.p_crit / (9 * beta)
        self.b = beta * alpha / (3 * self.p_crit)

    def get_state(self, p=None, t=None, v=None, check_sat=True):

        if (t is not None) and (v is not None):

            return FluidState(self, t, v)

        elif (t is not None) and (p is not None):

            return FluidState(self, t, self.v(t, p))

        elif (p is not None) and (v is not None):

            return FluidState(self, self.t(v, p, check_sat=check_sat), v)

        return FluidState(self)

    def get_sat_state(self, p=None, t=None, liquid=True):

        if (t is None) and (p is None):

            return FluidState(self)

        if t is not None:

            t_sat = t
            p_sat, z_l, z_v = self.p_sat(t)

        else:

            p_sat = p
            t_sat, z_l, z_v = self.t_sat(p)

        if liquid:

            v_sat = (z_l * self.r_spc * t_sat) / p_sat

        else:

            v_sat = (z_v * self.r_spc * t_sat) / p_sat

        return FluidState(self, t_sat, v_sat)

    def p(self, t, v, check_sat=True):

        if check_sat and t < self.t_crit:

            p_sat, z_l, z_v = self.p_sat(t)

            v_liq = z_l * self.r_spc * t / p_sat
            v_vap = z_v * self.r_spc * t / p_sat

            if v_liq < v < v_vap:

                return p_sat

        z_l, z_v = self.z(t=t, v=v)
        return z_l * self.r_spc * t / v

    def t(self, v, p, check_sat=True):

        if check_sat and p < self.p_crit:

            t_sat, z_l, z_v = self.t_sat(p)
            v_liq = z_l * self.r_spc * t_sat / p
            v_vap = z_v * self.r_spc * t_sat / p

            if v_liq < v < v_vap:
                return t_sat

        z_l, z_v = self.z(v=v, p=p)
        return p * v / (self.r_spc * z_l)

    def v(self, t, p):

        z_l, z_v = self.z(t=t, p=p)

        if not z_l == z_v:

            p_sat, z_l_sat, z_v_sat = self.p_sat(t)

            if p > p_sat:

                z = z_l

            else:

                z = z_v

        else:

            z = z_l

        return z * self.r_spc * t / p

    def p_sat(self, t):

        if t < self.t_crit:

            p_rels = [0, 1]
            p_curr, z_l, z_v = np.zeros(3)

            n = 0
            while n < 25:

                p_rel = np.mean(p_rels)
                p_curr = p_rel * self.p_crit

                err, z_l, z_v = self.__error_sat(t, p_curr)

                if err < 0:

                    p_rels[1] = p_rel

                else:

                    p_rels[0] = p_rel

                n += 1

            return p_curr, z_l, z_v

        return np.nan, np.nan, np.nan

    def t_sat(self, p):

        if p < self.p_crit:

            t_rels = [0, 1]
            t_curr, z_l, z_v = np.zeros(3)

            n = 0
            while n < 25:

                t_rel = np.mean(t_rels)
                t_curr = t_rel * self.t_crit

                err, z_l, z_v = self.__error_sat(t_curr, p)

                if err < 0:

                    t_rels[1] = t_rel

                else:

                    t_rels[0] = t_rel

                n += 1

            return t_curr, z_l, z_v

        elif p == self.p_crit:

            return self.t_crit, self.v_crit, self.v_crit

        return np.nan, np.nan, np.nan

    def z(self, p=None, t=None, v=None):

        if (t is not None) and (v is not None):

            p = self.r_spc * t / (v - self.b) - self.a / (np.sqrt(t) * v * (v + self.b))
            z_l = p * v / (self.r_spc * t)
            z_v = z_l

        elif (t is not None) and (p is not None):

            A = self.a * p / (np.sqrt(t) * self.r_spc ** 2 * t ** 2)
            B = self.b * p / (self.r_spc * t)

            res = np.roots([1, -1, A - B - B ** 2, -A * B])
            z_l, z_v = self.__get_real_res(res)

        elif (p is not None) and (v is not None):

            A = self.a / (v * (v + self.b))
            B = self.r_spc / (v - self.b)

            res = np.roots([B ** 2, -2 * B * p, p ** 2, - A ** 2])
            t_l, t_v = self.__get_real_res(res)
            (z_l, z_v) = p * v / (self.r_spc * np.array((t_l, t_v)))

        else:

            z_l = np.nan
            z_v = np.nan

        return z_l, z_v

    @staticmethod
    def __get_real_res(res):

        res_real = np.isreal(res)
        if np.sum(res_real) == 1:

            res_1 = np.real(np.sum(res, where=res_real))
            res_2 = res_1

        else:

            res_1 = np.min(np.real(res))
            res_2 = np.max(np.real(res))

        return res_1, res_2

    def __error_sat(self, t, p):

        A = self.a * p / (np.sqrt(t) * self.r_spc ** 2 * t ** 2)
        B = self.b * p / (self.r_spc * t)
        res = np.roots([1, -1, A - B - B ** 2, -A * B])

        z_l = np.min(np.real(res))
        z_v = np.max(np.real(res))

        f_l = z_l - 1 - np.log(z_l - B) - A / B * np.log((z_l + B) / z_l)
        f_v = z_v - 1 - np.log(z_v - B) - A / B * np.log((z_v + B) / z_v)

        return f_l - f_v, z_l, z_v


class FluidState:

    def __init__(self, fluid_solver: RKFluid, t=None, v=None, check_sat=True):

        self.fluid_solver = fluid_solver

        self.__t = t
        self.__v = v
        self.__check_sat = check_sat

        self.__p = None
        self.__s = None
        self.__h = None

        self.__r = None
        self.__cp = None
        self.__dpdt = None
        self.__dpdv = None
        self.__int_cp = None

        self.__init_saturation_condition()

    def __init_saturation_condition(self):

        self.__x = None
        self.__liquid_state = None
        self.__vapour_state = None

        if self.__check_sat:

            if self.__t < self.fluid_solver.t_crit:

                p_sat, z_l, z_v = self.fluid_solver.p_sat(self.__t)
                v_liq = z_l * self.fluid_solver.r_spc * self.__t / p_sat
                v_vap = z_v * self.fluid_solver.r_spc * self.__t / p_sat

                if v_liq < self.__v < v_vap:

                    self.__p = p_sat
                    self.__x = (self.__v - v_liq) / (v_vap - v_liq)
                    self.__liquid_state = FluidState(self.fluid_solver, t=self.__t, v=v_liq, check_sat=False)
                    self.__vapour_state = FluidState(self.fluid_solver, t=self.__t, v=v_vap, check_sat=False)

    def update_state(self, t=None, v=None, check_sat=True):

        self.__t = t
        self.__v = v
        self.__check_sat = check_sat

        self.__p = None
        self.__r = None
        self.__cp = None
        self.__dpdt = None
        self.__dpdv = None
        self.__int_cp = None

    @property
    def t(self):
        return self.__t

    @property
    def v(self):

        return self.__v

    @property
    def p(self):

        if self.__p is None:

            self.__p = self.fluid_solver.p(self.__t, self.__v, self.__check_sat)

        return self.__p

    @property
    def h(self):

        return 0.

    @property
    def s(self):

        return 0.

    @property
    def r(self):

        if self.__r is None:

            self.__r = - self.__t * self.dpdt ** 2 / self.dpdv

        return self.__r

    @property
    def cp(self):

        if self.bifase:

            return self.__liquid_state.cp * (1 - self.__x) + self.__vapour_state.cp * self.__x

        if self.__cp is None:

            self.__cp = self.fluid_solver.cp_ideal + self.r - self.fluid_solver.r_spc + self.int_cp

        return self.__cp

    @property
    def dpdt(self):

        if self.__dpdt is None:
            v = self.__v
            t = self.__t
            b = self.fluid_solver.b
            a = self.fluid_solver.a
            r = self.fluid_solver.r_spc

            self.__dpdt = r / (v - b) + 1 / 2 * a * (np.power(t, -3 / 2) / (v * (v + b)))

        return self.__dpdt

    @property
    def dpdv(self):

        if self.__dpdv is None:
            v = self.__v
            t = self.__t
            r = self.fluid_solver.r_spc
            b = self.fluid_solver.b
            a = self.fluid_solver.a

            self.__dpdv = a / np.sqrt(t) * (b + 2 * v) / (v * (v + b)) ** 2 - r * t / (v - b) ** 2

        return self.__dpdv

    @property
    def int_cp(self):

        if self.__int_cp is None:

            v = self.__v
            t = self.__t
            a = self.fluid_solver.a
            b = self.fluid_solver.b

            if v <= 0:

                self.__int_cp = np.inf
                return np.inf

            self.__int_cp = 3 / 4 * a / b * np.log((v + b) / b) / np.power(t, 3 / 2)

        return self.__int_cp

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


def calculate_vertical(fluid_eos, state_in, res_depth, downward=True, check_sat=False):
    def rk_overall_der(z, y):

        t_curr = y[0]
        v_curr = y[1]

        state = fluid_eos.get_state(t=t_curr, v=v_curr, check_sat=check_sat)

        dp = - scipy.constants.g / v_curr
        dt = state.r / state.cp / state.dpdt * dp
        dv = (1 + state.r / state.cp) / state.dpdv * dp

        return [dt, dv]

    if downward:

        res_depth = -abs(res_depth)

    else:

        res_depth = abs(res_depth)

    integrator = RK45(rk_overall_der, 0, [state_in.t, state_in.v], res_depth, rtol=1e-06, atol=1e-07)

    while integrator.status == 'running':

        integrator.step()

    output = integrator.y

    return fluid_eos.get_state(t=output[0], v=output[1], check_sat=False)


def evaluate_system(fluid_eos, in_state, depth_res, t_res):

    res_in = calculate_vertical(fluid_eos, in_state, depth_res, downward=True)
    res_out = fluid_eos.get_state(p=res_in.p, t=t_res, check_sat=False)
    sys_out = calculate_vertical(fluid_eos, res_out, depth_res, downward=False, check_sat=True)

    return [in_state, res_in, res_out, sys_out]
