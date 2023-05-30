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


class CubicEOS:

    def __init__(self, p_crit, t_crit, cp_ideal, m_molar, r_1=0, r_2=-1):

        self.r_spc = scipy.constants.R / m_molar

        self.p_crit = p_crit
        self.t_crit = t_crit
        self.cp_ideal = cp_ideal
        self.m_molar = m_molar

        self.__r_1 = r_1
        self.__r_2 = r_2

        self.__init_coefficients()

    def __init_coefficients(self):

        self.a_vdw = 27 / 64 * (self.r_spc * self.t_crit) ** 2 / self.p_crit
        self.b_vdw = self.r_spc * self.t_crit / (8 * self.p_crit)

        self.a_0 = self.a_vdw
        self.b = self.b_vdw
        self.z_crit = 1/3

        self.init_coefficient()

        self.v_crit = self.z_crit * self.r_spc * self.t_crit / self.p_crit

        self.__r12 = self.__r_1 * self.__r_2
        self.__r1r2 = self.__r_1 + self.__r_2

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

        f_l = self.fug(t, p, z_l)
        f_v = self.fug(t, p, z_v)

        if f_v > f_l:

            z = z_l

        else:

            z = z_v

        return z * self.r_spc * t / p

    def z(self, p=None, t=None, v=None):

        if (t is not None) and (v is not None):

            p = self.r_spc * t / (v - self.b) - self.a(t) / ((v - self.b * self.__r_1) * (v - self.b * self.__r_2))
            z_l = p * v / (self.r_spc * t)
            z_v = z_l

        elif (t is not None) and (p is not None):

            alpha = self.a(t) / (self.b * self.r_spc * t)
            beta = self.b * p / (self.r_spc * t)

            A_1 = beta * (self.__r1r2 + 1) + 1
            A_2 = beta * (beta * self.__r12 + alpha + self.__r1r2 * (beta + 1))
            A_3 = beta ** 2 * (self.__r12 * (beta + 1) + alpha)

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

                if err > 0:

                    t_rels[1] = t_rel

                else:

                    t_rels[0] = t_rel

                n += 1

            return t_curr, z_l, z_v

        elif p == self.p_crit:

            return self.t_crit, self.v_crit, self.v_crit

        return np.nan, np.nan, np.nan

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
        a_2 = (self.__r_1 ** 2 + self.__r_2 ** 2 + 3 * self.__r12) - gamma * self.__r1r2 + 3 * gamma ** 2
        a_3 = (gamma - self.__r_1) ** 3 * (gamma - self.__r_2) ** 3

        return a_0 * (a_1 - a_2 / a_3 * alpha)

    def fug(self, t, p, z):

        v = z * self.r_spc * t / p
        alpha = self.a(t) / (self.b * self.r_spc * t)
        beta = self.b * p / (self.r_spc * t)

        f_v = (np.log(v - self.b * self.__r_1) - np.log(v - self.b * self.__r_2)) / (self.__r_1 - self.__r_2)
        return z - 1 + alpha * f_v - np.log(z - beta)

    def a(self, t):

        return self.a_0 / np.sqrt(t)

    def iterate_t(self, p, v):

        def root_funct(t_rel):

            t = t_rel * self.t_crit
            return p - self.r_spc * t / (v - self.b) + self.a(t) / (v * (v + self.b))

        A = self.r_spc / (v - self.b_vdw)
        B = self.a_vdw / (v * (v + self.b_vdw))
        t_0 = (p + B) / A

        sol = scipy.optimize.root_scalar(root_funct, x0=t_0 / self.t_crit, x1=1)
        return sol.root * self.t_crit

    def init_coefficient(self):

        alpha = self.r_spc * self.t_crit
        beta = np.power(2, 1 / 3) - 1

        self.a_0 = alpha ** 2 * np.sqrt(self.t_crit) / self.p_crit / (9 * beta)
        self.b = beta * alpha / (3 * self.p_crit)
