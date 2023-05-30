from main_classes.cubic_eos import CubicEOS, get_real_res
import numpy as np

class RKEOS(CubicEOS):

    def a(self, t):

        return self.a_0 / np.sqrt(t)

    def da(self, t):

        return -self.a_0 / (2 * np.power(t, 3 / 2))

    def dda(self, t):

        return 3 * self.a_0 / (4 * np.power(t, 5 / 2))

    def init_coefficient(self):

        alpha = self.r_spc * self.t_crit
        beta = np.power(2, 1 / 3) - 1

        self.a_0 = alpha ** 2 * np.sqrt(self.t_crit) / self.p_crit / (9 * beta)
        self.b = beta * alpha / (3 * self.p_crit)

        self.z_crit = 1 / 3
        self.r_1 = 0
        self.r_2 = -1

    def iterate_t(self, p, v):

        A = self.a_0 / (v * (v + self.b))
        B = self.r_spc / (v - self.b)

        res = np.roots([B ** 2, -2 * B * p, p ** 2, - A ** 2])
        t_l, t_v = get_real_res(res)

        return t_l
