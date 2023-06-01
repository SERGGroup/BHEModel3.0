from main_classes.cubic_eos import CubicEOS, get_real_res
import numpy as np


class SRKEOS(CubicEOS):

    def a(self, t):

        return self.a_0 * (1 + self.__f_acntr * (1 - np.sqrt(t / self.t_crit)) ** 2)

    def da(self, t):

        t_r = np.sqrt(self.t_crit / t)
        return self.a_0 * self.__f_acntr / t * self.__f_acntr * t_r * (self.__f_acntr * (t_r - 1) - 1)

    def dda(self, t):

        t_r = np.sqrt(self.t_crit / t)
        return self.a_0 * self.__f_acntr * (self.__f_acntr + 1) / 2 * t_r / t ** 2

    def init_coefficient(self):

        self.z_crit = 1 / 3
        self.v_crit = self.z_crit * self.r_spc * self.t_crit / self.p_crit

        self.r_1 = 0
        self.r_2 = -1
        self.__f_acntr = 0.48 + 1.57 * self.acntr - 0.176 * self.acntr ** 2

        # v_crit_0 = self.v_crit
        self.a_0 = (np.power(2, 1/3) - 1) / 9 * (self.r_spc * self.t_crit) ** 2 / self.p_crit
        self.b = (np.power(2, 1/3) - 1) / 3 * self.r_spc * self.t_crit / self.p_crit

        self.a_0 = 0.42747 * (self.r_spc * self.t_crit) ** 2 / self.p_crit
        self.b = 0.08664 * self.r_spc * self.t_crit / self.p_crit

