from main_classes.eos.cubic_eos import CubicEOS
import numpy as np


class PREOS(CubicEOS):

    def a(self, t):

        return self.a_0 * (1 + self.__f_acntr * (1 - np.sqrt(t / self.t_crit))) ** 2

    def da(self, t):

        t_r = np.sqrt(self.t_crit / t)
        return self.a_0 * self.__f_acntr * t_r / t * (self.__f_acntr * (t_r - 1) - 1)

    def dda(self, t):

        t_r = np.sqrt(self.t_crit / t)
        return self.a_0 * self.__f_acntr * (self.__f_acntr + 1) * t_r / (2 * t ** 2)

    def init_coefficient(self):

        self.z_crit = 0.307401
        self.v_crit = self.z_crit * self.r_spc * self.t_crit / self.p_crit

        self.r_1 = -np.sqrt(2) - 1
        self.r_2 = np.sqrt(2) - 1
        self.__f_acntr = 0.37464 + 1.5422 * self.acntr - 0.26992 * self.acntr ** 2

        # v_crit_0 = self.v_crit
        self.a_0 = 0.45724 * (self.r_spc * self.t_crit) ** 2 / self.p_crit
        self.b = 0.07780 * self.r_spc * self.t_crit / self.p_crit