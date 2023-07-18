from main_classes.eos.cubic_eos import CubicEOS


class VdWEOS(CubicEOS):

    def a(self, t):
        return self.a_0

    def da(self, t):
        return 0.

    def dda(self, t):
        return 0.

    def init_coefficient(self):

        self.a_0 = self.a_vdw
        self.b = self.b_vdw

        self.z_crit = 3 / 8
        self.r_1 = 0
        self.r_2 = 0

    def iterate_t(self, p, v):

        A = self.r_spc / (v - self.b)
        B = self.a_0 / (v * (v + self.b))
        return (p + B) / A
