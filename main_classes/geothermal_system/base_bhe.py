from REFPROPConnector import ThermodynamicPoint
from scipy.optimize import minimize
from typing import Union, List
import scipy.constants as cst
from scipy.special import kn
import numpy as np


class ReservoirProperties:

    grad = 0.05  # [°C/m]
    c_rock = 902.67  # [J/(kg K)]
    rho_rock = 2600  # [kg/m^3]
    k_rock = 2.423  # [W/(m K)]
    pe = 1.  # [-]

    __corr = np.array([

        [0.6085, 0.2760, 0.2165],
        [1.3465, 0.4151, 0.7052],
        [0.3777, 0.2792, 0.2195],
        [0.3324, 2.8015, 2.6310],
        [0.2082, 0.0962, 0.5448]

    ])

    @property
    def alpha_rock(self) -> float:
        #   alpha_rocks -> rock thermal diffusivity in [m^2/s]
        return self.k_rock / (self.rho_rock * self.c_rock)

    @property
    def grad_km(self) -> float:
        return self.grad * 1E3

    @grad_km.setter
    def grad_km(self, grad_km: float):
        self.grad = grad_km / 1E3

    def evaluate_rel_resistance(self, times: [float], d: float) -> [float]:

        times_nd = self.get_nd_times(times, d)
        f = self.evaluate_f(times_nd)
        return d / (2 * self.k_rock * f)

    def evaluate_f(self, times_nd: [float]) -> [float]:

        params = self.__corr[:, 1] * np.power(self.pe, self.__corr[:, 0]) + self.__corr[:, 2]

        f = np.empty(times_nd.shape)
        f[:] = params[0]

        for i in range(int(np.floor((len(params) - 1) / 2))):
            f += kn(i, params[1 + i * 2] * np.power(times_nd, params[2 + i * 2]))

        return f

    def get_nd_times(self, times: [float], d: float) -> [float]:
        return np.array(times) * (4 * self.alpha_rock) / np.power(d, 2)

    def get_times(self, nd_times: [float], d: float) -> [float]:
        return np.array(nd_times) / (4 * self.alpha_rock) * np.power(d, 2)


class BHEGeometry:

    depth = 1000  # [m]
    l_horiz = 3000  # [m]
    d_well = 0.15  # [m]
    n_wells = 1

    @property
    def l_tot(self) -> float:
        return self.depth + self.l_horiz


class isobaricIntegral:

    n_it = 10
    n_test = 10

    def __init__(self, limiting_points: [ThermodynamicPoint]):

        self.__init_limiting_points(limiting_points)

        self.__solve_analytically = False
        self.__n_int_num = 20

        self.__delta = 0.0
        self.__dh_max = 0.0
        self.__int_params = []

    # -------------------- READ/WRITE PROPERTIES ---------------------------- #
    # Integral parameter must be recalculated if these properties changes

    def __init_limiting_points(self, limiting_points: [ThermodynamicPoint]):

        self.__limiting_point = limiting_points
        if limiting_points is not None:
            self.__tmp_point = limiting_points[0].duplicate()

        self.__ready = False

    @property
    def limiting_points(self) -> [ThermodynamicPoint]:

        return self.__limiting_point

    @limiting_points.setter
    def limiting_points(self, limiting_points: [ThermodynamicPoint]):

        self.__init_limiting_points(limiting_points)

    @property
    def solve_analytically(self) -> bool:

        return self.__solve_analytically

    @solve_analytically.setter
    def solve_analytically(self, new_solve_analytically: bool):

        if not self.__solve_analytically == new_solve_analytically:
            self.__ready = False

        self.__solve_analytically = new_solve_analytically

    @property
    def n_int_num(self) -> int:

        return self.__n_int_num

    @n_int_num.setter
    def n_int_num(self, new_n_int_num: int):

        if not self.__n_int_num == new_n_int_num:
            self.__ready = False

        self.__n_int_num = new_n_int_num

    # -------------------- READ ONLY PROPERTIES ---------------------------- #

    # If integral parameter have not been recalculated, accessing these
    # properties triggers the calculation

    @property
    def delta(self):

        if not self.__ready:
            self.__evaluate_param()

        return self.__delta

    @property
    def dh_max(self):

        if not self.__ready:
            self.__evaluate_param()

        return self.__dh_max

    @property
    def params(self):

        if not self.__ready:
            self.__evaluate_param()

        return self.__int_params

    # -------------------- BASE METHODS      ---------------------------- #

    # Used to identify the solution and to evaluate the integral

    def __evaluate_param(self):

        if self.__solve_analytically:

            self.__prepare_analytical_integral_param()

        else:

            self.__prepare_numeric_integral_param()

        self.__ready = True

    def identify_solutions(self, UAs: [float]) -> [float]:

        if self.__solve_analytically:

            return self.__identify_analytical_solution(UAs)

        else:

            return self.__identify_numerical_solution(UAs)

    def evaluate_integral(self, dh: [float]) -> [float]:

        if self.__solve_analytically:

            return self.__evaluate_analytical_integral(dh)

        else:

            return self.__evaluate_numerical_integral(dh)

    # -------------------- ANALYTICAL INTEGRATION OF 1/DT ---------------------------- #

    # The following functions perform the analytical Integration of 1/DT dh
    # (with __perform_parallel_bisection being used to retrieve the condition
    # for which the integral equals UA/(m_dot h_max))

    def __prepare_analytical_integral_param(self) -> [float]:

        h_0 = self.__limiting_point[0].get_variable("H")
        dh_max = self.__limiting_point[1].get_variable("H") - h_0

        h_1 = 1 / 3 * dh_max + h_0
        h_2 = 2 / 3 * dh_max + h_0

        t_rocks = self.__limiting_point[1].get_variable("T")
        dt_0 = t_rocks - self.__limiting_point[0].get_variable("T")
        p_0 = self.__limiting_point[0].get_variable("P")

        self.__tmp_point.set_variable("H", h_1)
        self.__tmp_point.set_variable("P", p_0)
        dt_1 = t_rocks - self.__tmp_point.get_variable("T")

        self.__tmp_point.set_variable("H", h_2)
        self.__tmp_point.set_variable("P", p_0)
        dt_2 = t_rocks - self.__tmp_point.get_variable("T")

        # a_1 not used in the program (the knowledge that, for dh = 0 DT=DT_max is used instead)
        # a_1 = (2 * dt_1 - 1 * dt_2 - 11 / 9 * dt_0) * 9 / 2

        a_2 = (4 * dt_2 - 5 * dt_1 + 2 * dt_0) * 9 / 2
        a_3 = (3 * dt_1 - 3 * dt_2 - dt_0) * 9 / 2

        a = -a_3
        b = a - a_2
        c = dt_0

        delta = b ** 2 - 4 * a * c
        self.__delta = delta
        self.__dh_max = dh_max

        if delta >= 0:

            h2 = (-b + np.sqrt(delta)) / (2 * a)
            h3 = (-b - np.sqrt(delta)) / (2 * a)
            A1 = 1 / (a * (1 - h2) * (1 - h3))
            A2 = - 1 / (a * (h2 - h3) * (1 - h2))
            A3 = 1 / (a * (h2 - h3) * (1 - h3))

            self.__int_params = [[1, h2, h3], [A1, A2, A3]]

        else:

            alpha = a + b + c
            beta = np.sqrt(-delta)
            gamma = 2 * (2 * a + b) / beta

            self.__int_params = [[b / c, a / c], [2 * alpha, beta / c, gamma]]

    def __evaluate_analytical_integral(self, dh: [float]) -> [float]:

        delta = self.delta
        int_params = self.__int_params

        if delta > 0:

            h = int_params[0]
            A = int_params[1]

            log_arg = 1

            for i in range(len(h)):
                log_arg *= np.power(h[i] / (h[i] - dh), A[i])

            return np.log(log_arg)

        else:

            dh_sqr = np.power(dh, 2)
            dh_rev_sqr = np.power(1 - dh, 2)
            a = int_params[0]
            par = int_params[1]

            log = np.log((1 + a[0] * dh + a[1] * dh_sqr) / dh_rev_sqr)
            atan = np.arctan((a[0] + 2 * a[1] * dh) * par[1]) - np.arctan(a[0] / par[1])
            return (log + par[2] * atan) / par[0]

    def __identify_analytical_solution(self, UAs: [float]) -> [float]:

        rel_x = np.linspace(start=0, stop=1, num=self.n_test)
        ranges = [{

            "range": [0, 1],
            "values": UAs

        }]

        for i in range(self.n_it):

            new_ranges = list()

            for curr_range in ranges:

                dhs = (curr_range["range"][1] - curr_range["range"][0]) * rel_x + curr_range["range"][0]
                int_res = self.__evaluate_analytical_integral(dhs)

                for j in range(len(int_res) - 1):

                    new_range = np.where(np.logical_and(UAs >= int_res[j], UAs < int_res[j + 1]))

                    if len(new_range[0]) > 0:
                        new_ranges.append({

                            "range": [dhs[j], dhs[j + 1]],
                            "values": UAs[new_range]

                        })

            ranges = new_ranges

        result_list = list()
        x_list = list()
        for curr_range in ranges:
            mean = (curr_range["range"][1] + curr_range["range"][0]) / 2
            result_list.extend(np.ones(len(curr_range["values"])) * mean)
            x_list.extend(curr_range["values"])

        x_list = np.array(x_list)
        result_list = np.array(result_list)

        result_arr = list()
        for UA in UAs:
            result_arr.append(np.mean(result_list[np.where(x_list == UA)]))

        return np.array(result_arr)

    # -------------------- NUMERICAL INTEGRATION OF 1/DT  ---------------------------- #

    # The following functions perform S numerical integration of 1/DT dh
    # (with __perform_parallel_interpolation being used to retrieve the condition
    # for which the integral equals UA/(m_dot h_max))

    def __init_h_rels(self, p_0, h_0, dh_max) -> [float]:

        h_rels = 1 - np.logspace(-3, 0, self.n_int_num)

        # --------------------------- CHECK SATURATION ----------------------------------------------------------
        # If needed, add the two saturation condition (and remove two-phase conditions for single component flow)

        if p_0 < self.__limiting_point[0].RPHandler.PC:

            h_vaps = list()

            for q in [0, 1]:

                try:

                    self.__tmp_point.set_variable("P", p_0)
                    self.__tmp_point.set_variable("Q", q)
                    h_vap = (self.__tmp_point.get_variable("H") - h_0) / dh_max

                    if 0 < h_vap < 1:
                        h_vaps.append(h_vap)

                except:

                    pass

            if len(h_vaps) == 2 and len(self.__limiting_point[0].RPHandler.fluids) == 1:
                # Remove two-phase conditions for single component flow
                h_vaps.sort()
                np.delete(h_rels, np.where(np.logical_and(h_vaps[0] <= h_rels, h_rels <= h_vaps[1])))

            # Append saturation enthalpies
            for h_vap in h_vaps:

                if h_vap not in h_rels:
                    np.insert(h_rels, 0, h_vap)

        h_rels.sort()
        return h_rels

    def __prepare_numeric_integral_param(self) -> [float]:

        h_0 = self.__limiting_point[0].get_variable("H")
        p_0 = self.__limiting_point[0].get_variable("P")

        dh_max = self.__limiting_point[1].get_variable("H") - h_0
        t_rocks = self.__limiting_point[1].get_variable("T")

        h_rels = self.__init_h_rels(p_0, h_0, dh_max)
        dh_rels = h_rels[1:] - h_rels[:-1]

        integral = np.zeros(self.n_int_num)
        dTs = np.zeros(self.n_int_num)
        dTs[0] = t_rocks - self.__limiting_point[0].get_variable("T")

        for i in range(len(dh_rels)):

            self.__tmp_point.set_variable("H", h_rels[i + 1] * dh_max + h_0)
            self.__tmp_point.set_variable("P", p_0)
            dTs[i + 1] = t_rocks - self.__tmp_point.get_variable("T")

            if dTs[i] == dTs[i + 1]:
                dtlm = dTs[i + 1]

            else:
                dtlm = np.log(dTs[i] / dTs[i + 1]) / (dTs[i] - dTs[i + 1])

            integral[i + 1] = integral[i] + dtlm * dh_rels[i]

        h_rels = np.append(h_rels, 1)
        dTs = np.append(dTs, 0)
        integral = np.append(integral, np.inf)

        self.__dh_max = dh_max
        self.__int_params = [h_rels, dTs, integral]

    def __evaluate_numerical_integral(self, dh: [float]) -> [float]:

        pass

    def __identify_numerical_solution(self, UAs: [float]) -> [float]:

        params = self.params

        h_rels = params[0]
        dTs = params[1]
        integral = params[2]
        dh_rels = h_rels[1:] - h_rels[:-1]

        results = np.empty(UAs.shape)
        results[:] = 0.

        for i in range(len(integral) - 1):

            j_values = np.where(np.logical_and(UAs >= integral[i], UAs < integral[i + 1]))

            if len(j_values[0]) > 0:

                ua_rel = (UAs[j_values] - integral[i]) / dh_rels[i]

                if dTs[i + 1] == dTs[i]:

                    rel_hs = ua_rel * dTs[i]

                else:

                    ddt = dTs[i + 1] - dTs[i]
                    rel_hs = (np.exp(ua_rel * ddt) - 1) / ddt * dTs[i]

                results[j_values] = rel_hs * dh_rels[i] + h_rels[i]

        return results


class BaseBHE:

    def __init__(

            self, input_point: ThermodynamicPoint,
            geometry: BHEGeometry = None, reservoir_properties: ReservoirProperties = None

    ):

        if geometry is None:
            geometry = BHEGeometry()

        if reservoir_properties is None:
            reservoir_properties = ReservoirProperties()

        self.geom = geometry
        self.res_prop = reservoir_properties

        self.__user_unit_system = input_point.RPHandler.unit_system
        si_intput = input_point  # .get_alternative_unit_system("MASS BASE SI")

        self.__ideal_points = [si_intput, si_intput.duplicate(), si_intput.duplicate(), si_intput.duplicate()]
        self.__real_points = [si_intput.duplicate(), si_intput.duplicate()]
        self.__tmp_point = si_intput.duplicate()

        self.integrator = isobaricIntegral([si_intput, si_intput])

    @property
    def ideal_exergy_efficiency(self) -> float:

        dh = self.__ideal_points[-1].get_variable("H") - self.__ideal_points[0].get_variable("H")
        ds = self.__ideal_points[-1].get_variable("S") - self.__ideal_points[0].get_variable("S")

        t_0 =  self.__ideal_points[0].get_variable("T")
        t_rocks = t_0 + self.res_prop.grad * self.geom.depth
        dex_rocks = dh * (1 - t_0 / t_rocks)
        return (dh - t_0 * ds) / dex_rocks

    @property
    def input_point(self) -> ThermodynamicPoint:

        return self.__ideal_points[0].get_alternative_unit_system(self.__user_unit_system)

    @input_point.setter
    def input_point(self, input_point):

        self.__ideal_points[0] = input_point.get_alternative_unit_system("MASS BASE SI")

    @property
    def ideal_points(self) -> list[ThermodynamicPoint]:

        return_list = list()
        for ideal_point in self.__ideal_points:
            return_list.append(ideal_point.get_alternative_unit_system(self.__user_unit_system))

        return return_list

    @property
    def real_points(self) -> list[ThermodynamicPoint]:

        return_list = [

            self.__ideal_points[0].get_alternative_unit_system(self.__user_unit_system),
            self.__ideal_points[1].get_alternative_unit_system(self.__user_unit_system),

        ]

        for real_point in self.__real_points:

            return_list.append(real_point.get_alternative_unit_system(self.__user_unit_system))

        return return_list

    def set_HX_condition(self):

        self.__evaluate_ideal_thermo()

    def __evaluate_ideal_thermo(self):

        dh = cst.g * self.geom.depth
        dt_max = self.res_prop.grad * self.geom.depth

        self.__ideal_points[1].set_variable("H", self.__ideal_points[0].get_variable("H") + dh)
        self.__ideal_points[1].set_variable("S", self.__ideal_points[0].get_variable("S"))

        self.__ideal_points[2].set_variable("T", self.__ideal_points[0].get_variable("T") + dt_max)
        self.__ideal_points[2].set_variable("P", self.__ideal_points[1].get_variable("P"))

        self.__ideal_points[3].set_variable("H", self.__ideal_points[2].get_variable("H") - dh)
        self.__ideal_points[3].set_variable("S", self.__ideal_points[2].get_variable("S"))

        self.integrator.limiting_points = [self.__ideal_points[1], self.__ideal_points[2]]

    def evaluate_HXG(self, times: Union[float, List[float], np.ndarray], m_dot: Union[float, List[float], np.ndarray]) -> dict:

        times, m_dot = np.meshgrid(np.array(times), np.array(m_dot))

        dh_max = self.integrator.dh_max
        UdAs = 1 / self.res_prop.evaluate_rel_resistance(times=times, d=self.geom.d_well) * np.pi * self.geom.d_well
        UAs = UdAs / dh_max * (self.geom.l_horiz * self.geom.n_wells / m_dot)
        dh_percs = self.integrator.identify_solutions(UAs)

        dhs = dh_percs * dh_max
        original_shape = dhs.shape
        dhs_flat = dhs.flatten()

        ds_flat = np.empty(dhs_flat.shape)
        dex_flat = np.empty(dhs_flat.shape)

        p_down_flat = np.empty(dhs_flat.shape)
        h_down_flat = np.empty(dhs_flat.shape)
        drho_flat = np.empty(dhs_flat.shape)
        t_down_flat = np.empty(dhs_flat.shape)

        ds_flat[:] = np.nan
        dex_flat[:] = np.nan
        p_down_flat[:] = np.nan
        h_down_flat[:] = np.nan
        drho_flat[:] = np.nan
        t_down_flat[:] = np.nan

        for i in range(len(dhs_flat)):

            self.__real_points[0].set_variable("H", self.__ideal_points[1].get_variable("H") + dhs_flat[i])
            self.__real_points[0].set_variable("P", self.__ideal_points[1].get_variable("P"))

            ds = self.__real_points[0].get_variable("S") - self.__ideal_points[1].get_variable("S")

            if ds > 0:
                ds_flat[i] = self.__real_points[0].get_variable("S") - self.__ideal_points[1].get_variable("S")
                dex_flat[i] = dhs_flat[i] - self.__ideal_points[0].get_variable("T") * ds_flat[i]
                drho_flat[i] = self.__ideal_points[1].get_variable("rho") - self.__real_points[0].get_variable("rho")

                t_down_flat[i] = self.__real_points[0].get_variable("T")
                p_down_flat[i] = self.__real_points[0].get_variable("P")
                h_down_flat[i] = self.__real_points[0].get_variable("H")

        return {

            "time": times,
            "m_dot": m_dot,
            "dh_perc": dh_percs,
            "dh": dhs,
            "ds": ds_flat.reshape(original_shape),
            "dex": dex_flat.reshape(original_shape),
            "drho": drho_flat.reshape(original_shape),
            "T_down": t_down_flat.reshape(original_shape),
            "P_down": p_down_flat.reshape(original_shape),
            "H_down": h_down_flat.reshape(original_shape),

        }

    def evaluate_surface_condition(self, HXG_result: dict, indices: Union[List[int], np.ndarray] = None) -> dict:

        p_downs = HXG_result["P_down"]
        h_downs = HXG_result["H_down"]
        original_shape = p_downs.shape

        p_flat = p_downs.flatten()
        h_flat = h_downs.flatten()

        p_out_flat = np.empty(p_flat.shape)
        t_out_flat = np.empty(p_flat.shape)
        dp_flat = np.empty(p_flat.shape)
        dt_flat = np.empty(p_flat.shape)

        p_out_flat[:] = np.nan
        t_out_flat[:] = np.nan
        dp_flat[:] = np.nan
        dt_flat[:] = np.nan

        dh_vert = cst.g * self.geom.depth

        if indices is None:
            indices = range(len(p_flat))

        for i in indices:

            if not np.isnan(p_flat[i]):

                self.__real_points[0].set_variable("P", p_flat[i])
                self.__real_points[0].set_variable("H", h_flat[i])

                self.__real_points[1].set_variable("H", self.__real_points[0].get_variable("H") - dh_vert)
                self.__real_points[1].set_variable("S", self.__real_points[0].get_variable("S"))

                p_out = self.__real_points[1].get_variable("P")
                t_out = self.__real_points[1].get_variable("T")
                rho_out = self.__real_points[1].get_variable("rho")

                if t_out > 0 and p_out > 0 and rho_out > 0:

                    p_out_flat[i] = self.__real_points[1].get_variable("P")
                    t_out_flat[i] = self.__real_points[1].get_variable("T")
                    dp_flat[i] = self.__real_points[1].get_variable("T") - self.__ideal_points[0].get_variable("T")
                    dp_flat[i] = self.__real_points[1].get_variable("T") - self.__ideal_points[0].get_variable("T")

        HXG_result.update({

            "p_out": p_out_flat.reshape(original_shape),
            "T_out": t_out_flat.reshape(original_shape),
            "dt_overall": dt_flat.reshape(original_shape),
            "dp_overall": dp_flat.reshape(original_shape),

        })

        return HXG_result

    def evaluate_direct_expansion(self, surface_result: dict, indices: Union[List[int], np.ndarray] = None) -> dict:

        p_surfs = surface_result["p_out"]
        h_surfs = surface_result["T_out"]
        original_shape = p_surfs.shape

        p_flat = p_surfs.flatten()
        h_flat = h_surfs.flatten()

        dh_turb_flat = np.empty(p_flat.shape)
        dh_cool_flat = np.empty(p_flat.shape)

        dh_turb_flat[:] = np.nan
        dh_cool_flat[:] = np.nan

        if indices is None:
            indices = range(len(p_flat))

        for i in indices:

            if not np.isnan(p_flat[i]):

                self.__real_points[1].set_variable("P", p_flat[i])
                self.__real_points[1].set_variable("H", h_flat[i])

                self.__tmp_point.set_variable("P", self.__ideal_points[0].get_variable("P"))
                self.__tmp_point.set_variable("S", self.__real_points[1].get_variable("S"))

                dh_turb_flat[i] = self.__real_points[1].get_variable("H") - self.__tmp_point.get_variable("H")
                dh_cool_flat[i] = self.__tmp_point.get_variable("H") - self.__ideal_points[0].get_variable("H")

        surface_result.update({

            "dh_turb": dh_turb_flat.reshape(original_shape),
            "dh_cool": dh_cool_flat.reshape(original_shape)}

        )

        return surface_result


class economicEvaluator:

    Le = 20
    i_rate = 0.04
    om_ratio = 0.05
    hy = 8000
    c_el = 0.075

    def __init__(self, thermo_evaluator: BaseBHE):
        self.thermo = thermo_evaluator

    @property
    def alpha(self):
        return (1 - (1 + self.i_rate) ** (-self.Le)) / self.i_rate

    @property
    def beta(self):
        return (1 + self.alpha * self.om_ratio) / (self.alpha * self.hy)

    def evaluate_LCOx(self, xs: [str], len_hor: float, m_dot: float) -> [float]:

        # self.thermo.set_HX_condition()
        self.thermo.geom.l_horiz = len_hor

        l_overall = self.thermo.geom.l_tot
        d_well = self.thermo.geom.d_well
        c_well = 1.15 * 1.05 * 2.86 * (0.105 * l_overall ** 2 + 1776 * l_overall * d_well + 2.735E5)
        self.thermo.evaluate_HXG([10 * 3.154e+7], m_dot=m_dot)

        result_list = list()
        for x in xs:

            if hasattr(self.thermo, x):

                result_list.append((c_well * self.beta) / getattr(self.thermo, x))

            else:

                result_list.append(np.nan)

        return result_list

    def optimize_LOCx(self, xs: str) -> [float]:

        x0 = 300
        LCOx_0 = self.evaluate_LCOx([xs], len_hor=x0, m_dot=1)[0]

        def func(x):

            self.thermo.geom.l_horiz = x[0]*x0

            l_overall = self.thermo.geom.l_tot
            d_well = self.thermo.geom.d_well
            c_well = 1.15 * 1.05 * 2.86 * (0.105 * l_overall ** 2 + 1776 * l_overall * d_well + 2.735E5)
            self.thermo.evaluate_HXG([10 * 3.154e+7], m_dot=1)

            LCOx = (c_well * self.beta) / getattr(self.thermo, xs)
            print("{} - {}".format(x, LCOx))

            return LCOx / LCOx_0

        print(minimize(func, np.array([1])))
