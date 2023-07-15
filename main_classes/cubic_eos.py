from main_classes.support.p_sat_calculation import p_sat_rel, t_sat_rel
from scipy.optimize import root_scalar
from abc import ABC, abstractmethod
from copy import deepcopy
from ctypes import Union
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

    """

        Base class representing a Cubic Equation of State
        It cannot be implemented as it is, but it provides the subclasses with most of the calculation procedures.

        Procedures implemented in this class:
            - EoS solution (root finding of the cubic equation, Saturation Condition Identification, physical root selection)

            - State Class initialization

        Methods that must be implemented in subclasses:
            - init_coefficient(): method that must set the correct values of a_0, v_crit and b starting from t_crit and p_crit

            - method connected with the attraction term (a(t), da(t), dda(t)) in which the function a(t) should be implemented

        NOTE: cp_ideal can be provided either as a float or as a list. In the latter case, the list should contain
        the parameters for a polinomial regression (eg: cp_ideal=[cp0, cp1, cp2] implies that the class will
        calculate the ideal gas heat capacity as: cp = cp0 + cp1 * T + cp2 * T^2)

        NOTE: default values in the initialization are for C02

        :param p_crit: critical pressure [Pa]
        :param t_crit: critical temperature [K]
        :param cp_ideal: ideal gas specific heat capacity [J/kg K] or list containing params for a regression
        :param m_molar: molar mass [kg/mol]
        :param acntr: acentricity factor (omega) - if required by the EoS [-]

    """

    def __init__(

            self, p_crit=7380000., t_crit=304.1,
            cp_ideal=845.85, m_molar=0.04401,
            acntr=0.239

    ):

        self.r_spc = scipy.constants.R / m_molar

        self.p_crit = p_crit
        self.t_crit = t_crit
        self.__cp_ideal = cp_ideal
        self.m_molar = m_molar
        self.acntr = acntr

        if type(self.__cp_ideal) == list:

            # If a list of coefficient is provided, it has to be
            # reversed for numerical reason (see the self.cp_ideal(t) method)
            self.__cp_ideal.reverse()

        self.__init_coefficients()

    def __init_coefficients(self):

        """

            Initialize the coefficient to be used for the calculation. It also call init_coefficient() that must be
            implemented in subclasses and should initialize self.a_0, self.b, self.z_crit, self.r_1 and self.r_2
            according to the selected EoS.

        """

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

        """

            Return a state class initialized according to the provided state variable, the user should provide two of
            them (i.e. fluid.get_state(t=300, p=1E5) or fluid.get_state(v=0.05, p=1E5))

            :param p: pressure in [Pa]
            :param t: temperature in [K]
            :param v: specific volume in [m^3/kg]

            :return: state class

        """

        if (t is not None) and (v is not None):

            return FluidState(self, t, v)

        elif (t is not None) and (p is not None):

            return FluidState(self, t, self.v(t, p))

        elif (p is not None) and (v is not None):

            return FluidState(self, self.t(v, p), v)

        return FluidState(self)

    def get_sat_state(self, p=None, t=None, which="both"):

        """

            Return a state class referring to the saturation state identified by the provided state variable,
            the user should provide one of them (i.e. fluid.get_sat_state(t=300) or fluid.get_sat_state(p=1E5)). the
            "which" parameter can be used to force the method to return the state corresponding to the saturated liquid,
            the saturated vapour or both of them (by default it returns both)

            :param p: pressure in [Pa]
            :param t: temperature in [K]
            :param which: if "both" the programs return both the saturated liquid and vapour state (default), other possibilities are "liq", "Liquid", "l" for returning the liquid only or "vap", "vapour", "v" for returning the vapour. the input IS NOT case sensitive

            :return: state class

        """
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

        """

            Return the predicted pressure given the provided temperature and specific volume

            :param t: temperature in [K]
            :param v: specific volume in [m^3/kg]

            :return: p: pressure in [Pa]

        """

        if t < self.t_crit:

            """
                
                If the provided temperature is less than the critical one, check if the specific volume provided lies 
                between the saturation condition. If its so it returns the saturation pressure, otherwise return the 
                result of the EoS. 
                
                As the p_sat calculation takes time, if an approximation is available the code uses it in order to 
                check the likelihood for the specific point to be in the saturation condition. If so, it iterates the 
                EoS to indentify the exact saturation pressure otherwise it simply returns the EoS prediction
                            
            """

            # Check for approximate prediction
            p_sat, z_l, z_v = self.p_sat(t, get_approximate=True)

            if p_sat is None:

                # Approximation are not allowed for the selected EoS (None has been returned),
                # the code iterates p_sat in any condition

                p_sat, z_l, z_v = self.p_sat(t, get_approximate=False)
                v_liq = z_l * self.r_spc * t / p_sat
                v_vap = z_v * self.r_spc * t / p_sat

                if v_liq < v < v_vap:

                    # If the provided specific volume lies between the saturation
                    # condition, the code returns p_sat

                    return p_sat

            else:

                # The approximation exists and has
                # returned the predicted saturation condition

                v_liq = z_l * self.r_spc * t / p_sat
                v_vap = z_v * self.r_spc * t / p_sat

                if v_liq * 0.98 < v < v_vap * 1.02:

                    # The specific volume is within or near the
                    # saturation range, p_sat is iterated

                    p_sat, z_l, z_v = self.p_sat(t, get_approximate=False)
                    v_liq = z_l * self.r_spc * t / p_sat
                    v_vap = z_v * self.r_spc * t / p_sat

                    if v_liq < v < v_vap:

                        # If the provided specific volume lies between the saturation
                        # condition, the code returns p_sat

                        return p_sat

        # If the execution has reached this point, it means that
        # the specific volume is outside the saturation range, hence
        # the direct result of the EoS is returned

        z_l, z_v = self.z(t=t, v=v)
        return z_l * self.r_spc * t / v

    def t(self, v, p):

        """

            Return the predicted temperature given the provided pressure and specific volume

            :param p: pressure in [Pa]
            :param v: specific volume in [m^3/kg]

            :return: t: temperature in [K]

        """

        if p < self.p_crit:

            """

                If the provided pressure is less than the critical one, check if the specific volume provided lies 
                between the saturation condition. If its so it returns the saturation temperature, otherwise return the 
                result of the EoS. 

                As the t_sat calculation takes time, if an approximation is available the code uses it in order to 
                check the likelihood for the specific point to be in the saturation condition. If so, it iterates the 
                EoS to indentify the exact saturation pressure otherwise it simply returns the EoS prediction

            """

            # Check for approximate prediction
            t_sat, z_l, z_v = self.t_sat(p, get_approximate=True)

            if t_sat is None:

                # Approximation are not allowed for the selected EoS (None has been returned),
                # the code iterates t_sat in any condition

                t_sat, z_l, z_v = self.t_sat(p, get_approximate=False)
                v_liq = z_l * self.r_spc * t_sat / p
                v_vap = z_v * self.r_spc * t_sat / p

                if v_liq < v < v_vap:

                    # If the provided specific volume lies between the saturation
                    # condition, the code returns t_sat

                    return t_sat

            else:

                # The approximation exists and has
                # returned the predicted saturation condition

                v_liq = z_l * self.r_spc * t_sat / p
                v_vap = z_v * self.r_spc * t_sat / p

                if v_liq * 0.98 < v < v_vap * 1.02:

                    # The specific volume is within or near the
                    # saturation range, p_sat is iterated

                    t_sat, z_l, z_v = self.t_sat(p, get_approximate=False)
                    v_liq = z_l * self.r_spc * t_sat / p
                    v_vap = z_v * self.r_spc * t_sat / p

                    if v_liq < v < v_vap:

                        # If the provided specific volume lies between the saturation
                        # condition, the code returns t_sat

                        return t_sat

        # If the execution has reached this point, it means that
        # the specific volume is outside the saturation range, hence
        # the direct result of the EoS is returned

        z_l, z_v = self.z(v=v, p=p)
        return p * v / (self.r_spc * z_l)

    def v(self, t, p):

        """

            Return the predicted specific volume given the provided pressure and temperature

            :param t: temperature in [K]
            :param p: pressure in [Pa]

            :return: v: specific volume in [m^3/kg]

        """

        # The method identify physical roots of the EoS. If only one root exist z_l == z_v, hence the corresponding
        # specific volume can be directly returned. Otherwise, the method evaluates the fugacity of both the conditions
        # and returns the condition with the smallest one.

        # EoS root identification
        z_l, z_v = self.z(t=t, p=p)

        if z_l == z_v:

            # Only one physical root
            z = z_l

        else:

            # Two roots, evaluate fugacity
            f_l = self.fug(t, p, z_l)
            f_v = self.fug(t, p, z_v)

            # Return the condition with lower fugacity
            if f_v > f_l:

                z = z_l

            else:

                z = z_v

        # Evaluate the specific volume starting
        # from the compressibility factor

        return z * self.r_spc * t / p

    def z(self, p=None, t=None, v=None):

        """

            Evaluate the compressibility factor from the EoS according to the provided state variable, the user should
            provide two of them (i.e. fluid.z(t=300, p=1E5) or fluid.z(v=0.05, p=1E5)). If multiple solution are
            possible it returns different values for z_l and z_v, otherwise it returns the same value for both.

            :param p: pressure in [Pa]
            :param t: temperature in [K]
            :param v: specific volume in [m^3/kg]

            :return: z_l, z_v (if different, compressibility factors for liquid and vapour phase)

        """

        if (t is not None) and (v is not None):

            # If T and v are provided it directly evaluate p from the EoS
            # (only 1 solution exists)

            p = self.r_spc * t / (v - self.b) - self.a(t) / ((v - self.b * self.r_1) * (v - self.b * self.r_2))
            z_l = p * v / (self.r_spc * t)
            z_v = z_l

        elif (t is not None) and (p is not None):

            # If T and p are provided it evaluates the
            # root of the resulting cubic equation for v
            alpha = self.a(t) / (self.b * self.r_spc * t)
            beta = self.b * p / (self.r_spc * t)

            A_1 = beta * (self.r1r2 + 1) + 1
            A_2 = beta * (beta * self.r12 + alpha + self.r1r2 * (beta + 1))
            A_3 = beta ** 2 * (self.r12 * (beta + 1) + alpha)

            # np.roots solve the cubic equation and returns
            # both the real and the complex solutions
            res = np.roots([1, -A_1, A_2, -A_3])

            # Get_real_res analyze the roots and returns the physical ones
            z_l, z_v = get_real_res(res)

        elif (p is not None) and (v is not None):

            # If v and p are provided it iterates the EoS to identify the resulting temperature,
            # For simpler EoS, iterate_t(p, v) can be overwritten in subclasses to allow for explicit solutions
            # (only 1 solution exists)

            t = self.iterate_t(p, v)
            z_l = p * v / (self.r_spc * t)
            z_v = z_l

        else:

            # Not enough state variables have been provided
            # np.nan returned
            z_l = np.nan
            z_v = np.nan

        return z_l, z_v

    def p_sat(self, t, get_approximate=False):

        """

            Evaluate the saturation pressure in [Pa]. If "get_approximate" is True, the method tries to evaluate
            the saturation pressure using the provided correlation (None is returned if no interpolation
            is possible for the selected subclass)

            :param t: temperature in [K]
            :param get_approximate: if true tries to approximate the saturation condition
            :return: p_sat, z_l, z_v ---- (p_sat in [Pa])

        """

        return self.__get_sat(t, get_approximate=get_approximate, get_p=True)

    def t_sat(self, p, get_approximate=False):
        """

            Evaluate the saturation temperature in [K]. If "get_approximate" is True, the method tries to evaluate
            the saturation temperature using the provided correlation (None is returned if no interpolation
            is possible for the selected subclass)

            :param p: pressure in [Pa]
            :param get_approximate: if true tries to approximate the saturation condition
            :return: t_sat, z_l, z_v ---- (t_sat in [K])

        """

        return self.__get_sat(p, get_approximate=get_approximate, get_p=False)

    def __get_sat(self, y, get_approximate=False, get_p=True):

        """

            Evaluate the saturation condition given the provided input. If "get_approximate" is True, the method tries
            to evaluate the saturation condition using the provided correlation (None is returned if no interpolation
            is possible for the selected subclass). If "get_p" is true, the method consider y to be the saturation
            temperature and evaluates p_sat (the opposite otherwise).

            :param y: pressure in [Pa] or temperature in [K] (depending on "get_p")
            :param get_approximate: if true tries to approximate the saturation condition
            :param get_p: if True y is considered to be T_sat and the function return p_sat (otherwise is the opposite)

            :return: t_sat or p_sat, z_l, z_v ---- (t_sat in [K], p_sat in [Pa])

        """
        if get_approximate:

            # Run the method that approximate p_sat or t_sat
            # --> (currently disabled) <--
            # return self.__approx_sat(y, get_p=get_p)
            return None, None, None

        else:

            # Run the standard iteration
            return self.__iterate_sat(y, get_p=get_p)

    def __iterate_sat(self, y, get_p):

        """

            Standard iteration procedure for the saturation condition using the bisection method

            :param y: pressure in [Pa] or temperature in [K] (depending on "get_p")
            :param get_p: if True y is considered to be T_sat and the function return p_sat (otherwise is the opposite)

            :return: t_sat or p_sat, z_l, z_v ---- (t_sat in [K], p_sat in [Pa])

        """

        if get_p:

            # "get_p" is true:
            #
            #   y is the saturation temperature
            #   x is the unknown saturation pressure
            #
            # The critical values are initialized accordingly

            x_crit = self.p_crit
            y_crit = self.t_crit

        else:

            # "get_p" is false:
            #
            #   y is the saturation pressure
            #   x is the unknown saturation temperature
            #
            # The critical values are initialized accordingly

            x_crit = self.t_crit
            y_crit = self.p_crit

        # Check if the provided saturation condition is allowable
        # (if it is above the critical point it cannot be used)
        if y < y_crit:

            # If it can be used the bisection method is initialized. The reduced form of the unknown saturation value
            # must be confined between 0 and 1. The reduced form being the ratio between the variable and its
            # critical value (e.g. p_rel = p / p_crit)
            x_rels = [0, 1]
            x_curr, z_l, z_v = np.zeros(3)

            n = 0
            # The iteration begins (25 steps are considered enough to achieve an acceptable tolerance)
            # TODO: implement a better tolerance
            while n < 25:

                # get the guess value
                x_rel = np.mean(x_rels)
                x_curr = x_rel * x_crit

                # evaluate the error
                # (different depending on "get_p")
                if get_p:

                    err, z_l, z_v = self.__error_sat(t=y, p=x_curr)

                else:

                    # reverse the error while looking for t_sat in order to
                    # maintain the same bisection update method:
                    #
                    # In any case, the error returned by "self.__error_sat" is positive if
                    # the current condition is on the VAPOUR side of the diagram (which case
                    # the pressure is too low or the temperature is too high)
                    #
                    # In order to make the error positive when the guess value is
                    # too low, the result given by "self.__error_sat" is inverted if we
                    # are looking for the saturation temperature

                    err, z_l, z_v = self.__error_sat(t=x_curr, p=y)
                    err = -err

                if err < 0:

                    # if error < 0, the guess value is too high hence it is set as new high limit of the range
                    # (see __error_sat for clarification)

                    x_rels[1] = x_rel

                else:

                    # if error > 0, the guess value is too low hence it is set as new low limit of the range
                    # (see __error_sat for clarification)

                    x_rels[0] = x_rel

                n += 1

            # once the iteration is completed, return the results

            return x_curr, z_l, z_v

        elif y == y_crit:

            # If the provided saturation condition is equal to the
            # critical point, use the critical condition as a result

            return x_crit, self.z_crit, self.z_crit

        return np.nan, np.nan, np.nan

    def __approx_sat(self, y, get_p):

        """

            Approximation of the saturation condition using a polynomial regression

            :param y: pressure in [Pa] or temperature in [K] (depending on "get_p")
            :param get_p: if True y is considered to be T_sat and the function return p_sat (otherwise is the opposite)

            :return: t_sat or p_sat, z_l, z_v ---- (t_sat in [K], p_sat in [Pa])

        """

        cs = self.sat_coefficients

        if self.sat_coefficients is not None:

            # The "self.sat_coefficients" has been implemented in the subclass hence
            # the approximation can start

            if get_p:

                # Call the calculation function (imported from the support module)
                x = p_sat_rel(y / self.t_crit, cs) * self.p_crit

                # Evaluate the resulting compressibility factors
                z_l, z_v = self.z(t=y, p=x)

            else:

                # Call the calculation function (imported from the support module)
                x = t_sat_rel(y / self.p_crit, cs) * self.t_crit

                # Evaluate the resulting compressibility factors
                z_l, z_v = self.z(t=x, p=y)

            return x, z_l, z_v

        else:

            # The "self.sat_coefficients" has NOT been implemented in the subclass hence
            # the approximation returns None and ask the calling method to proceed with
            # the standard iteration

            return None, None, None

    def __error_sat(self, t, p):

        """

            Return the error committed in considering t and p to be the saturation condition.
            (NOTE: By convention, the error is POSITIVE if the current condition is on the VAPOUR SIDE of the diagram)

            :param t: temperature in [K]
            :param p: pressure in [Pa]

            :return: error, z_l, z_v

        """

        # Identify the root of the EoS
        z_l, z_v = self.z(t=t, p=p)

        if not z_l == z_v:

            # If the EoS has two real roots, the error committed is the difference between the two fugacities
            # (NOTE: The error is POSITIVE if the current condition is on the VAPOUR SIDE of the diagram)

            f_l = self.fug(t, p, z_l)
            f_v = self.fug(t, p, z_v)

            return f_l - f_v, z_l, z_v

        else:

            # If the EoS has one real root the method use the xi parameter to understand it the current condition is on
            # the vapour or on the liquid side of the diagram.
            # (see the documentation for additional information)
            #
            # (NOTE: The error is POSITIVE if the current condition is on the VAPOUR SIDE of the diagram)

            v = z_l * self.r_spc * t / p
            error = -(self.a(t) / (v - self.b) - self.a(self.t_crit) / (self.v_crit - self.b))

            return error, z_l, z_v

    def fug(self, t, p, z):

        """

            Return the fugacity of the selected condition

            :param t: temperature in [K]
            :param p: pressure in [Pa]
            :param z: compressibility factor [-]

            :return: fugacity [-]

        """

        # See the documentation for the analytic formula

        v = z * self.r_spc * t / p
        alpha = self.a(t) / (self.b * self.r_spc * t)
        beta = self.b * p / (self.r_spc * t)

        f_v = (np.log(v - self.b * self.r_1) - np.log(v - self.b * self.r_2)) / (self.r_1 - self.r_2)
        return z - 1 + alpha * f_v - np.log(z - beta)

    def cp_ideal(self, t) -> float:

        """

            Return the ideal gas heat capacity for the given temperature

            :param t: temperature in [K]
            :return: ideal gas heat capacity [J/(kg * K)]

        """

        if type(self.__cp_ideal) == list:

            # If "self.__cp_ideal" is a list, a polinomial regression has to be evaluated.
            # The polinomial is evaluated in the form:
            #
            #       y = ((a_0 * T + a_1) * T + a_2) * T + a_3 ...
            #
            # For this reason, the order of the coefficients has to be reversed if compared
            # to the standard polinomial representation

            cp_ideal = 0.
            # noinspection PyTypeChecker
            for cp_coefficient in self.__cp_ideal:

                cp_ideal = cp_ideal * t + cp_coefficient

        else:

            # If "self.__cp_ideal" is a float, no regression is needed
            # and the value is directly returned
            cp_ideal = self.__cp_ideal

        return cp_ideal

    @property
    def cp_coefficients(self):

        """

            Ideal gas heat capacity coefficients

        """

        if type(self.__cp_ideal) == list:

            self.__cp_ideal.reverse()
            param = deepcopy(self.__cp_ideal)
            self.__cp_ideal.reverse()
            return param

        return self.__cp_ideal

    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #
    #                                 FUNCTION THAT MUST BE UPDATED DEPENDING ON THE EOS
    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #

    @abstractmethod
    def a(self, t):
        """

            Return the attraction term of the EoS
            TO BE IMPLEMENTED IN SUBCLASS

            :param t: temperature in [K]
            :return: attraction term of the EoS [m^5/(kg * s^2)]

        """
        return self.a_0 / np.sqrt(t)

    @abstractmethod
    def da(self, t):
        """

            Return the derivative of the attraction term of the EoS
            TO BE IMPLEMENTED IN SUBCLASS

            :param t: temperature in [K]
            :return: derivative of the attraction term of the EoS [m^5/(kg * s^2 * K)]

        """
        return -self.a_0 / (2 * np.power(t, 3/2))

    @abstractmethod
    def dda(self, t):

        """

            Return the second derivative of the attraction term of the EoS
            TO BE IMPLEMENTED IN SUBCLASS

            :param t: temperature in [K]
            :return: second derivative of the attraction term of the EoS [m^5/(kg * s^2 * K^2)]

        """

        return 3 * self.a_0 / (4 * np.power(t, 5/2))

    @abstractmethod
    def init_coefficient(self):

        """

            Function that must be used by the subclass to initialize
            coefficients to be used during the calculations
            TO BE IMPLEMENTED IN SUBCLASS

        """
        alpha = self.r_spc * self.t_crit
        beta = np.power(2, 1 / 3) - 1

        self.a_0 = alpha ** 2 * np.sqrt(self.t_crit) / self.p_crit / (9 * beta)
        self.b = beta * alpha / (3 * self.p_crit)

        self.z_crit = 1/3
        self.r_1 = 0
        self.r_2 = -1

    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #
    #                                 FUNCTION THAT CAN BE UPDATED DEPENDING ON THE EOS
    # <--------------------------------------------------------------------------------------------------------------> #
    # <--------------------------------------------------------------------------------------------------------------> #

    def iterate_t(self, p, v):

        """

            Function that iterate the EoS to find the temperature given
            the pressure and the specific volume. Can be overwritten in
            subclasses if analytic solution are possible

            :param p: pressure in [Pa]
            :param v: specific volume in [m^3/kg]

            :return: t: temperature in [K]

        """

        def root_funct(t_rel):

            t = t_rel * self.t_crit
            return p - self.r_spc * t / (v - self.b) + self.a(t) / (v * (v + self.b))

        # Initialize the calculation with the Van der Walls solution
        A = self.r_spc / (v - self.b_vdw)
        B = self.a_vdw / (v * (v + self.b_vdw))
        t_0 = (p + B) / A

        # Find the root of the function
        sol = root_scalar(root_funct, x0=t_0 / self.t_crit, x1=1)
        return sol.root * self.t_crit

    @property
    def sat_coefficients(self):

        return None


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
        self.__cp_ideal = None

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

    def update_state(self, t=None, v=None, p=None):

        if (t is not None) and (v is not None):

            self.__t = t
            self.__v = v

        elif t is None:

            self.__t = self.fluid_solver.t(v, p)
            self.__v = v

        else:

            self.__t = t
            self.__v = self.fluid_solver.v(t, p)

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

                h_ideal = self.cp_ideal * self.__t
                h_dep = self.int_t_dpdt + self.p * self.__v - self.fluid_solver.r_spc * self.__t
                self.__h = h_ideal + h_dep

        return self.__h

    @property
    def s(self):

        if self.__s is None:

            if self.bifase:

                self.__s = self.__liquid_state.s * (1 - self.__x) + self.__vapour_state.s * self.__x

            else:

                s_ideal = self.cp_ideal * np.log(self.__t) - self.fluid_solver.r_spc * np.log(self.p)
                s_dep = - self.int_rv
                self.__s = s_ideal + s_dep

        return self.__s

    @property
    def r(self):

        if self.__r is None:

            self.__r = (- self.__t * self.dpdt ** 2 / self.dpdv)

        return self.__r

    @property
    def cp(self):

        if self.bifase:

            return np.inf

        if self.__cp is None:

            self.__cp = self.cp_ideal + self.r - self.fluid_solver.r_spc + self.int_t_ddpddt

        return self.__cp

    @property
    def cp_ideal(self):

        if self.__cp_ideal is None:

            self.__cp_ideal = self.fluid_solver.cp_ideal(self.t)

        return self.__cp_ideal

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

            b = self.fluid_solver.b
            r1 = self.fluid_solver.r_1
            r2 = self.fluid_solver.r_2

            a_0 = 1 / (self.eta - 1) ** 2
            a_1 = 2 * self.eta + (r1 + r2)
            a_2 = ((self.eta - r1) * (self.eta - r2))**2
            a_3 = self.fluid_solver.r_spc * self.t / b ** 2

            self.__dpdv = (a_1 / a_2 * self.alpha - a_0) * a_3

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
            r_1 = self.fluid_solver.r_1
            r_2 = self.fluid_solver.r_2

            try:

                if r_1 == r_2:

                    self.__u = 1 / (- self.fluid_solver.b * (eta - r_1))

                else:

                    self.__u = np.log((eta - r_1) / (eta - r_2)) / (self.fluid_solver.b * (r_1 - r_2))

            except:

                self.__u = np.inf

        return self.__u

    @property
    def g(self):

        if self.__g is None:

            a_0 = self.fluid_solver.r_spc * self.__t / self.p
            self.__g = self.fluid_solver.r_spc * np.log(a_0 / (self.__v - self.fluid_solver.b))

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

            pass

        return self
