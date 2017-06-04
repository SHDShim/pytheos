"""
Todo's
"""
from scipy import constants
import numpy as np
from scipy.optimize import brenth
from uncertainties import unumpy as unp
from ..eqn_bm3 import bm3_p
from ..eqn_vinet import vinet_p
from ..eqn_kunc import kunc_p
from ..eqn_therm import alphakt_pth
from ..eqn_therm_constq import constq_pth
from ..eqn_therm_Tange import tange_pth
from ..eqn_therm_Speziale import speziale_pth
from ..eqn_hugoniot import hugoniot_p, hugoniot_t
from ..eqn_jamieson import hugoniot_p_nlin
from ..eqn_therm_Dorogokupets2007 import dorogokupets2007_pth
from ..eqn_therm_Dorogokupets2015 import dorogokupets2015_pth
from ..eqn_anharmonic import zharkov_panh
from ..eqn_electronic import zharkov_pel, tsuchiya_pel
from ..conversion import vol_uc2mol, vol_mol2uc

func_st = {'bm3': bm3_p, 'vinet': vinet_p, 'kunc': kunc_p}
func_th = {'constq': constq_pth, 'tange': tange_pth,
           'speziale': speziale_pth, 'dorogokupets2007': dorogokupets2007_pth,
           'dorogokupets2015': dorogokupets2015_pth, 'alphakt': alphakt_pth}
func_el = {'zharkov': zharkov_pel, 'tsuchiya': tsuchiya_pel}
func_anh = {'zharkov': zharkov_panh}


class MGEOS(object):
    """
    All EOS following the Mie-Gruneisen fomulation
    """

    def __init__(self, n, z, params_st=None, params_th=None, params_anh=None,
                 params_el=None, eqn_st='bm3', eqn_th=None,
                 eqn_anh=None, eqn_el=None, t_ref=300.,
                 three_r=3. * constants.R, reference=None):
        """
        :param params_st: elastic parameters for static EOS in an OrderedDict
            [v0 in A^3, k0 in GPa, k0p]
        :param params_th: thermal parameters for thermal EOS in an OrderedDict.
            The component can differ depending on the equation used.
        :param params_anh: anharmonic parameters for anharmonic correction in
            an OrderedDict.  The component can differ depending on the
            equation used.
        :param params_el: electronic parameters for electronic correction in
            an OrderedDict. The component can differ depending on the
            equation used.
        :param eqn_st: equation type for the static EOS. 'bm3', 'vinet', or
            'kunc'
        :param eqn_th: equation type for the thermal EOS. 'constq', 'tange',
            'speziale', 'dorogokupets2007', 'dorogokupets2015', 'alphakt'
        :param eqn_anh: equation type for anharmonic correction. 'zharkov',
        :param eqn_el: equation type for electonic correction. 'zharkov',
            'tsuchiya'
        :param t_ref: reference temperature, 300 K
        :param three_r: 3 times gas constant.
            Jamieson modified this value to compensate for mismatches
        :param reference: reference for the EOS
        """
        self.params_st = params_st
        self.params_th = params_th
        self.params_anh = params_anh
        self.params_el = params_el
        self.eqn_st = eqn_st
        self.eqn_th = eqn_th
        self.eqn_el = eqn_el
        self.eqn_anh = eqn_anh
        self.n = n
        self.z = z
        self.three_r = three_r
        self.t_ref = t_ref
        self.reference = reference
        self.force_norm = False

    def print_reference(self):
        """
        show reference for the EOS
        """
        print("Ref: ", self.reference)

    def print_equations(self):
        """
        show equations used for the EOS
        """
        print("P_static: ", self.eqn_st)
        print("P_thermal: ", self.eqn_th)
        print("P_anharmonic: ", self.eqn_anh)
        print("P_electronic: ", self.eqn_el)

    def print_parameters(self):
        """
        show thermoelastic parameters for the EOS
        """
        print("Static: ", self.params_st)
        print("Thermal: ", self.params_th)
        print("Anharmonic: ", self.params_anh)
        print("Electronic: ", self.params_el)

    def _set_params(self, p):
        """
        change parameters in OrderedDict to list with or without uncertainties

        :param p: parameters in OrderedDict
        :return: parameters in list
        :note: internal function
        """
        if self.force_norm:
            params = [value.n for key, value in p.items()]
        else:
            params = [value for key, value in p.items()]
        return params

    def cal_pst(self, v):
        """
        calculate static pressure at 300 K.

        :param v: unit-cell volume in A^3
        :return: static pressure at t_ref (=300 K) in GPa
        """
        params = self._set_params(self.params_st)
        return func_st[self.eqn_st](v, *params)

    def cal_pth(self, v, temp):
        """
        calculate thermal pressure

        :param v: unit-cell volume in A^3
        :param temp: temperature in K
        :return: thermal pressure in GPa
        """
        if (self.eqn_th is None) or (self.params_th is None):
            return np.zeros_like(v)
        params = self._set_params(self.params_th)
        return func_th[self.eqn_th](v, temp, *params,
                                    self.n, self.z, t_ref=self.t_ref,
                                    three_r=self.three_r)

    def cal_pel(self, v, temp):
        """
        calculate pressure from electronic contributions

        :param v: unit-cell volume in A^3
        :param temp: temperature in K
        :return: pressure in GPa
        """
        if (self.eqn_el is None) or (self.params_el is None):
            return np.zeros_like(v)
        params = self._set_params(self.params_el)
        return func_el[self.eqn_el](v, temp, *params,
                                    self.n, self.z, t_ref=self.t_ref,
                                    three_r=self.three_r)

    def cal_panh(self, v, temp):
        """
        calculate pressure from anharmonic contributions

        :param v: unit-cell volume in A^3
        :param temp: temperature in K
        :return: pressure in GPa
        """
        if (self.eqn_anh is None) or (self.params_anh is None):
            return np.zeros_like(v)
        params = self._set_params(self.params_anh)
        return func_anh[self.eqn_anh](v, temp, *params,
                                      self.n, self.z, t_ref=self.t_ref,
                                      three_r=self.three_r)

    def cal_p(self, v, temp):
        """
        calculate total pressure at given volume and temperature

        :param v: unit-cell volume in A^3
        :param temp: temperature in K
        :return: pressure in GPa
        :note: 2017/05/10 temp must be numpy array.  If not, such as list,
            create an error.
        """
        if not isinstance(temp, np.ndarray):
            temp = np.asarray(temp)
        if not isinstance(v, np.ndarray):
            v = np.asarray(v)
        return self.cal_pst(v) + self.cal_pth(v, temp) + \
            self.cal_pel(v, temp) + self.cal_panh(v, temp)

    def cal_v(self, p, temp, min_strain=0.2, max_strain=1.0):
        """
        calculate unit-cell volume at given pressure and temperature

        :param p: pressure in GPa
        :param temp: temperature in K
        :param min_strain: minimum strain searched for volume root
        :param max_strain: maximum strain searched for volume root
        :return: unit-cell volume in A^3
        :note: 2017/05/10 I found wrap function is not compatible with
            OrderedDict. So I convert unp array to np array.
        """
        v0 = self.params_st['v0'].nominal_value
        self.force_norm = True
        pp = unp.nominal_values(p)
        ttemp = unp.nominal_values(temp)

        def _cal_v_single(pp, ttemp):
            if (pp <= 1.e-5) and (ttemp == 300.):
                return v0

            def f_diff(v, ttemp, pp):
                return self.cal_p(v, ttemp) - pp
            # print(f_diff(v0 * 0.3, temp, p))
            v = brenth(f_diff, v0 * max_strain, v0 * min_strain,
                       args=(ttemp, pp))
            return v
        f_vu = np.vectorize(_cal_v_single)
        v = f_vu(pp, ttemp)
        self.force_norm = False
        return v


class JHEOS(object):
    """
    Jamieson's hugoniot EOS.  The equations are from Jamieson 1982
    """

    def __init__(self, n, z, mass, params_hugoniot, params_therm,
                 three_r=3. * constants.R, c_v=0.0, nonlinear=False,
                 reference=None, t_ref=300.):
        """
        :param n: number of elements in a chemical formula
        :param z: number of formula unit in a unit cell
        :param m: molar mass in gram
        :param params_hugoniot: hugoniot parameters.
            [rho0 in g/cm3, c0 in km/s, s] for linear case and
            [rho0 in g/cm3, a in km/s, b, c in s/km] for non-linear case.
            See Jamieson for detail.
        :param params_therm: thermal parameters for the constq equations
            [v0 in A3, gamma0, q, theta0 in K]
        :param three_r: 3 times gas constant.
            Jamieson modified this value to compensate for mismatches
        :param c_v: heat capacity for hugoniot temperature calculation
        :param nonlinear: nonlinear Us-Up fit particularly for Jamieson's gold
            scale
        :param reference: reference for the EOS
        :param t_ref: reference temperature, 300 K
        """
        self.params_hugoniot = params_hugoniot
        self.params_therm = params_therm
        self.n = n
        self.z = z
        self.mass = mass
        self.ef_v = 0.001
        self.ef_T = 0.05
        self.three_r = three_r
        self.nonlinear = nonlinear
        self.c_v = c_v
        self.reference = reference
        self.t_ref = t_ref
        self.force_norm = False

    def _set_params(self, p):
        """
        change parameters in OrderedDict to list with or without uncertainties

        :param p: parameters in OrderedDict
        :return: parameters in list
        :note: internal function
        """
        if self.force_norm:
            params = [value.n for key, value in p.items()]
        else:
            params = [value for key, value in p.items()]
        return params

    def print_reference(self):
        """
        show reference for the EOS
        """
        print("Ref: ", self.reference)

    def print_equations(self):
        """
        show equations used for the EOS
        """
        print("P_static: Jamieson")
        print("P_thermal: Constq")
        print("P_anharmonic: None")
        print("P_electronic: None")

    def print_parameters(self):
        """
        show thermoelastic parameters for the EOS
        """
        print("Static: ", self.params_hugoniot)
        print("Thermal: ", self.params_therm)
        print("Anharmonic: None")
        print("Electronic: None")

    def _get_rho(self, v):
        """
        convert unit-cell volume in A^3 to density in g/cm^3

        :param v: unit cell volume in A^3
        :return: density in g/cm^3

        :note: internal function
        """
        v_mol = vol_uc2mol(v, self.z)  # in m^3
        rho = self.mass / v_mol * 1.e-6  # in g/cm^3
        return rho

    def _get_v(self, rho):
        """
        convert density in g/cm^3 to unit-cell volume

        :param rho: density in g/cm^3
        :return: unit-cell volume in A^3

        :note: internal function
        """
        v_mol = self.mass / rho  # in cm^3
        v = vol_mol2uc(v_mol * 1.e-6, self.z)
        print(v)
        return v

    def cal_pst(self, v):
        """
        calculate static pressure at 300 K.

        :param v: unit-cell volume in A^3
        :return: static pressure at t_ref (=300 K) in GPa
        """
        p_h = self._hugoniot_p(v)
        p_th_h = self._hugoniot_pth(v)
        p_st = p_h - p_th_h
        return p_st

    def _hugoniot_p(self, v):
        """
        calculate static pressure at 300 K.

        :param v: unit-cell volume in A^3
        :return: static pressure at t_ref (=300 K) in GPa
        """
        rho = self._get_rho(v)
        params = self._set_params(self.params_hugoniot)
        if self.nonlinear:
            return hugoniot_p_nlin(rho, *params)
        else:
            return hugoniot_p(rho, *params)

    def _hugoniot_t(self, v):
        """
        calculate Hugoniot temperature

        :param v: unit-cell volume in A^3
        :return: temperature in K
        :note: 2017/05/10, It is intentional that I call hugoniot_t
        instead of hugoniot_t_nlin for the nonlinear case.
        The reason is the hugoniot_t_nlin (and its equivalent part in mdaap)
        has numerical problem.  In fact, I found no problem to match the
        Jamieson table values with the use of linear thermal part instead of
        non-linear version.
        """
        rho = self._get_rho(v)
        params_h = self._set_params(self.params_hugoniot)
        params_t = self._set_params(self.params_therm)
        if self.nonlinear:
            return hugoniot_t(rho, *params_h[:-1], *params_t[1:],
                              self.n, self.mass, three_r=self.three_r,
                              c_v=self.c_v)
        else:
            return hugoniot_t(rho, *params_h, *params_t[1:],
                              self.n, self.mass, three_r=self.three_r,
                              c_v=self.c_v)

    def cal_pth(self, v, temp):
        """
        calculate thermal pressure

        :param v: unit-cell volume in A^3
        :param temp: temperature in K
        :return: thermal pressure in GPa
        """
        params_t = self._set_params(self.params_therm)
        return constq_pth(v, temp, *params_t, self.n, self.z,
                          t_ref=self.t_ref, three_r=self.three_r)

    def _hugoniot_pth(self, v):
        """
        calculate thermal pressure along hugoniot

        :param v: unit-cell volume in A^3
        :return: thermal pressure along hugoniot in GPa
        """
        temp = self._hugoniot_t(v)
        return self.cal_pth(v, temp)

    def cal_p(self, v, temp):
        """
        calculate total pressure at given volume and temperature

        :param v: unit-cell volume in A^3
        :param temp: temperature in K
        :return: pressure in GPa
        """
        return self.cal_pst(v) + self.cal_pth(v, temp)

    def cal_v(self, p, temp, min_strain=0.3, max_strain=1.0):
        """
        calculate unit-cell volume at given pressure and temperature

        :param p: pressure in GPa
        :param temp: temperature in K
        :param min_strain: minimum strain searched for volume root
        :param max_strain: maximum strain searched for volume root
        :return: unit-cell volume in A^3
        :note: 2017/05/10 I found wrap function is not compatible with
            OrderedDict. So I convert unp array to np array.
        """
        v0 = self.params_therm['v0'].nominal_value
        self.force_norm = True
        pp = unp.nominal_values(p)
        ttemp = unp.nominal_values(temp)

        def _cal_v_single(pp, ttemp):
            if (pp <= 1.e-5) and (ttemp == 300.):
                return v0

            def f_diff(v, ttemp, pp):
                return self.cal_p(v, ttemp) - pp
            # print(f_diff(v0 * 0.3, temp, p))
            v = brenth(f_diff, v0 * max_strain, v0 * min_strain,
                       args=(ttemp, pp))
            return v
        f_vu = np.vectorize(_cal_v_single)
        v = f_vu(pp, ttemp)
        self.force_norm = False
        return v
