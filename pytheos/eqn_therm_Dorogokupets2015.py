import numpy as np
from uncertainties import unumpy as unp
from scipy import constants
from .etc import isuncertainties
from .conversion import vol_uc2mol
from .eqn_therm_Dorogokupets2007 import altshuler_grun, altshuler_debyetemp


def dorogokupets2015_pth(v, temp, v0, gamma0, gamma_inf, beta, theta01,
                         m1, theta02, m2, n, z, t_ref=300.,
                         three_r=3. * constants.R):
    """
    calculate thermal pressure for Dorogokupets 2015 EOS

    :param v: unit-cell volume in A^3
    :param temp: temperature in K
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param gamma_inf: Gruneisen parameter at infinite pressure
    :param beta: volume dependence of Gruneisen parameter
    :param theta01: Debye temperature at 1 bar in K
    :param m1: weighting factor, see Dorogokupets 2015 for detail
    :param theta02: Debye temperature at 1 bar in K
    :param m2: weighting factor, see Dorogokupets 2015 for detail
    :param n: number of elements in a chemical formula
    :param z: number of formula unit in a unit cell
    :param three_r: 3 times gas constant.
        Jamieson modified this value to compensate for mismatches
    :param t_ref: reference temperature, 300 K
    :return: thermal pressure in GPa
    """
    # x = v / v0
    # a = a0 * np.power(x, m)
    v_mol = vol_uc2mol(v, z)
    gamma = altshuler_grun(v, v0, gamma0, gamma_inf, beta)
    theta1 = altshuler_debyetemp(v, v0, gamma0, gamma_inf, beta, theta01)
    theta2 = altshuler_debyetemp(v, v0, gamma0, gamma_inf, beta, theta02)
    if isuncertainties([v, temp, v0, gamma0, gamma_inf, beta,
                        theta01, m1, theta02, m2]):
        term_h1 = m1 / (m1 + m2) * three_r * n * gamma / v_mol * \
            (theta1 / (unp.exp(theta1 / temp) - 1.))
        term_h2 = m2 / (m1 + m2) * three_r * n * gamma / v_mol * \
            (theta2 / (unp.exp(theta2 / temp) - 1.))
        term_h1_ref = m1 / (m1 + m2) * three_r * n * gamma / v_mol * \
            (theta1 / (unp.exp(theta1 / t_ref) - 1.))
        term_h2_ref = m2 / (m1 + m2) * three_r * n * gamma / v_mol * \
            (theta2 / (unp.exp(theta2 / t_ref) - 1.))
    else:
        term_h1 = m1 / (m1 + m2) * three_r * n * gamma / v_mol * \
            (theta1 / (np.exp(theta1 / temp) - 1.))
        term_h2 = m2 / (m1 + m2) * three_r * n * gamma / v_mol * \
            (theta2 / (np.exp(theta2 / temp) - 1.))
        term_h1_ref = m1 / (m1 + m2) * three_r * n * gamma / v_mol * \
            (theta1 / (np.exp(theta1 / t_ref) - 1.))
        term_h2_ref = m2 / (m1 + m2) * three_r * n * gamma / v_mol * \
            (theta2 / (np.exp(theta2 / t_ref) - 1.))
    p_th = term_h1 * 1.e-9 + term_h2 * 1.e-9
    p_th_ref = term_h1_ref * 1.e-9 + term_h2_ref * 1.e-9
    return (p_th - p_th_ref)
