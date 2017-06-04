import numpy as np
from uncertainties import unumpy as unp
from .eqn_debye import debye_E
from .conversion import vol_uc2mol
import scipy.constants as constants
from .etc import isuncertainties


def altshuler_grun(v, v0, gamma0, gamma_inf, beta):
    """
    calculate Gruneisen parameter for Altshuler equation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param gamma_inf: Gruneisen parameter at infinite pressure
    :param beta: volume dependence of Gruneisen parameter
    :return: Gruneisen parameter
    """
    x = v / v0
    return gamma_inf + (gamma0 - gamma_inf) * np.power(x, beta)


def altshuler_debyetemp(v, v0, gamma0, gamma_inf, beta, theta0):
    """
    calculate Debye temperature for Altshuler equation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param gamma_inf: Gruneisen parameter at infinite pressure
    :param beta: volume dependence of Gruneisen parameter
    :param theta0: Debye temperature at 1 bar in K
    :return: Debye temperature in K
    """
    x = v / v0
    if isuncertainties([v, v0, gamma0, gamma_inf, beta, theta0]):
        theta = theta0 * np.power(x, -1. * gamma_inf) *\
            unp.exp((gamma0 - gamma_inf) / beta * (1. - np.power(x, beta)))
    else:
        theta = theta0 * np.power(x, -1. * gamma_inf) *\
            np.exp((gamma0 - gamma_inf) / beta * (1. - np.power(x, beta)))
    return theta


def dorogokupets2007_pth(v, temp, v0, gamma0, gamma_inf, beta, theta0, n, z,
                         three_r=3. * constants.R, t_ref=300.):
    """
    calculate thermal pressure for Dorogokupets 2007 EOS

    :param v: unit-cell volume in A^3
    :param temp: temperature in K
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param gamma_inf: Gruneisen parameter at infinite pressure
    :param beta: volume dependence of Gruneisen parameter
    :param theta0: Debye temperature at 1 bar in K
    :param n: number of elements in a chemical formula
    :param z: number of formula unit in a unit cell
    :param three_r: 3 times gas constant.
        Jamieson modified this value to compensate for mismatches
    :param t_ref: reference temperature, 300 K
    :return: thermal pressure in GPa
    """
    v_mol = vol_uc2mol(v, z)
    # x = v_mol / v0_mol
    gamma = altshuler_grun(v, v0, gamma0, gamma_inf, beta)
    theta = altshuler_debyetemp(v, v0, gamma0, gamma_inf, beta, theta0)

    def f(t):
        xx = theta / t
        debye = debye_E(xx)
        Eth = three_r * n * t * debye
        return (gamma / v_mol * Eth) * 1.e-9
    return f(temp) - f(t_ref)
