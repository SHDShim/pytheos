"""
2017/05/03 I believe now I settle down with the uncertainties issue
through `isuncertainties` check.
"""
import numpy as np
from uncertainties import unumpy as unp
from scipy import constants
from .conversion import vol_uc2mol
from .eqn_debye import debye_E
from .etc import isuncertainties


def tange_grun(v, v0, gamma0, a, b):
    """
    calculate Gruneisen parameter for the Tange equation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param a: volume-independent adjustable parameters
    :param b: volume-independent adjustable parameters
    :return: Gruneisen parameter
    """
    x = v / v0
    return gamma0 * (1. + a * (np.power(x, b) - 1.))


def tange_debyetemp(v, v0, gamma0, a, b, theta0):
    """
    calculate Debye temperature for the Tange equation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param a: volume-independent adjustable parameters
    :param b: volume-independent adjustable parameters
    :param theta0: Debye temperature at 1 bar in K
    :return: Debye temperature in K
    """
    x = v / v0
    gamma = tange_grun(v, v0, gamma0, a, b)
    if isuncertainties([v, v0, gamma0, a, b, theta0]):
        theta = theta0 * np.power(x, (-1. * (1. - a) * gamma0)) * \
            unp.exp((gamma0 - gamma) / b)
    else:
        theta = theta0 * np.power(x, (-1. * (1. - a) * gamma0)) * \
            np.exp((gamma0 - gamma) / b)
    return theta


def tange_pth(v, temp, v0, gamma0, a, b, theta0, n, z,
              t_ref=300., three_r=3. * constants.R):
    """
    calculate thermal pressure for the Tange equation

    :param v: unit-cell volume in A^3
    :param temp: temperature in K
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param a: volume-independent adjustable parameters
    :param b: volume-independent adjustable parameters
    :param theta0: Debye temperature at 1 bar in K
    :param n: number of atoms in a formula unit
    :param z: number of formula unit in a unit cell
    :param t_ref: reference temperature
    :param three_r: 3R in case adjustment is needed
    :return: thermal pressure in GPa
    """
    v_mol = vol_uc2mol(v, z)
    gamma = tange_grun(v, v0, gamma0, a, b)
    theta = tange_debyetemp(v, v0, gamma0, a, b, theta0)
    xx = theta / temp
    debye = debye_E(xx)
    if t_ref == 0.:
        debye0 = 0.
    else:
        xx0 = theta / t_ref
        debye0 = debye_E(xx0)
    Eth0 = three_r * n * t_ref * debye0
    Eth = three_r * n * temp * debye
    delEth = Eth - Eth0
    p_th = (gamma / v_mol * delEth) * 1.e-9
    return p_th
