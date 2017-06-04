"""
2017/05/03 I believe now I settle down with the uncertainties issue
through `isuncertainties` check.
"""
import numpy as np
from uncertainties import unumpy as unp
import uncertainties as uct
from scipy import constants
from scipy.integrate import quad
from .conversion import vol_uc2mol
from .eqn_debye import debye_E
from .etc import isuncertainties


def speziale_grun(v, v0, gamma0, q0, q1):
    """
    calculate Gruneisen parameter for the Speziale equation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param q0: logarithmic derivative of Gruneisen parameter
    :param q1: logarithmic derivative of Gruneisen parameter
    :return: Gruneisen parameter
    """
    if isuncertainties([v, v0, gamma0, q0, q1]):
        gamma = gamma0 * unp.exp(q0 / q1 * ((v / v0) ** q1 - 1.))
    else:
        gamma = gamma0 * np.exp(q0 / q1 * ((v / v0) ** q1 - 1.))
    return gamma


def speziale_debyetemp(v, v0, gamma0, q0, q1, theta0):
    """
    calculate Debye temperature for the Speziale equation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param q0: logarithmic derivative of Gruneisen parameter
    :param q1: logarithmic derivative of Gruneisen parameter
    :param theta0: Debye temperature at 1 bar in K
    :return: Debye temperature in K
    """
    if isuncertainties([v, v0, gamma0, q0, q1, theta0]):
        f_vu = np.vectorize(uct.wrap(integrate_gamma),
                            excluded=[1, 2, 3, 4, 5, 6])
        integ = f_vu(v, v0, gamma0, q0, q1, theta0)
        theta = unp.exp(unp.log(theta0) - integ)
    else:
        f_v = np.vectorize(integrate_gamma, excluded=[1, 2, 3, 4, 5, 6])
        integ = f_v(v, v0, gamma0, q0, q1, theta0)
        theta = np.exp(np.log(theta0) - integ)
    return theta


def integrate_gamma(v, v0, gamma0, q0, q1, theta0):
    """
    internal function to calculate Debye temperature

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param q0: logarithmic derivative of Gruneisen parameter
    :param q1: logarithmic derivative of Gruneisen parameter
    :param theta0: Debye temperature at 1 bar in K
    :return: Debye temperature in K
    """
    def f_integrand(v):
        gamma = gamma0 * np.exp(q0 / q1 * ((v / v0) ** q1 - 1.))
        return gamma / v

    theta_term = quad(f_integrand, v0, v)[0]
    return theta_term


def speziale_pth(v, temp, v0, gamma0, q0, q1, theta0, n, z, t_ref=300.,
                 three_r=3. * constants.R):
    """
    calculate thermal pressure for the Speziale equation

    :param v: unit-cell volume in A^3
    :param temp: temperature in K
    :param v0: unit-cell volume in A^3 at 1 bar
    :param gamma0: Gruneisen parameter at 1 bar
    :param q0: logarithmic derivative of Gruneisen parameter
    :param q1: logarithmic derivative of Gruneisen parameter
    :param theta0: Debye temperature at 1 bar in K
    :param n: number of atoms in a formula unit
    :param z: number of formula unit in a unit cell
    :param t_ref: reference temperature
    :param three_r: 3R in case adjustment is needed
    :return: thermal pressure in GPa
    """
    v_mol = vol_uc2mol(v, z)
    gamma = speziale_grun(v, v0, gamma0, q0, q1)
    theta = speziale_debyetemp(v, v0, gamma0, q0, q1, theta0)
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
