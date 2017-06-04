"""
2017/05/18 At the moment the functions here are not being revealed
    to the users.  I need to check if any of these functions are being
    called from other functions in pytheos.
"""
import scipy.constants as constants
import numpy as np
from scipy.integrate import odeint
import uncertainties as uct
from uncertainties import unumpy as unp
from .eqn_therm_constq import constq_pth
from .eqn_hugoniot import hugoniot_p, hugoniot_t
from .conversion import vol_uc2mol
from .eqn_debye import debye_E
from .etc import isuncertainties

"""
def _get_rho(self, v):
    v_mol = v / self.z * constants.Avogadro / 1.e24
    rho = self.mass / v_mol  # g/cm^3
    return rho

def _get_v(self, rho):
    v_mol = self.mass / rho
    v = v_mol * self.z / constants.Avogadro * 1.e24
    return v
"""


def jamieson_pst(v, v0, c0, s, gamma0, q, theta0, n, z, mass, c_v,
                 three_r=3. * constants.R, t_ref=300.):
    """
    calculate static pressure at 300 K from Hugoniot data using the constq
    formulation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param c0: velocity at 1 bar in km/s
    :param s: slope of the velocity change
    :param gamma0: Gruneisen parameter at 1 bar
    :param q: logarithmic derivative of Gruneisen parameter
    :param theta0: Debye temperature in K
    :param n: number of elements in a chemical formula
    :param z: number of formula unit in a unit cell
    :param mass: molar mass in gram
    :param c_v: heat capacity
    :param three_r: 3 times gas constant.
        Jamieson modified this value to compensate for mismatches
    :param t_ref: reference temperature, 300 K
    :return: static pressure in GPa
    :note: 2017/05/18 I am unsure if this is actually being used in pytheos
    """
    rho = mass / vol_uc2mol(v, z) * 1.e-6
    rho0 = mass / vol_uc2mol(v0, z) * 1.e-6
    p_h = hugoniot_p(rho, rho0, c0, s)
    p_th_h = jamieson_pth(v, v0, c0, s, gamma0, q, theta0, n, z, mass, c_v,
                          three_r=three_r, t_ref=t_ref)
    p_st = p_h - p_th_h
    return p_st


def jamieson_pth(v, v0, c0, s, gamma0, q, theta0, n, z, mass, c_v,
                 three_r=3. * constants.R, t_ref=300.):
    """
    calculate thermal pressure from Hugoniot data using the constq
    formulation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param c0: velocity at 1 bar in km/s
    :param s: slope of the velocity change
    :param gamma0: Gruneisen parameter at 1 bar
    :param q: logarithmic derivative of Gruneisen parameter
    :param theta0: Debye temperature in K
    :param n: number of elements in a chemical formula
    :param z: number of formula unit in a unit cell
    :param mass: molar mass in gram
    :param c_v: heat capacity
    :param three_r: 3 times gas constant.
        Jamieson modified this value to compensate for mismatches
    :param t_ref: reference temperature, 300 K
    :return: static pressure in GPa
    :note: 2017/05/18 I am unsure if this is actually being used in pytheos
    """
    rho = mass / vol_uc2mol(v, z) * 1.e-6
    rho0 = mass / vol_uc2mol(v0, z) * 1.e-6
    temp = hugoniot_t(rho, rho0, c0, s, gamma0, q, theta0, n, mass,
                      three_r=three_r, t_ref=t_ref, c_v=c_v)
    pth = constq_pth(v, temp, v0, gamma0, q, theta0, n, z, t_ref=t_ref,
                     three_r=three_r)
    return pth


def hugoniot_p_nlin(rho, rho0, a, b, c):
    """
    calculate pressure along a Hugoniot throug nonlinear equations
    presented in Jameison 1982

    :param rho: density in g/cm^3
    :param rho0: density at 1 bar in g/cm^3
    :param a: prefactor for nonlinear fit of Hugoniot data
    :param b: prefactor for nonlinear fit of Hugoniot data
    :param c: prefactor for nonlinear fit of Hugoniot data
    :return: pressure along Hugoniot in GPa
    """
    eta = 1. - (rho0 / rho)
    Up = np.zeros_like(eta)
    if isuncertainties([rho, rho0, a, b, c]):
        Up[eta != 0.] = ((b * eta - 1.) + unp.sqrt(
            np.power((1. - b * eta), 2.) - 4. * np.power(eta, 2.) * a * c)) /\
            (-2. * eta * c)
    else:
        Up[eta != 0.] = ((b * eta - 1.) + np.sqrt(
            np.power((1. - b * eta), 2.) - 4. * np.power(eta, 2.) * a * c)) /\
            (-2. * eta * c)
    Us = a + Up * b + Up * Up * c
    Ph = rho0 * Up * Us
    return Ph


"""
The lines below do not work.  The problem is in dUp_delta.  But without this
approach, I can still use the linear case to make similar value as Jamieson.
In fact, in mdaap this part is there but never been called because of Internal
bug.  At the moment, I do not want to struggle with this anymore as the value
I got is decent.
def _dT_h_delta_nlin(T_in_kK, eta, k, threenk, c_v):

    rho0 = k[0]  # g/m^3
    a = k[1]
    b = k[2]
    c = k[3]
    gamma0 = k[4]  # no unit
    q = k[5]  # no unit
    theta0_in_kK = k[6]  # K, see Jamieson 1983 for detail
    rho = rho0 / (1. - eta)
    Ph = hugoniot_p_nlin(rho, rho0, a, b, c)  # in [GPa]
    Up = ((b * eta - 1.) + np.sqrt(
        np.power((1. - b * eta), 2.) - 4. * np.power(eta, 2.) * a * c)) / \
        (-2. * eta * c)
    Us = Up * Up * c + Up * b + a
    dUp_delta = Us / (1. - b * eta - 2. * c * eta * Up)
    dUsdUp = b + 2. * c * Up
    dPhdelta_H = rho0 * Us * dUp_delta + rho0 * Up * dUsdUp * dUp_delta

    # calculate Cv
    gamma = gamma0 * np.power((1. - eta), q)
    theta_in_kK = theta0_in_kK * np.exp((gamma0 - gamma) / q)
    x = theta_in_kK / T_in_kK
    debye3 = debye_E(x)
    if c_v == 0.:
        c_v = threenk * (4. * debye3 - 3. * x / (np.exp(x) - 1.))  # [J/g/K]
    # calculate dYdX
    dYdX = (gamma / (1. - eta) * T_in_kK) + (dPhdelta_H * eta - Ph) / \
        (2. * c_v * rho0)
    return dYdX


def hugoniot_t_nlin_single(rho, rho0, a, b, c, gamma0, q, theta0, n, mass,
                           three_r=3. * constants.R, t_ref=300., c_v=0.):
    eta = 1. - rho0 / rho
    if abs(eta) <= 1.e-3:
        return 300.
    threenk = three_r / mass * n  # [J/mol/K] / [g/mol] = [J/g/K]
    k = [rho0, a, b, c, gamma0, q, theta0 / 1.e3]
    t_h = odeint(_dT_h_delta_nlin, t_ref / 1.e3, [0., eta],
                 args=(k, threenk, c_v), full_output=1)
    temp_h = np.squeeze(t_h[0][1])
    return temp_h * 1.e3


def hugoniot_t_nlin(rho, rho0, a, b, c, gamma0, q, theta0, n, mass,
                    three_r=3. * constants.R, t_ref=300., c_v=0.):
    if isuncertainties([rho, rho0, a, b, c, gamma0, q, theta0]):
        print('uncertainties case')
        f_v = np.vectorize(uct.wrap(hugoniot_t_nlin_single),
                           excluded=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    else:
        f_v = np.vectorize(hugoniot_t_nlin_single, excluded=[1, 2, 3, 4, 5, 6,
                                                             7, 8, 9, 10, 11,
                                                             12])
    return f_v(rho, rho0, a, b, c, gamma0, q, theta0, n, mass,
               three_r=three_r, t_ref=t_ref, c_v=c_v)
#    return np.squeeze(T_h[0]) * 1.e3
"""
