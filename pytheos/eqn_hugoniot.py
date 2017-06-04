import numpy as np
from scipy.integrate import odeint
import scipy.constants as constants
import uncertainties as uct
from scipy.optimize import brenth
from .eqn_debye import debye_E
from .etc import isuncertainties


def hugoniot_p(rho, rho0, c0, s):
    """
    calculate pressure along a Hugoniot

    :param rho: density in g/cm^3
    :param rho0: density at 1 bar in g/cm^3
    :param c0: velocity at 1 bar in km/s
    :param s: slope of the velocity change
    :return: pressure in GPa
    """
    eta = 1. - (rho0 / rho)
    Ph = rho0 * c0 * c0 * eta / np.power((1. - s * eta), 2.)
    return Ph


def _dT_h_delta(T_in_kK, eta, k, threenk, c_v):
    """
    internal function for calculation of temperature along a Hugoniot

    :param T_in_kK: temperature in kK scale, see Jamieson for detail
    :param eta: = 1 - rho0/rho
    :param k: = [rho0, c0, s, gamma0, q, theta0]
    :param threenk: see the definition in Jamieson 1983,
        it is a correction term mostly for Jamieson gold scale
    :param c_v: manual input of Cv value,
        if 0 calculated through Debye function
    :return: eta derivative of temperature
    """
    rho0 = k[0]  # g/m^3
    gamma0 = k[3]  # no unit
    q = k[4]  # no unit
    theta0_in_kK = k[5]  # K, see Jamieson 1983 for detail
    rho = rho0 / (1. - eta)
    c0 = k[1]  # km/s
    s = k[2]  # no unit
    dPhdelta_H = rho0 * c0 * c0 * (1. + s * eta) / \
        np.power((1. - s * eta), 3.)
    # [g/cm^3][km/s]^2 = 1e9[kg m^2/s^2] = [GPa]
    Ph = hugoniot_p(rho, rho0, c0, s)  # in [GPa]
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
    # print('dYdX', dYdX)
    return dYdX


def hugoniot_t_single(rho, rho0, c0, s, gamma0, q, theta0, n, mass,
                      three_r=3. * constants.R, t_ref=300., c_v=0.):
    """
    internal function to calculate pressure along Hugoniot

    :param rho: density in g/cm^3
    :param rho0: density at 1 bar in g/cm^3
    :param c0: velocity at 1 bar in km/s
    :param s: slope of the velocity change
    :param gamma0: Gruneisen parameter at 1 bar
    :param q: logarithmic derivative of Gruneisen parameter
    :param theta0: Debye temperature in K
    :param n: number of elements in a chemical formula
    :param mass: molar mass in gram
    :param three_r: 3 times gas constant.
        Jamieson modified this value to compensate for mismatches
    :param t_ref: reference temperature, 300 K
    :param c_v: heat capacity, see Jamieson 1983 for detail
    :return: temperature along hugoniot
    """
    eta = 1. - rho0 / rho
    if eta == 0.0:
        return 300.
    threenk = three_r / mass * n  # [J/mol/K] / [g/mol] = [J/g/K]
    k = [rho0, c0, s, gamma0, q, theta0 / 1.e3]
    t_h = odeint(_dT_h_delta, t_ref / 1.e3, [0., eta],
                 args=(k, threenk, c_v), full_output=1)
    temp_h = np.squeeze(t_h[0][1])
    return temp_h * 1.e3


def hugoniot_t(rho, rho0, c0, s, gamma0, q, theta0, n, mass,
               three_r=3. * constants.R, t_ref=300., c_v=0.):
    """
    calculate temperature along a hugoniot

    :param rho: density in g/cm^3
    :param rho0: density at 1 bar in g/cm^3
    :param c0: velocity at 1 bar in km/s
    :param s: slope of the velocity change
    :param gamma0: Gruneisen parameter at 1 bar
    :param q: logarithmic derivative of Gruneisen parameter
    :param theta0: Debye temperature in K
    :param n: number of elements in a chemical formula
    :param mass: molar mass in gram
    :param three_r: 3 times gas constant.
        Jamieson modified this value to compensate for mismatches
    :param t_ref: reference temperature, 300 K
    :param c_v: heat capacity, see Jamieson 1983 for detail
    :return: temperature along hugoniot
    """
    if isuncertainties([rho, rho0, c0, s, gamma0, q, theta0]):
        f_v = np.vectorize(uct.wrap(hugoniot_t_single),
                           excluded=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    else:
        f_v = np.vectorize(hugoniot_t_single, excluded=[1, 2, 3, 4, 5, 6,
                                                        7, 8, 9, 10, 11])
    return f_v(rho, rho0, c0, s, gamma0, q, theta0, n, mass,
               three_r=three_r, t_ref=t_ref, c_v=c_v)
#    return np.squeeze(T_h[0]) * 1.e3


def hugoniot_rho(p, rho0, c0, s, min_strain=0.01):
    """
    calculate density in g/cm^3 from a hugoniot curve

    :param p: pressure in GPa
    :param rho0: density at 1 bar in g/cm^3
    :param c0: velocity at 1 bar in km/s
    :param s: slope of the velocity change
    :param min_strain: defining minimum v/v0 value to search volume for
    :return: density in g/cm^3
    :note: this is a wrapper function by vectorizing hugoniot_rho_single
    """
    if isuncertainties([p, rho0, c0, s]):
        f_u = np.vectorize(uct.wrap(hugoniot_rho_single),
                           excluded=[1, 2, 3, 4])
        return f_u(p, rho0, c0, s, min_strain=min_strain)
    else:
        f_v = np.vectorize(hugoniot_rho_single, excluded=[1, 2, 3, 4])
        return f_v(p, rho0, c0, s, min_strain=min_strain)


def hugoniot_rho_single(p, rho0, c0, s, min_strain=0.01):
    """
    calculate density in g/cm^3 from a hugoniot curve

    :param p: pressure in GPa
    :param rho0: density at 1 bar in g/cm^3
    :param c0: velocity at 1 bar in km/s
    :param s: slope of the velocity change
    :param min_strain: defining minimum v/v0 value to search volume for
    :return: density in g/cm^3
    """
    if p <= 1.e-5:
        return rho0

    def f_diff(rho):
        return hugoniot_p(rho, rho0, c0, s) - p
    rho = brenth(f_diff, rho0, rho0 / min_strain)
    return rho
