"""
2017/05/03
Solve uncertainties problem by adding new function `isuncertainties`
"""
import numpy as np
import uncertainties as uct
from uncertainties import unumpy as unp
from scipy.optimize import brenth
from scipy.misc import derivative
from .etc import isuncertainties


def vinet_p(v, v0, k0, k0p):
    """
    calculate pressure from vinet equation

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :return: pressure in GPa
    """
    # unumpy.exp works for both numpy and unumpy
    # so I set uncertainty default.
    # if unumpy.exp is used for lmfit, it generates an error
    return cal_p_vinet(v, [v0, k0, k0p],
                       uncertainties=isuncertainties([v, v0, k0, k0p]))


def cal_p_vinet(v, k, uncertainties=True):
    """
    calculate pressure from vinet equation

    :param v: volume at different pressures
    :param k: [v0, k0, k0p]
    :return: static pressure
    :note: internal function
    """
    v0 = k[0]
    k0 = k[1]
    k0p = k[2]
    x = np.power(v / v0, 1. / 3.)
    f1 = (1. - x) / (np.power(x, 2.))
    # np to unp for exp
    if uncertainties:
        f2 = unp.exp(1.5 * (k0p - 1.) * (1. - x))
    else:
        f2 = np.exp(1.5 * (k0p - 1.) * (1. - x))
    p = 3. * k0 * f1 * f2
    return p


def vinet_v_single(p, v0, k0, k0p, min_strain=0.01):
    """
    find volume at given pressure using brenth in scipy.optimize
    this is for single p value, not vectorized

    :param p: pressure in GPa
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :param min_strain: defining minimum v/v0 value to search volume for
    :return: unit cell volume at high pressure in A^3
    """
    if p <= 1.e-5:
        return v0

    def f_diff(v, v0, k0, k0p, p):
        return vinet_p(v, v0, k0, k0p) - p
    v = brenth(f_diff, v0, v0 * min_strain, args=(v0, k0, k0p, p))
    return v


def vinet_v(p, v0, k0, k0p, min_strain=0.01):
    """
    find volume at given pressure

    :param p: pressure in GPa
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :param min_strain: defining minimum v/v0 value to search volume for
    :return: unit cell volume at high pressure in A^3
    :note: wrapper function vetorizing vinet_v_single
    """
    if isuncertainties([p, v0, k0, k0p]):
        f_u = np.vectorize(uct.wrap(vinet_v_single), excluded=[1, 2, 3, 4])
        return f_u(p, v0, k0, k0p, min_strain=min_strain)
    else:
        f_v = np.vectorize(vinet_v_single, excluded=[1, 2, 3, 4])
        return f_v(p, v0, k0, k0p, min_strain=min_strain)


def cal_v_vinet(p, k):
    """
    calculate volume from vinet equation. wrapper for vinet_v

    :param p: pressure in GPa
    :param k: [v0, k0, k0p]
    :return: unit cell volume at high pressure in A^3
    :note: internal function
    """
    return vinet_v(p, k[0], k[1], k[2])


def vinet_k(p, v0, k0, k0p, numerical=False):
    """
    calculate bulk modulus, wrapper for cal_k_vinet
    cannot handle uncertainties

    :param p: pressure in GPa
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :return: bulk modulus at high pressure in GPa
    """
    f_u = uct.wrap(cal_k_vinet)
    return f_u(p, [v0, k0, k0p])


def vinet_dPdV(v, v0, k0, k0p, precision=1.e-5):
    """
    calculate dP/dV for numerical calculation of bulk modulus
    according to test this differs from analytical result by 1.e-5

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :param precision: precision for numerical calc (default = 1.e-5 * v0)
    :return: dP/dV
    """
    def f_scalar(v, v0, k0, k0p, precision=1.e-5):
        return derivative(vinet_p, v, args=(v0, k0, k0p), dx=v0 * precision)
    f_v = np.vectorize(f_scalar, excluded=[1, 2, 3, 4])
    return f_v(v, v0, k0, k0p, precision=precision)


def vinet_k_num(v, v0, k0, k0p, precision=1.e-5):
    """
    calculate bulk modulus numerically from volume, not pressure
    according to test this differs from analytical result by 1.e-5

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :param precision: precision for numerical calc (default = 1.e-5 * v0)
    :return: dP/dV
    """
    return -1. * v * vinet_dPdV(v, v0, k0, k0p, precision=precision)


def cal_k_vinet(p, k):
    """
    calculate bulk modulus in GPa

    :param p: pressure in GPa
    :param k: [v0, k0, k0p]
    :return: bulk modulus at high pressure in GPa
    """
    v = cal_v_vinet(p, k)
    return cal_k_vinet_from_v(v, k[0], k[1], k[2])


def cal_k_vinet_from_v(v, v0, k0, k0p):
    """
    calculate bulk modulus in GPa

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :return: bulk modulus at high pressure in GPa
    """
    x = v / v0
    y = np.power(x, 1. / 3.)
    eta = 1.5 * (k0p - 1.)
    k = k0 * np.power(y, -2.) * (1. + (eta * y + 1.) * (1. - y)) * \
        unp.exp((1. - y) * eta)
    return k
