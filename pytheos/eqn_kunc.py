"""
2017/05/03 kunc_k is not available yet because I did not derive analytical
form yet.
"""
import numpy as np
import uncertainties as uct
from uncertainties import unumpy as unp
from scipy.optimize import brenth
from scipy.misc import derivative
from .etc import isuncertainties


def kunc_p(v, v0, k0, k0p, order=5):
    """
    calculate Kunc EOS
    see Dorogokupets 2015 for detail

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :param order: order for the Kunc equation
    :return: pressure in GPa
    """
    return cal_p_kunc(v, [v0, k0, k0p], order=order,
                      uncertainties=isuncertainties([v, v0, k0, k0p]))


def cal_p_kunc(v, k, order=5, uncertainties=True):
    """
    calculate Kunc EOS,
    see Dorogokupets2015 for functional form

    :param v: unit-cell volume in A^3
    :param k: [v0, k0, k0p]
    :param order: order for the Kunc equation
    :param uncertainties: use of uncertainties package
    :return: pressure in GPa
    :note: internal function
    """
    v0 = k[0]
    k0 = k[1]
    k0p = k[2]
    x = np.power(v / v0, 1. / 3.)
    f1 = (1. - x) / (np.power(x, order))
    if uncertainties:
        f2 = unp.exp((1.5 * k0p - order + 0.5) * (1. - x))
    else:
        f2 = np.exp((1.5 * k0p - order + 0.5) * (1. - x))
    p = 3. * k0 * f1 * f2
    return p


def kunc_v_single(p, v0, k0, k0p, order=5, min_strain=0.01):
    """
    find volume at given pressure using brenth in scipy.optimize

    :param p: pressure in GPa
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :param order: order of Kunc function
    :param min_strain: defining minimum v/v0 value to search volume for
    :return: unit-cell volume at high pressure in GPa
    """
    if p <= 1.e-5:
        return v0

    def f_diff(v, v0, k0, k0p, p, order=order):
        return kunc_p(v, v0, k0, k0p, order=order) - p
    v = brenth(f_diff, v0, v0 * min_strain, args=(v0, k0, k0p, p, order))
    return v


def kunc_v(p, v0, k0, k0p, order=5, min_strain=0.01):
    """
    find volume at given pressure using brenth in scipy.optimize

    :param p: pressure in GPa
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :param order: order of Kunc function
    :param min_strain: defining minimum v/v0 value to search volume for
    :return: unit-cell volume at high pressure in GPa
    :note: a wrapper function vectorizing kunc_v_single
    """
    if isuncertainties([p, v0, k0, k0p]):
        f_u = np.vectorize(uct.wrap(kunc_v_single), excluded=[1, 2, 3, 4, 5])
        return f_u(p, v0, k0, k0p, order=order, min_strain=min_strain)
    else:
        f_v = np.vectorize(kunc_v_single, excluded=[1, 2, 3, 4, 5])
        return f_v(p, v0, k0, k0p, order=order, min_strain=min_strain)


"""
def cal_v_kunc(p, k, order=5):
    calculate volume for Kunc EOS

    Parameters
    ----------
    p = pressure
    k = [v0, k0, k0p]
    p_ref = reference pressure
    order = Kunc order

    Returns
    -------
    volume
    return kunc_v(p, k[0], k[1], k[2], order=order)
"""


def kunc_dPdV(v, v0, k0, k0p, order=5, precision=1.e-5):
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
    def f_scalar(v, v0, k0, k0p, order=order, precision=1.e-5):
        return derivative(kunc_p, v, args=(v0, k0, k0p, order),
                          dx=v0 * precision)
    f_v = np.vectorize(f_scalar, excluded=[1, 2, 3, 4, 5])
    return f_v(v, v0, k0, k0p, order=order, precision=precision)


def kunc_k_num(v, v0, k0, k0p, order=5, precision=1.e-5):
    """
    calculate bulk modulus numerically from volume, not pressure
    according to test this differs from analytical result by 1.e-5

    :param v: unit-cell volume in A^3
    :param v0: unit-cell volume in A^3 at 1 bar
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at reference conditions
    :param precision: precision for numerical calc (default = 1.e-5 * v0)
    :return: bulk modulus
    """
    return -1. * v * kunc_dPdV(v, v0, k0, k0p, order=order,
                               precision=precision)
