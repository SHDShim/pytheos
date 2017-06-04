"""
2017/05/03 No need to put isuncertainties because there is no particularly
function to be problematic to the functions here.
"""
import numpy as np
import uncertainties as uct
from scipy.optimize import brenth
from scipy.misc import derivative
from .etc import isuncertainties


def bm3_p(v, v0, k0, k0p, p_ref=0.0):
    """
    calculate pressure from 3rd order Birch-Murnathan equation

    :param v: volume at different pressures
    :param v0: volume at reference conditions
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at different conditions
    :param p_ref: reference pressure (default = 0)
    :return: pressure
    """
    return cal_p_bm3(v, [v0, k0, k0p], p_ref=p_ref)


def cal_p_bm3(v, k, p_ref=0.0):
    """
    calculate pressure from 3rd order Birch-Murnaghan equation

    :param v: volume at different pressures
    :param k: [v0, k0, k0p]
    :param p_ref: reference pressure, default = 0.
    :return: static pressure
    """
    vvr = v / k[0]
    p = (p_ref - 0.5 * (3. * k[1] - 5. * p_ref) * (1. - vvr**(-2. / 3.)) +
         9. / 8. * k[1] * (k[2] - 4. + 35. / 9. * p_ref / k[1]) *
         (1. - vvr**(-2. / 3.))**2.) * vvr**(-5. / 3.)
    return p


def bm3_v_single(p, v0, k0, k0p, p_ref=0.0, min_strain=0.01):
    """
    find volume at given pressure using brenth in scipy.optimize
    this is for single p value, not vectorized
    this cannot handle uncertainties

    :param p: pressure
    :param v0: volume at reference conditions
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at different conditions
    :param p_ref: reference pressure (default = 0)
    :param min_strain: minimum strain value to find solution (default = 0.01)
    :return: volume at high pressure
    """
    if p <= 1.e-5:
        return v0

    def f_diff(v, v0, k0, k0p, p, p_ref=0.0):
        return bm3_p(v, v0, k0, k0p, p_ref=p_ref) - p
    v = brenth(f_diff, v0, v0 * min_strain, args=(v0, k0, k0p, p, p_ref))
    return v


def bm3_v(p, v0, k0, k0p, p_ref=0.0, min_strain=0.01):
    """
    find volume at given pressure using brenth in scipy.optimize

    :param p: pressure
    :param v0: volume at reference conditions
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at different conditions
    :param p_ref: reference pressure (default = 0)
    :param min_strain: minimum strain value to find solution (default = 0.01)
    :return: volume at high pressure
    """
    if isuncertainties([p, v0, k0, k0p]):
        f_u = np.vectorize(uct.wrap(bm3_v_single), excluded=[1, 2, 3, 4, 5])
        return f_u(p, v0, k0, k0p, p_ref=p_ref, min_strain=min_strain)
    else:
        f_v = np.vectorize(bm3_v_single, excluded=[1, 2, 3, 4, 5])
        return f_v(p, v0, k0, k0p, p_ref=p_ref, min_strain=min_strain)


def cal_v_bm3(p, k):
    """
    calculate volume from bm3 equation. wrapper for bm3_v

    :param p: pressure
    :param k: [v0, k0, k0p]
    :return: volume at high pressure
    """
    return bm3_v(p, k[0], k[1], k[2])


def bm3_k(p, v0, k0, k0p):
    """
    calculate bulk modulus, wrapper for cal_k_bm3
    cannot handle uncertainties

    :param p: pressure
    :param v0: volume at reference conditions
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at different conditions
    :return: bulk modulus at high pressure
    """
    return cal_k_bm3(p, [v0, k0, k0p])


def bm3_dPdV(v, v0, k0, k0p, precision=1.e-5):
    """
    calculate dP/dV for numerical calculation of bulk modulus
    according to test this differs from analytical result by 1.e-5

    :param v: volume
    :param v0: volume at reference conditions
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at different conditions
    :param precision: precision for numerical calculation (default = 1.e-5*v0)
    :return: dP/dV
    """
    def f_scalar(v, v0, k0, k0p, precision=precision):
        return derivative(bm3_p, v, args=(v0, k0, k0p), dx=v0 * precision)
    f_v = np.vectorize(f_scalar, excluded=[1, 2, 3, 4])
    return f_v(v, v0, k0, k0p, precision=precision)


def bm3_k_num(v, v0, k0, k0p, precision=1.e-5):
    """
    calculate bulk modulus numerically from volume, not pressure
    according to test this differs from analytical result by 1.e-5

    :param v: volume
    :param v0: volume at reference conditions
    :param k0: bulk modulus at reference conditions
    :param k0p: pressure derivative of bulk modulus at different conditions
    :param precision: precision for numerical calculation (default = 1.e-5*v0)
    :return: dP/dV
    """
    return -1. * v * bm3_dPdV(v, v0, k0, k0p, precision=precision)


def cal_k_bm3(p, k):
    """
    calculate bulk modulus

    :param p: pressure
    :param k: [v0, k0, k0p]
    :return: bulk modulus at high pressure
    """
    v = cal_v_bm3(p, k)
    return cal_k_bm3_from_v(v, k)


def cal_k_bm3_from_v(v, k):
    v0 = k[0]
    k0 = k[1]
    k0p = k[2]

    f = 0.5 * (np.power(v0 / v, 2. / 3.) - 1.)
    a4 = 9. / 2. * v0 * k0 * (k0p - 4.)
    bulk_modulus = 1. / 3. / v0 * np.power((1. + 2. * f), (5. / 2.)) * \
        (5. * f * (3. * v0 * k0 + a4 * f) +
         (1. + 2. * f) * (3. * v0 * k0 + 2. * a4 * f))
    return bulk_modulus


def bm3_g(p, v0, g0, g0p, k0, k0p):
    """
    calculate shear modulus at given pressure.
    not fully tested with mdaap.

    :param p: pressure
    :param v0: volume at reference condition
    :param g0: shear modulus at reference condition
    :param g0p: pressure derivative of shear modulus at reference condition
    :param k0: bulk modulus at reference condition
    :param k0p: pressure derivative of bulk modulus at reference condition
    :return: shear modulus at high pressure
    """
    return cal_g_bm3(p, [g0, g0p], [v0, k0, k0p])


def cal_g_bm3(p, g, k):
    """
    calculate shear modulus at given pressure

    :param p: pressure
    :param g: [g0, g0p]
    :param k: [v0, k0, k0p]
    :return: shear modulus at high pressure
    """
    v = cal_v_bm3(p, k)
    v0 = k[0]
    k0 = k[1]
    kp = k[2]
    g0 = g[0]
    gp = g[1]
    f = 0.5 * ((v / v0)**(-2. / 3.) - 1.)
    return (1. + 2. * f)**(5. / 2.) * (g0 + (3. * k0 * gp - 5. * g0) * f +
                                       (6. * k0 * gp - 24. * k0 - 14. * g0 +
                                        9. / 2. * k0 * kp) * f**2.)


def bm3_small_f(v, v0):
    """
    calculate finite strain little f for linearized form
    not fully tested

    :param v: volume
    :param v0: volume at reference conditions
    :return: small f
    """
    return cal_small_f(v, v0)


def bm3_big_F(p, v, v0):
    """
    calculate big F for linearlized form
    not fully tested

    :param p:
    :param f:
    :return:
    """
    f = bm3_small_f(v, v0)
    return cal_big_F(p, f)


def cal_small_f(v, v0):
    """
    calculate finite strain little f for linearized form
    not fully tested

    :param v: volume
    :param v0: volume at reference condition
    :return: small f
    """
    return 0.5 * (np.power(v0 / v, 2. / 3.) - 1.)


def cal_big_F(p, f):
    """
    calculate finite strain big F for linearized form
    not fully tested

    :param p: pressure
    :param f: small f
    :return: big F
    """
    return p / (3. * f * np.power((1. + 2. * f), 2.5))
