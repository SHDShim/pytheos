from scipy import constants
from .conversion import vol_uc2mol
import numpy as np


def zharkov_pel(v, temp, v0, e0, g, n, z, t_ref=300.,
                three_r=3. * constants.R):
    """
    calculate electronic contributions in pressure for the Zharkov equation
    the equation can be found in Sokolova and Dorogokupets 2013

    :param v: unit-cell volume in A^3
    :param temp: temperature in K
    :param v0: unit-cell volume in A^3 at 1 bar
    :param e0: parameter in K-1 for the Zharkov equation
    :param g: parameter for the Zharkov equation
    :param n: number of atoms in a formula unit
    :param z: number of formula unit in a unit cell
    :param t_ref: reference temperature, 300 K
    :param three_r: 3 times gas constant
    :return: electronic contribution in GPa
    """
    v_mol = vol_uc2mol(v, z)
    x = v / v0
#    a = a0 * np.power(x, m)

    def f(t):
        return three_r * n / 2. * e0 * np.power(x, g) * np.power(t, 2.) * \
            g / v_mol * 1.e-9
    return f(temp) - f(t_ref)


def tsuchiya_pel(v, temp, v0, a, b, c, d, n, z, three_r=3. * constants.R,
                 t_ref=300.):
    """
    calculate electronic contributions in pressure for the Tsuchiya equation

    :param v: unit-cell volume in A^3
    :param temp: temperature in K
    :param v0: unit-cell volume in A^3 at 1 bar
    :param a: parameter for the Tsuchiya equation
    :param b: parameter for the Tsuchiya equation
    :param c: parameter for the Tsuchiya equation
    :param d: parameter for the Tsuchiya equation
    :param n: number of atoms in a formula unit
    :param z: number of formula unit in a unit cell
    :param t_ref: reference temperature, 300 K
    :param three_r: 3 times gas constant
    :return: electronic contribution in GPa
    :note: n, z, three_r are not used but in there for consistency
        with other electronic contribution equations
    """
    def f(temp):
        return a + b * temp + c * np.power(temp, 2.) + d * np.power(temp, 3.)
    return f(temp) - f(t_ref)
