import numpy as np
from scipy import constants
from .conversion import vol_uc2mol


def zharkov_panh(v, temp, v0, a0, m, n, z, t_ref=300.,
                 three_r=3. * constants.R):
    """
    calculate pressure from anharmonicity for Zharkov equation
    the equation is from Dorogokupets 2015

    :param v: unit-cell volume in A^3
    :param temp: temperature in K
    :param v0: unit-cell volume in A^3 at 1 bar
    :param a0: parameter in K-1 for the Zharkov equation
    :param m: parameter for the Zharkov equation
    :param n: number of elements in a chemical formula
    :param z: number of formula unit in a unit cell
    :param three_r: 3 times gas constant
    :return: anharmonic contribution for pressure in GPa
    """
    v_mol = vol_uc2mol(v, z)
    x = v / v0
    a = a0 * np.power(x, m)

    def f(t):
        return three_r * n / 2. * a * m / v_mol * np.power(t, 2.) * 1.e-9

    return f(temp) - f(t_ref)
