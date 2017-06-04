
from scipy import constants


def alphakt_pth(v, temp, v0, alpha0, k0, n, z, t_ref=300.,
                three_r=3. * constants.R):
    """
    calculate thermal pressure from thermal expansion and bulk modulus

    :param v: unit-cell volume in A^3
    :param temp: temperature in K
    :param v0: unit-cell volume in A^3 at 1 bar
    :param alpha0: thermal expansion parameter at 1 bar in K-1
    :param k0: bulk modulus in GPa
    :param n: number of atoms in a formula unit
    :param z: number of formula unit in a unit cell
    :param t_ref: reference temperature
    :param three_r: 3R in case adjustment is needed
    :return: thermal pressure in GPa
    """
    return alpha0 * k0 * (temp - t_ref)
