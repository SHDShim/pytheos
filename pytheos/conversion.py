import scipy.constants as constants
import numpy as np


def velocities_to_moduli(rho, v_phi, v_s):
    """
    convert velocities to moduli
    mainly to support Burnman operations

    :param rho: density in kg/m^3
    :param v_phi: bulk sound speed in m/s
    :param v_s: shear velocity in m/s
    :return: K_s and G
    """
    return v_phi * v_phi * rho, v_s * v_s * rho


def moduli_to_velocities(rho, K_s, G):
    """
    convert moduli to velocities
    mainly to support Burnman operations

    :param rho: density in kg/m^3
    :param v_phi: adiabatic bulk modulus in Pa
    :param v_s: shear modulus in Pa
    :return: bulk sound speed and shear velocity
    """
    return np.sqrt(K_s / rho), np.sqrt(G / rho)


def vol_uc2mol(v_uc, z_uc):
    """
    convert unit cell volume in A^3 to molar volume in m^3/mol

    :param v_uc: unit-cell volume in A^3
    :param z_uc: number of formula unit in a unit cell (z)
    :return: molar volume
    """
    return v_uc * 1.e-30 * constants.N_A / z_uc


def vol_mol2uc(v_mol, z_uc):
    """
    convert unit cell volume in A^3 to molar volume in m^3/mol

    :param v_mol: molar volume in m^3/mol
    :param z_uc: number of formula unit in a unit cell (z)
    :return: unit-cell volume in A^3
    """
    return v_mol * 1.e30 / constants.N_A * z_uc
