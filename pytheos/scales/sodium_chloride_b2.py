from collections import OrderedDict
import uncertainties as uct
import periodictable as ptable
import scipy.constants as constants
from ..conversion import vol_uc2mol, vol_mol2uc
from .objs import MGEOS

v_ref_dorogokupets = 40.73302331777446  # 24.53 in cm3/mol
v_ref_fei = 41.35
n = 2.
z = 1.
ef_v = 0.001
ef_temp = 0.05
#mass = ptable.formula("NaCl").mass * 1.e-3  # to kg
#v0_mol = vol_uc2mol(v_ref, z) # 24.53
#rho0 = mass / v0_mol  # kg/m^3

class Dorogokupets2007(MGEOS):
    """
    Dorogokupets and Dewaele. 2007. HPR 27, 431+
    I can reproduce their table.
    However, table for Fei2007 in this paper does not seem to be correct even for 300 K isotherm.
    """

    def __init__(self, v0=v_ref_dorogokupets):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(29.72, 0.0)),
                                 ('k0p', uct.ufloat(5.14, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(1.64, 0.0)),
                                 ('gamma_inf', uct.ufloat(1.23, 0.0)),
                                 ('beta', uct.ufloat(6.83, 0.0)),
                                 ('theta0', uct.ufloat(270., 0.0))])
        params_anh = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                  ('a0', uct.ufloat(-24.e-6, 0.0)),
                                  ('m', uct.ufloat(7.02, 0.0))])
        params_el = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('e0', uct.ufloat(0.e-6, 0.0)),
                                 ('g', uct.ufloat(0.0, 0.0))])
        reference = 'Dorogokupets and Dewaele. 2007. HPR 27, 431+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       params_anh=params_anh, params_el=params_el,
                       eqn_st='vinet', eqn_th='dorogokupets2007',
                       eqn_anh='zharkov', eqn_el='zharkov',
                       reference=reference)

class Fei2007vinet(MGEOS):
    """
    Fei et al. 2007 PNAS 104, 9182+
    I can reproduce their figure.
    """

    def __init__(self, v0=v_ref_fei):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(26.86, 2.9)),
                                 ('k0p', uct.ufloat(5.25, 0.26))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(1.7, 0.0)),
                                 ('q', uct.ufloat(0.5, 0.3)),
                                 ('theta0', uct.ufloat(290., 0.0))])
        reference = 'Fei et al. 2007 PNAS 104, 9182+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='vinet', eqn_th='constq', reference=reference)


class Fei2007bm3(MGEOS):
    """
    Fei et al. 2007 PNAS 104, 9182+
    I can reproduce their figure.
    """

    def __init__(self, v0=v_ref_fei):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(30.69, 2.9)),
                                 ('k0p', uct.ufloat(4.33, 0.26))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(1.7, 0.0)),
                                 ('q', uct.ufloat(0.5, 0.3)),
                                 ('theta0', uct.ufloat(290., 0.0))])
        reference = 'Fei et al. 2007 PNAS 104, 9182+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='constq', reference=reference)