from collections import OrderedDict
import uncertainties as uct
import periodictable as ptable
import scipy.constants as constants
from ..conversion import vol_uc2mol, vol_mol2uc
from .objs import MGEOS

v_ref = 179.44  # 27.015 in cm3/mol
n = 2.
z = 4.
ef_v = 0.001
ef_temp = 0.05
mass = ptable.formula("NaCl").mass * 1.e-3  # to kg
v0_mol = vol_uc2mol(v_ref, z)
rho0 = mass / v0_mol  # kg/m^3

class Dorogokupets2007(MGEOS):
    """
    Dorogokupets and Dewaele. 2007. HPR 27, 431+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(23.83, 0.0)),
                                 ('k0p', uct.ufloat(5.09, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(1.64, 0.0)),
                                 ('gamma_inf', uct.ufloat(1.12, 0.0)),
                                 ('beta', uct.ufloat(4.36, 0.0)),
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
