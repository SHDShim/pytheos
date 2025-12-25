from collections import OrderedDict
import uncertainties as uct
import periodictable as ptable
import scipy.constants as constants
from ..conversion import vol_uc2mol, vol_mol2uc
from .objs import MGEOS

v_ref = 88.967  # 
n = 1.
z = 4.
ef_v = 0.001
ef_temp = 0.05
mass = ptable.formula("NaCl").mass * 1.e-3  # to kg
v0_mol = vol_uc2mol(v_ref, z) # 
rho0 = mass / v0_mol  # kg/m^3

class Fei2007vinet(MGEOS):
    """
    Fei et al. 2007 PNAS 104, 9182+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(1.16, 0.14)),
                                 ('k0p', uct.ufloat(8.23, 0.31))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.05, 0.0)),
                                 ('q', uct.ufloat(0.6, 0.3)),
                                 ('theta0', uct.ufloat(75.1, 0.0))])
        reference = 'Fei et al. 2007 PNAS 104, 9182+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='vinet', eqn_th='constq', reference=reference)


class Fei2007bm3(MGEOS):
    """
    Fei et al. 2007 PNAS 104, 9182+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(1.43, 0.14)),
                                 ('k0p', uct.ufloat(8.02, 0.31))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.05, 0.0)),
                                 ('q', uct.ufloat(0.6, 0.3)),
                                 ('theta0', uct.ufloat(75.1, 0.0))])
        reference = 'Fei et al. 2007 PNAS 104, 9182+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='constq', reference=reference)