from collections import OrderedDict
import uncertainties as uct
import periodictable as ptable
import scipy.constants as constants
from ..conversion import vol_uc2mol, vol_mol2uc
from .objs import MGEOS, JHEOS

v_ref = 74.698  # 3.9231**3  # default v0
n = 2.
z = 4.
ef_v = 0.001
ef_temp = 0.05
mass = ptable.formula("MgO").mass * 1.e-3  # to kg
v0_mol = vol_uc2mol(v_ref, z)
rho0 = mass / v0_mol  # kg/m^3


class Jamieson1982(JHEOS):
    """
    Jamieson et al. 1982. High pressure research in geophysics.
    """

    def __init__(self, v0=v_ref):
        mass_shock = mass * 1.e3  # to mass in g
        three_r = 1.23754 / (3. * n * constants.R / mass_shock) * 3. *\
            constants.R
        rho0 = 3.585  # g/cm^3 from Jamieson
        params_hugoniot = OrderedDict([('rho0', uct.ufloat(rho0, 0.0)),
                                       ('c0', uct.ufloat(6.597, 0.0)),
                                       ('s', uct.ufloat(1.369, 0.0))])
        v0 = vol_mol2uc(mass / (rho0 * 1.e3), z)
        params_therm = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                    ('gamma0', uct.ufloat(1.32, 0.0)),
                                    ('q', uct.ufloat(1.0, 0.0)),
                                    ('theta0', uct.ufloat(760., 0.0))])
        reference = 'Jamieson et al. 1982. High pressure research in \
            geophysics.'
        JHEOS.__init__(self, n, z, mass_shock, params_hugoniot, params_therm,
                       three_r=three_r, reference=reference)


class Zha2000(MGEOS):
    """
    Zha et al. 2000. PNAS 97, 13494+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(160.2, 0.0)),
                                 ('k0p', uct.ufloat(4.03, 0.0))])
        reference = 'Zha et al. 2000. PNAS 97, 13494+'
        MGEOS.__init__(self, n, z, params_st=params_st,
                       eqn_st='bm3', reference=reference)


class Ye2017(MGEOS):
    """
    Ye et al. 2017. JGR 10.1002/2016JB013811
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(160.3, 0.0)),
                                 ('k0p', uct.ufloat(4.109, 0.022))])
        reference = 'Ye et al. 2017. JGR 10.1002/2016JB013811'
        MGEOS.__init__(self, n, z, params_st=params_st,
                       eqn_st='vinet', reference=reference)


class Speziale2001(MGEOS):
    """
    Speziale et al. 2001. JGR 106, 515+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(160.2, 0.0)),
                                 ('k0p', uct.ufloat(3.99, 0.01))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(1.524, 0.03)),
                                 ('q0', uct.ufloat(1.65, 0.4)),
                                 ('q1', uct.ufloat(11.8, 0.2)),
                                 ('theta0', uct.ufloat(773., 0.0))])
        reference = 'Speziale et al. 2001. JGR 106, 515+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='speziale', reference=reference)


class Tange2009(MGEOS):
    """
    Tange et al. 2009. JGR 114, B03208+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(160.63, 0.18)),
                                 ('k0p', uct.ufloat(4.367, 0.013))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(1.442, 0.015)),
                                 ('a', uct.ufloat(0.138, 0.019)),
                                 ('b', uct.ufloat(5.4, 1.1)),
                                 ('theta0', uct.ufloat(761., 13.0))])
        reference = 'Tange, 2009, ?'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='vinet', eqn_th='tange', reference=reference)


class Dorogokupets2007(MGEOS):
    """
    Dorogokupets and Dewaele. 2007. HPR 27, 431+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(160.3, 0.0)),
                                 ('k0p', uct.ufloat(4.18, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(1.50, 0.0)),
                                 ('gamma_inf', uct.ufloat(0.75, 0.0)),
                                 ('beta', uct.ufloat(2.96, 0.0)),
                                 ('theta0', uct.ufloat(760., 0.0))])
        params_anh = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                  ('a0', uct.ufloat(-14.9e-6, 0.0)),
                                  ('m', uct.ufloat(5.12, 0.0))])
        params_el = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('e0', uct.ufloat(0.e-6, 0.0)),
                                 ('g', uct.ufloat(0.0, 0.0))])
        reference = 'Dorogokupets and Dewaele. 2007. HPR 27, 431+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       params_anh=params_anh, params_el=params_el,
                       eqn_st='vinet', eqn_th='dorogokupets2007',
                       eqn_anh='zharkov', eqn_el='zharkov',
                       reference=reference)


class Dorogokupets2015(MGEOS):
    """
    Dorogokupets et al. 2015. RGG 56, 172+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(160.3, 0.0)),
                                 ('k0p', uct.ufloat(4.25, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(1.53, 0.0)),
                                 ('gamma_inf', uct.ufloat(0.624, 0.0)),
                                 ('beta', uct.ufloat(2.115, 0.0)),
                                 ('theta01', uct.ufloat(747., 0.0)),
                                 ('m1', uct.ufloat(1.5, 0.0)),
                                 ('theta02', uct.ufloat(399., 0.0)),
                                 ('m2', uct.ufloat(1.5, 0.0))])
        params_anh = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                  ('a0', uct.ufloat(-15.9e-6, 0.0)),
                                  ('m', uct.ufloat(4.48, 0.0))])
        params_el = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('e0', uct.ufloat(0.e-6, 0.0)),
                                 ('g', uct.ufloat(0.0, 0.0))])
        reference = 'Dorogokupets et al. 2015. RGG 56, 172+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       params_anh=params_anh, params_el=params_el,
                       eqn_st='kunc', eqn_th='dorogokupets2015',
                       eqn_anh='zharkov', eqn_el='zharkov',
                       reference=reference)
