from collections import OrderedDict
import uncertainties as uct
from scipy import constants
import periodictable as ptable
from ..conversion import vol_uc2mol
from .objs import MGEOS, JHEOS
from ..conversion import vol_mol2uc

v_ref = 3.9231**3  # default v0
n = 1.
z = 4.
ef_v = 0.001
ef_temp = 0.05
mass = ptable.formula("Pt").mass * 1.e-3
v0_mol = vol_uc2mol(v_ref, z)
rho0 = mass / v0_mol  # kg/m^3


class Fei2004(MGEOS):
    """
    Fei et al. 2004 PEPI 143-144, 515+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(273., 3.0)),
                                 ('k0p', uct.ufloat(4.8, 0.03))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.69, 0.03)),
                                 ('q', uct.ufloat(0.5, 0.5)),
                                 ('theta0', uct.ufloat(230., 0.0))])
        reference = 'Fei et al. 2004 PEPI 143-144, 515+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='constq', reference=reference)


class Fei2007vinet(MGEOS):
    """
    Fei et al. 2007 PNAS 104, 9182+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(277., 0.0)),
                                 ('k0p', uct.ufloat(5.08, 0.02))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.72, 0.03)),
                                 ('q', uct.ufloat(0.5, 0.5)),
                                 ('theta0', uct.ufloat(230., 0.0))])
        reference = 'Fei et al. 2007 PNAS 104, 9182+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='vinet', eqn_th='constq', reference=reference)


class Fei2007bm3(MGEOS):
    """
    Fei et al. 2007 PNAS 104, 9182+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(277., 0.0)),
                                 ('k0p', uct.ufloat(4.95, 0.02))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.72, 0.03)),
                                 ('q', uct.ufloat(0.5, 0.5)),
                                 ('theta0', uct.ufloat(230., 0.0))])
        reference = 'Fei et al. 2007 PNAS 104, 9182+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='constq', reference=reference)


class Dorfman2012(MGEOS):
    """
    Dorfman et al. 2012, JGR 117, B08210
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(277., 0.0)),
                                 ('k0p', uct.ufloat(5.43, 0.02))])
        reference = 'Dorfman et al. 2012, JGR 117, B08210'
        MGEOS.__init__(self, n, z, params_st=params_st,
                       eqn_st='vinet', reference=reference)


class Ye2017(MGEOS):
    """
    Ye et al. 2017. JGR 10.1002/2016JB013811
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(277.3, 0.0)),
                                 ('k0p', uct.ufloat(5.226, 0.033))])
        reference = 'Ye et al. 2017. JGR 10.1002/2016JB013811'
        MGEOS.__init__(self, n, z, params_st=params_st,
                       eqn_st='vinet', reference=reference)


class Yokoo2009(MGEOS):
    """
    Yokoo et al. 2009. PRB 80, 104114
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(276.4, 0.0)),
                                 ('k0p', uct.ufloat(5.12, 0.1))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.63, 0.0)),
                                 ('a', uct.ufloat(0.39, 0.08)),
                                 ('b', uct.ufloat(5.2, 1.1)),
                                 ('theta0', uct.ufloat(230., 0.0))])
        params_el = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('a', uct.ufloat(0.011316, 0.)),
                                 ('b', uct.ufloat(-5.6486e-7, 0.0)),
                                 ('c', uct.ufloat(2.67e-7, 0.0)),
                                 ('d', uct.ufloat(-2.8531e-11, 0.0))])
        reference = 'Yokoo et al. 2009. PRB 80, 104114'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       params_el=params_el,
                       eqn_st='bm3', eqn_th='tange', eqn_el='tsuchiya',
                       reference=reference)


class Jamieson1982(JHEOS):
    """
    Jamieson et al. 1982. High pressure research in geophysics.
    """

    def __init__(self, v0=v_ref):
        mass_shock = mass * 1.e3  # to mass in g
        three_r = 0.12786 / (3. * n * constants.R / (mass_shock)) *\
            3. * constants.R
        rho0 = 21.4449  # g/cm^3
        params_hugoniot = OrderedDict([('rho0', uct.ufloat(rho0, 0.0)),
                                       ('c0', uct.ufloat(3.574, 0.0)),
                                       ('s', uct.ufloat(1.582, 0.0))])
        v0 = vol_mol2uc(mass / (rho0 * 1.e3), z)
        params_therm = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                    ('gamma0', uct.ufloat(2.40, 0.0)),
                                    ('q', uct.ufloat(1.0, 0.0)),
                                    ('theta0', uct.ufloat(200., 0.0))])
        reference = 'Jamieson et al. 1982. High pressure research in \
            geophysics.'
        JHEOS.__init__(self, n, z, mass_shock, params_hugoniot, params_therm,
                       three_r=three_r, reference=reference)


class Holmes1989(MGEOS):
    """
    Holmes et al. 1989. JAP 66, 2962+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(266., 0.0)),
                                 ('k0p', uct.ufloat(5.81, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('alpha0', uct.ufloat(0.261e-4, 0.0)),
                                 ('k0', uct.ufloat(266., 0.0))])
        reference = 'Holmes et al. 1989. JAP 66, 2962+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='vinet', eqn_th='alphakt', reference=reference)


class Dorogokupets2007(MGEOS):
    """
    Dorogokupets and Dewaele. 2007. HPR 27, 431+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(277.3, 0.0)),
                                 ('k0p', uct.ufloat(5.12, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.82, 0.0)),
                                 ('gamma_inf', uct.ufloat(1.83, 0.0)),
                                 ('beta', uct.ufloat(8.11, 0.0)),
                                 ('theta0', uct.ufloat(220., 0.0))])
        params_anh = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                  ('a0', uct.ufloat(-166.9e-6, 0.0)),
                                  ('m', uct.ufloat(4.32, 0.0))])
        params_el = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('e0', uct.ufloat(260.e-6, 0.0)),
                                 ('g', uct.ufloat(2.4, 0.0))])
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
                                 ('k0', uct.ufloat(275.0, 0.0)),
                                 ('k0p', uct.ufloat(5.43, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.77, 0.0)),
                                 ('gamma_inf', uct.ufloat(0.43, 0.0)),
                                 ('beta', uct.ufloat(2.26, 0.0)),
                                 ('theta01', uct.ufloat(184., 0.0)),
                                 ('m1', uct.ufloat(1.5, 0.0)),
                                 ('theta02', uct.ufloat(137., 0.0)),
                                 ('m2', uct.ufloat(1.5, 0.0))])
        params_anh = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                  ('a0', uct.ufloat(0.e-6, 0.0)),
                                  ('m', uct.ufloat(0., 0.0))])
        params_el = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('e0', uct.ufloat(79.e-6, 0.0)),
                                 ('g', uct.ufloat(0.26, 0.0))])
        reference = 'Dorogokupets et al. 2015. RGG 56, 172+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       params_anh=params_anh, params_el=params_el,
                       eqn_st='vinet', eqn_th='dorogokupets2015',
                       eqn_anh='zharkov', eqn_el='zharkov',
                       reference=reference)
