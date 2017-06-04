from collections import OrderedDict
import uncertainties as uct
from scipy import constants
import periodictable as ptable
from ..conversion import vol_uc2mol, vol_mol2uc
from .objs import MGEOS, JHEOS

v_ref = 4.07860 ** 3  # default v0
n = 1.
z = 4.
ef_v = 0.001
ef_temp = 0.05
mass = ptable.formula("Au").mass * 1.e-3
v0_mol = vol_uc2mol(v_ref, z)
rho0 = mass / v0_mol  # kg/m^3


class Jamieson1982L(JHEOS):
    """
    Jamieson et al. 1982. High pressure research in geophysics.
    Fit C in table 2.
    """

    def __init__(self, v0=v_ref):
        mass_shock = mass * 1.e3  # to mass in g
        three_r = 0.12500 / (3. * n * constants.R / mass_shock) *\
            3. * constants.R
        rho0 = 19.2827
        params_hugoniot = OrderedDict([('rho0', uct.ufloat(rho0, 0.0)),
                                       ('a', uct.ufloat(2.975, 0.0)),
                                       ('b', uct.ufloat(1.896, 0.0)),
                                       ('c', uct.ufloat(-0.309, 0.0))])
        v0 = vol_mol2uc(mass / (rho0 * 1.e3), z)
        params_therm = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                    ('gamma0', uct.ufloat(3.215, 0.0)),
                                    ('q', uct.ufloat(1.0, 0.0)),
                                    ('theta0', uct.ufloat(170., 0.0))])
        reference = 'Jamieson et al. 1982. High pressure research in \
                    geophysics. Fit C in table 2.'
        JHEOS.__init__(self, n, z, mass_shock, params_hugoniot, params_therm,
                       three_r=three_r, reference=reference, nonlinear=True)


class Jamieson1982H(JHEOS):
    """
    Jamieson et al. 1982. High pressure research in geophysics.
    Fit A in table 2.
    """

    def __init__(self, v0=v_ref):
        mass_shock = mass * 1.e3  # to mass in g
        three_r = 0.12500 / (3. * n * constants.R / mass_shock) *\
            3. * constants.R
        rho0 = 19.2827
        params_hugoniot = OrderedDict([('rho0', uct.ufloat(rho0, 0.0)),
                                       ('c0', uct.ufloat(3.071, 0.0)),
                                       ('s', uct.ufloat(1.536, 0.0))])
        v0 = vol_mol2uc(mass / (rho0 * 1.e3), z)
        params_therm = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                    ('gamma0', uct.ufloat(3.215, 0.0)),
                                    ('q', uct.ufloat(1.0, 0.0)),
                                    ('theta0', uct.ufloat(170., 0.0))])
        reference = 'Jamieson et al. 1982. High pressure research in \
                    geophysics. Fit A in table 2.'
        JHEOS.__init__(self, n, z, mass_shock, params_hugoniot, params_therm,
                       three_r=three_r, reference=reference)


class Heinz1984(MGEOS):
    """
    Heinz and Jeanloz. 1984. JAP 55, 885+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(166.65, 5.0)),
                                 ('k0p', uct.ufloat(5.4833, 0.43))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.95, 0.43)),
                                 ('q', uct.ufloat(1.7, 0.7)),
                                 ('theta0', uct.ufloat(170., 0.0))])
        reference = 'Heinz and Jeanloz. 1984. JAP 55, 885+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='constq', reference=reference)


class Tsuchiya2003(MGEOS):
    """
    Tsuchiya 2003 JGR 108. 2462+

    Original k0 = 166.7, k0p = 6.12.  However, I cannot reproduce their table
    exactly unless I change the values to 166.1.

    If reproduce_table=True, the adjusted k0 will be used for reproducing
    the table value down to the first number after decimal point.
    """

    def __init__(self, v0=v_ref, reproduce_table=False):
        if reproduce_table:
            k0 = 166.1
            k0p = 6.12
        else:
            k0 = 166.7
            k0p = 6.12
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(k0, 0.0)),
                                 ('k0p', uct.ufloat(k0p, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(3.16, 0.0)),
                                 ('q', uct.ufloat(2.15, 0.0)),
                                 ('theta0', uct.ufloat(180., 0.0))])
        reference = 'Tsuchiya 2003 JGR 108. 2462+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='vinet', eqn_th='constq', reference=reference)


class Fei2007vinet(MGEOS):
    """
    Fei et al. 2007 PNAS 104, 9182+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(167., 0.0)),
                                 ('k0p', uct.ufloat(6.00, 0.02))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.97, 0.03)),
                                 ('q', uct.ufloat(0.6, 0.3)),
                                 ('theta0', uct.ufloat(170., 0.0))])
        reference = 'Fei et al. 2007 PNAS 104, 9182+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='vinet', eqn_th='constq', reference=reference)


class Fei2007bm3(MGEOS):
    """
    Fei et al. 2007 PNAS 104, 9182+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(167., 0.0)),
                                 ('k0p', uct.ufloat(5.77, 0.02))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.97, 0.03)),
                                 ('q', uct.ufloat(0.6, 0.3)),
                                 ('theta0', uct.ufloat(170., 0.0))])
        reference = 'Fei et al. 2007 PNAS 104, 9182+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='constq', reference=reference)


class Fei2004(MGEOS):
    """
    Fei et al. 2004 PEPI 143-144, 515+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(167., 3.0)),
                                 ('k0p', uct.ufloat(5.00, 0.2))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.97, 0.03)),
                                 ('q', uct.ufloat(0.7, 0.3)),
                                 ('theta0', uct.ufloat(170., 0.0))])
        reference = 'Fei et al. 2004 PEPI 143-144, 515+'
        three_r = 0.125 / 0.12664 * 3. * constants.R

        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='constq', reference=reference,
                       three_r=three_r)


class Shim2002(MGEOS):
    """
    Shim et al 2002 EPSL 203, 729+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(167., 3.0)),
                                 ('k0p', uct.ufloat(5.00, 0.2))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.97, 0.05)),
                                 ('q', uct.ufloat(1.0, 0.1)),
                                 ('theta0', uct.ufloat(170., 0.0))])
        reference = 'Shim et al 2002 EPSL 203, 729+'
        three_r = 0.125 / 0.12664 * 3. * constants.R
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       eqn_st='bm3', eqn_th='constq', reference=reference,
                       three_r=three_r)


class Dorfman2012(MGEOS):
    """
    Dorfman et al. 2012, JGR 117, B08210
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(167., 0.0)),
                                 ('k0p', uct.ufloat(5.88, 0.02))])
        reference = 'Dorfman et al. 2012, JGR 117, B08210'
        MGEOS.__init__(self, n, z, params_st=params_st,
                       eqn_st='vinet', reference=reference)


class Ye2017(MGEOS):
    """
    Ye et al. 2017. JGR 10.1002/2016JB013811
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(167., 0.0)),
                                 ('k0p', uct.ufloat(5.897, 0.022))])
        reference = 'Ye et al. 2017. JGR 10.1002/2016JB013811'
        MGEOS.__init__(self, n, z, params_st=params_st,
                       eqn_st='vinet', reference=reference)


class Yokoo2009(MGEOS):
    """
    Yokoo et al. 2009. PRB 80, 104114.
    In original publication k0p=5.79, but the table fits the best when
    I use k0p=5.749.
    """

    def __init__(self, v0=v_ref, reproduce_table=False):
        if reproduce_table:
            k0p = 5.749
        else:
            k0p = 5.79
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(167.5, 0.0)),
                                 ('k0p', uct.ufloat(k0p, 0.1))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.96, 0.0)),
                                 ('a', uct.ufloat(0.45, 0.09)),
                                 ('b', uct.ufloat(4.2, 0.6)),
                                 ('theta0', uct.ufloat(170., 0.0))])
        params_el = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('a', uct.ufloat(-.00021606, 0.)),
                                 ('b', uct.ufloat(-4.3795e-6, 0.0)),
                                 ('c', uct.ufloat(1.4526e-8, 0.0)),
                                 ('d', uct.ufloat(7.8072e-14, 0.0))])
        reference = 'Yokoo et al. 2009. PRB 80, 104114'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       params_el=params_el,
                       eqn_st='bm3', eqn_th='tange', eqn_el='tsuchiya',
                       reference=reference)


class Dorogokupets2007(MGEOS):
    """
    Dorogokupets and Dewaele. 2007. HPR 27, 431+
    """

    def __init__(self, v0=v_ref):
        params_st = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('k0', uct.ufloat(167.0, 0.0)),
                                 ('k0p', uct.ufloat(5.90, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.89, 0.0)),
                                 ('gamma_inf', uct.ufloat(1.54, 0.0)),
                                 ('beta', uct.ufloat(4.36, 0.0)),
                                 ('theta0', uct.ufloat(170., 0.0))])
        params_anh = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                  ('a0', uct.ufloat(0.e-6, 0.0)),
                                  ('m', uct.ufloat(0.0, 0.0))])
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
                                 ('k0', uct.ufloat(167.0, 0.0)),
                                 ('k0p', uct.ufloat(5.9, 0.0))])
        params_th = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('gamma0', uct.ufloat(2.918, 0.0)),
                                 ('gamma_inf', uct.ufloat(0.66, 0.0)),
                                 ('beta', uct.ufloat(2.406, 0.0)),
                                 ('theta01', uct.ufloat(178., 0.0)),
                                 ('m1', uct.ufloat(1.5, 0.0)),
                                 ('theta02', uct.ufloat(84., 0.0)),
                                 ('m2', uct.ufloat(1.5, 0.0))])
        params_anh = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                  ('a0', uct.ufloat(0.e-6, 0.0)),
                                  ('m', uct.ufloat(0., 0.0))])
        params_el = OrderedDict([('v0', uct.ufloat(v0, ef_v)),
                                 ('e0', uct.ufloat(6.1e-6, 0.0)),
                                 ('g', uct.ufloat(0.66, 0.0))])
        reference = 'Dorogokupets et al. 2015. RGG 56, 172+'
        MGEOS.__init__(self, n, z, params_st=params_st, params_th=params_th,
                       params_anh=params_anh, params_el=params_el,
                       eqn_st='vinet', eqn_th='dorogokupets2015',
                       eqn_anh='zharkov', eqn_el='zharkov',
                       reference=reference)
