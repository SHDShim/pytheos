"""
PythEOS
=======

PythEOS is an open source python toolbox for equation of state (EOS).
It includes:
- a range of different equations
- EOS fitting tools
- a range of pressure scales

Some important motivations are:
- reproducible error propagation for EOS
- flexible EOS fitting with a wide range of combinations of static, thermal, electronic, and anharmonic effects
- reliable conversion of pressure scale

**Notations**

:param n: number of elements in a chemical formula
:param z: number of formula unit in a unit cell
:param m: molar mass in gram

:param p: pressure in GPa
:param v: unit-cell volume in A^3
:param temp: temperature in K
:param rho: density in g/cm^3

:param k0: bulk modulus at reference conditions
:param k0p: pressure derivative of bulk modulus at reference conditions

:param c0: velocity at 1 bar in km/s
:param s: slope of the velocity change
:param rho0: density at 1 bar in g/cm^3
:param v0: unit-cell volume in A^3 at 1 bar
:param e0: parameter in K-1 for the Zharkov equation
:param g: parameter for the Zharkov equation
:param gamma0: Gruneisen parameter at 1 bar
:param q: logarithmic derivative of Gruneisen parameter
:param theta0: Debye temperature in K

:param t_ref: reference temperature, 300 K
:param three_r: 3 times gas constant.
    Jamieson modified this value to compensate for mismatches
:param c_v: heat capacity

:param min_strain: defining minimum v/v0 value to search volume for

:param params_hugoniot: hugoniot parameters.
    [rho0 in g/cm3, c0 in km/s, s] for linear case and
    [rho0 in g/cm3, a in km/s, b, c in s/km] for non-linear case.
    See Jamieson for detail.
:param params_therm: thermal parameters for the constq equations
    [v0 in A3, gamma0, q, theta0 in K]

:param c_v: heat capacity for hugoniot temperature calculation
:param nonlinear: nonlinear Us-Up fit particularly for Jamieson's gold
    scale
:param reference: reference for the EOS

"""
from .eqn_bm3 import bm3_p, bm3_v, bm3_k, bm3_g, \
    bm3_small_f, bm3_big_F, bm3_k_num  # , bm3_p_u, bm3_v_u, bm3_k_u
from .eqn_vinet import vinet_p, vinet_v, vinet_k, \
    vinet_k_num  # ,  vinet_p_u, vinet_v_u, vinet_k_u
from .eqn_kunc import kunc_p, kunc_v, kunc_k_num  # , kunc_p_u, kunc_v_u
from .eqn_hugoniot import hugoniot_p, hugoniot_t, hugoniot_rho
from .eqn_jamieson import jamieson_pst, jamieson_pth
from .eqn_debye import debye_E
from .eqn_therm_constq import constq_grun, \
    constq_debyetemp, constq_pth
from .eqn_therm_Tange import tange_grun, tange_debyetemp, \
    tange_pth
from .eqn_therm_Speziale import speziale_grun, speziale_debyetemp, \
    speziale_pth
from .eqn_therm_Dorogokupets2007 import altshuler_grun, altshuler_debyetemp,\
    dorogokupets2007_pth
from .eqn_therm_Dorogokupets2015 import dorogokupets2015_pth
from .eqn_anharmonic import zharkov_panh
from .eqn_electronic import zharkov_pel
from .fit_static import BM3Model, VinetModel, KuncModel
from .fit_thermal import ConstqModel, Dorogokupets2007Model, \
    Dorogokupets2015Model, SpezialeModel, TangeModel
from .fit_electronic import ZharkovElecModel
from .fit_anharmonic import ZharkovAnhModel
from .conversion import vol_uc2mol
from . import plot
from .scales import gold
from .scales import platinum
from .scales import periclase
