import lmfit
from .eqn_therm_constq import constq_pth
from .eqn_therm_Speziale import speziale_pth
from .eqn_therm_Tange import tange_pth
from .eqn_therm_Dorogokupets2007 import dorogokupets2007_pth
from .eqn_therm_Dorogokupets2015 import dorogokupets2015_pth


class ConstqModel(lmfit.Model):
    """
    lmfit Model class for Constant Q model fitting
    """

    def __init__(self, n, z, independent_vars=['v', 'temp'],
                 param_names=['v0', 'gamma0', 'q', 'theta0'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param n: number of elements in a chemical formula
        :param z: number of formula unit in a unit cell
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3 and temperature in K
        :param param_names: define parameter names, v0, gamma0, q, theta0
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(ConstqModel, self).__init__(constq_pth, n=n, z=z, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('gamma0', min=0.)
        self.set_param_hint('q')
        self.set_param_hint('theta0', min=0.)


class SpezialeModel(lmfit.Model):
    """
    lmfit Model class for Speziale model fitting
    """

    def __init__(self, n, z, independent_vars=['v', 'temp'],
                 param_names=['v0', 'gamma0', 'q0', 'q1', 'theta0'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param n: number of elements in a chemical formula
        :param z: number of formula unit in a unit cell
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3 and temperature in K
        :param param_names: define parameter names, v0, gamma0, q0, q1, theta0
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(SpezialeModel, self).__init__(speziale_pth, n=n, z=z, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('gamma0', min=0.)
        self.set_param_hint('q0')
        self.set_param_hint('q1')
        self.set_param_hint('theta0', min=0.)


class TangeModel(lmfit.Model):
    """
    lmfit Model class for Tange model fitting
    """

    def __init__(self, n, z, independent_vars=['v', 'temp'],
                 param_names=['v0', 'gamma0', 'a', 'b', 'theta0'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param n: number of elements in a chemical formula
        :param z: number of formula unit in a unit cell
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3 and temperature in K
        :param param_names: define parameter names, v0, gamma0, a, b, theta0
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(TangeModel, self).__init__(tange_pth, n=n, z=z, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('gamma0', min=0.)
        self.set_param_hint('a')
        self.set_param_hint('b')
        self.set_param_hint('theta0', min=0.)


class Dorogokupets2007Model(lmfit.Model):
    """
    lmfit Model class for Dorogokupets2007 model fitting
    """

    def __init__(self, n, z, independent_vars=['v', 'temp'],
                 param_names=['v0', 'gamma0', 'gamma_inf', 'beta', 'theta0'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param n: number of elements in a chemical formula
        :param z: number of formula unit in a unit cell
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3 and temperature in K
        :param param_names: define parameter names, v0, gamma0, gamma_inf,
            beta, theta0
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(Dorogokupets2007Model, self).__init__(dorogokupets2007_pth, n=n,
                                                    z=z, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('gamma0', min=0.)
        self.set_param_hint('gamma_inf', min=0.)
        self.set_param_hint('beta')
        self.set_param_hint('theta0', min=0.)


class Dorogokupets2015Model(lmfit.Model):
    """
    lmfit Model class for Dorogokupets2015 model fitting
    """

    def __init__(self, n, z,  independent_vars=['v', 'temp'],
                 param_names=['v0', 'gamma0', 'gamma_inf', 'beta', 'theta01',
                              'm1', 'theta02', 'm2'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param n: number of elements in a chemical formula
        :param z: number of formula unit in a unit cell
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3 and temperature in K
        :param param_names: define parameter names, v0, gamma0, gamma_inf,
            beta, theta01, m1, theta02, m2
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(Dorogokupets2015Model, self).__init__(dorogokupets2015_pth, n=n,
                                                    z=z, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('gamma0', min=0.)
        self.set_param_hint('gamma_inf', min=0.)
        self.set_param_hint('beta')
        self.set_param_hint('theta01', min=0.)
        self.set_param_hint('m1')
        self.set_param_hint('theta02', min=0.)
        self.set_param_hint('m2')
