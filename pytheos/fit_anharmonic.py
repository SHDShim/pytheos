import lmfit
from .eqn_anharmonic import zharkov_panh


class ZharkovAnhModel(lmfit.Model):
    """
    lmfit Model class for Zharkov anharmonic fitting
    """

    def __init__(self, n, z, independent_vars=['v', 'temp'],
                 param_names=['v0', 'a0', 'm'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param n: number of elements in a chemical formula
        :param z: number of formula unit in a unit cell
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3 and temperature in K
        :param param_names: define parameter names, v0, a0, m
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(ZharkovAnhModel, self).__init__(zharkov_panh, n=n, z=z, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('a0')
        self.set_param_hint('m')
