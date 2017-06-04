import lmfit
from .eqn_bm3 import bm3_p
from .eqn_vinet import vinet_p
from .eqn_kunc import kunc_p


class BM3Model(lmfit.Model):
    """
    lmfit Model class for BM3 fitting
    """

    def __init__(self, independent_vars=['v'], param_names=['v0', 'k0', 'k0p'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3
        :param param_names: define parameter names, v0, k0, k0p
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(BM3Model, self).__init__(bm3_p, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('k0', min=0.)
        self.set_param_hint('k0p', min=0.)

    # not supported
    # def guess(self, data, x=None, negative=False, **kwargs):
    #    pars = guess_from_peak(self, data, x, negative)
    #    return update_param_vals(pars, self.prefix, **kwargs)

#    __init__.__doc__ = COMMON_INIT_DOC
#    guess.__doc__ = COMMON_GUESS_DOC


class VinetModel(lmfit.Model):
    """
    lmfit Model class for Vinet fitting
    """

    def __init__(self, independent_vars=['v'], param_names=['v0', 'k0', 'k0p'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3
        :param param_names: define parameter names, v0, k0, k0p
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(VinetModel, self).__init__(vinet_p, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('k0', min=0.)
        self.set_param_hint('k0p', min=0.)


class KuncModel(lmfit.Model):
    """
    lmfit Model class for Kunc fitting
    """

    def __init__(self, independent_vars=['v'], param_names=['v0', 'k0', 'k0p'],
                 prefix='', missing=None, name=None, **kwargs):
        """
        :param independent_vars: define independent variables for lmfit
            unit-cell volume in A^3
        :param param_names: define parameter names, v0, k0, k0p
        :param prefix: see lmfit
        :param missing: see lmfit
        :param name: see lmfit
        :param kwargs: see lmfit,
            particularly useful to define order for Kunc function
        """
        kwargs.update({'prefix': prefix, 'missing': missing,
                       'independent_vars': independent_vars,
                       'param_names': param_names})
        super(KuncModel, self).__init__(kunc_p, **kwargs)
        self.set_param_hint('v0', min=0.)
        self.set_param_hint('k0', min=0.)
        self.set_param_hint('k0p', min=0.)
