"""
Todo's
"""
import matplotlib.pyplot as plt
from uncertainties import unumpy as unp
import uncertainties as uct
import numpy as np


def static_fit_result(fit_result, v_residual=None,
                      v_label='Unit-cell volume $(\mathrm{\AA}^3)$',
                      figsize=(5, 5), height_ratios=(3, 1), ms_data=8,
                      p_err=None, v_err=None, pdf_filen=None,
                      title='Fit result'):
    """
    plot static compressional curve fitting result

    :param fit_result: lmfit result object, see example jnb file for detail
    :param v_residual: manual input of fit residual
    :param v_label: label for volume axis
    :param figsize: figure size
    :param fitline: manual input of fit line, array of (v_fitline, p_fitline)
        with uncertainties
    :param height_ratios: height ratio between the main and residue plots
    :param ms_data: marker size for data points
    :param p_err: pressure error bar
    :param v_err: volume error bar
    :param pdf_filen: name of pdf output file
    :param title: title of the figure
    :return: None
    """
    # basic figure setup
    f, ax = plt.subplots(2, 1, sharex=True, figsize=figsize,
                         gridspec_kw={'height_ratios': height_ratios})
    for ax_i in ax:
        ax_i.tick_params(direction='in')
    # read data to plot
    v_data = fit_result.userkws['v']
    p_data = fit_result.data
    p_datafit = fit_result.best_fit
    v0 = uct.ufloat(fit_result.params['v0'].value,
                    fit_result.params['v0'].stderr)
    k0 = uct.ufloat(fit_result.params['k0'].value,
                    fit_result.params['k0'].stderr)
    k0p = uct.ufloat(fit_result.params['k0p'].value,
                     fit_result.params['k0p'].stderr)
    # I don't know why but fit_result.residual is not p_data - p_fit
    # p_residual = fit_result.residual
    # setup fitline and uncertainties
    v_fitline = np.linspace(v0.n, min(v_data), 1000)
    p_fit = fit_result.model.func(v_fitline, v0, k0, k0p)
    p_fitline = unp.nominal_values(p_fit)
    # I am unsure about why I have to divide by 3 in order to fit the
    # eval_uncertainty result from lmfit.  For now, I do this to be
    # consistent
    del_p = unp.std_devs(p_fit) / 3.
    # v_fitline = v_data
    # p_fitline = p_datafit
    # del_p = fit_result.eval_uncertainty()
    # ax[0].plot(p_data, v_data, 'ko', ms=ms_data, label='Data')
    ax[0].errorbar(p_data, v_data, xerr=p_err, yerr=v_err, fmt='ko',
                   ms=ms_data, mec='w', capsize=0, elinewidth=0.5,
                   label='Data')
    ax[0].plot(p_fitline, v_fitline, 'k-', label='Fit')
    ax[0].fill_betweenx(v_fitline, p_fitline - del_p,
                        p_fitline + del_p, color="k", alpha=0.2)
    if v_residual is None:
        # ax[1].plot(p_data, p_data - fit_result.best_fit, 'ko', ms=ms_data)
        ax[1].errorbar(p_data, p_data - p_datafit, yerr=p_err,
                       fmt='ko', ms=ms_data, mec='w',
                       capsize=0, elinewidth=0.5)
        ax[1].set_ylabel('$P_{obs} - P_{fit}$')
        ax[1].fill_between(p_fitline, -1. * del_p,
                           del_p, color="k", alpha=0.2)
    else:
        # ax[1].plot(p_data, v_data - v_fit, 'ko', ms=ms_data)
        ax[1].errorbar(p_data, v_residual, yerr=v_err,
                       fmt='ko', ms=ms_data, mec='w',
                       capsize=0, elinewidth=0.5)
        ax[1].set_ylabel('$V_{obs} - V_{fit}$')
    ax[0].legend()
    ax[1].axhline(0, c='k', ls='--')
    ax[1].set_xlabel('Pressure (GPa)')
    ax[0].set_ylabel(v_label)
    ax[0].set_title(title)
    plt.tight_layout()
    if pdf_filen is not None:
        f.savefig(pdf_filen)
