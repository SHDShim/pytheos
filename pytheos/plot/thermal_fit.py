"""
Todo's
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import copy
import numpy as np
import uncertainties as uct
from uncertainties import unumpy as unp
from ..etc import isuncertainties

c_map = cm.jet  # mpl.cm.inferno #mpl.cm.CMRmap_r # mpl.cm.brg #mpl.cm.gray_r #


def thermal_data(data, figsize=(12, 4), ms_data=50,
                 v_label='Unit-cell volume $(\mathrm{\AA}^3)$',
                 pdf_filen=None, title='P-V-T data'):
    """
    plot P-V-T data before fitting

    :param data: {'p': unumpy array, 'v': unumpy array, 'temp': unumpy array}
    :param eoscurves: {'v': unumpy array, '300': unumpy array
            at the temperature ....}
    :param v_label: label for volume axis
    :param figsize: figure size
    :param ms_data: marker size for data points
    :param pdf_filen: name of pdf output file
    :param title: title of the figure
    :return: None
    """
    # basic figure setup
    f, ax = plt.subplots(1, 2, figsize=figsize, sharex=True)

    # read data to plot
    if isuncertainties([data['p'], data['v'], data['temp']]):
        p = unp.nominal_values(data['p'])
        v = unp.nominal_values(data['v'])
        temp = unp.nominal_values(data['temp'])
        sp = unp.std_devs(data['p'])
        sv = unp.std_devs(data['v'])
        stemp = unp.std_devs(data['temp'])
        ax[0].errorbar(p, v, xerr=sp, yerr=sv, marker=' ',
                       c='k', ms=0, mew=0, linestyle='None',
                       capsize=0, lw=0.5, zorder=1)
        ax[1].errorbar(p, temp, xerr=sp, yerr=stemp, marker=' ',
                       c='k', ms=0, mew=0, linestyle='None',
                       capsize=0, lw=0.5, zorder=1)
    else:
        p = data['p']
        v = data['v']
        temp = data['temp']
    points = ax[0].scatter(p, v, marker='o', s=ms_data, c=temp,
                           cmap=c_map, vmin=300., vmax=temp.max(), zorder=2)
    points = ax[1].scatter(p, temp, marker='o', s=ms_data, c=temp,
                           cmap=c_map, vmin=300., vmax=temp.max(), zorder=2)

    ax[0].set_xlabel('Pressure (GPa)')
    ax[1].set_xlabel('Pressure (GPa)')
    ax[0].set_ylabel(v_label)
    ax[1].set_ylabel('Temperature (K)')
    f.suptitle(title)
    # the parameters are the specified position you set
    position = f.add_axes([0.92, 0.11, .01, 0.75])
    f.colorbar(points, orientation="vertical", cax=position)
    # position.text(150., 0.5, 'Temperature (K)', fontsize=10,
    # rotation=270, va='center')
    if pdf_filen is not None:
        f.savefig(pdf_filen)


c_map = cm.jet


def thermal_fit_result(fit_result, v_residual=None,
                       v_label='Unit-cell volume $(\mathrm{\AA}^3)$',
                       temp_fitline=np.asarray(
                           [300., 1000., 1500., 2000., 2500., 3000.]),
                       figsize=(5, 5), height_ratios=(3, 1), ms_data=50,
                       p_err=None, v_err=None, cbar_loc=(0.99, 0.1, .01, 0.82),
                       pdf_filen=None, title='Fit result'):
    """
    plot P-V-T EOS curve fitting result

    :param fit_result: lmfit result object, see example jnb file for detail
    :param v_label: label for volume axis
    :param temp_fitline: temperatures to calculate isothermal compression
        curves, default = [300., 1000., 1500., 2000., 2500., 3000.]
    :param figsize: figure size, default = (7,7)
    :param height_ratios: height ratio between the main and residue plots,
        default = (3,1)
    :param ms_data: marker size for data points
    :param p_err: pressure error bar
    :param v_err: volume error bar
    :param cbar_loc: location of color bar
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
    temp_data = fit_result.userkws['temp']
    p_data = fit_result.data
    p_datafit = fit_result.best_fit
    v0 = uct.ufloat(fit_result.params['st_v0'].value,
                    fit_result.params['st_v0'].stderr)
    sm = plt.cm.ScalarMappable(cmap=c_map,
                               norm=plt.Normalize(
                                   vmin=300., vmax=temp_data.max()))
    a = sm.to_rgba(temp_fitline)
    v_fitline = np.linspace(v0.n, min(v_data), 1000)
    fitmodel_copy = copy.deepcopy(fit_result)
    for a_i, temp_i in zip(a, temp_fitline):
        p_fitline = fitmodel_copy.eval(v=v_fitline,
                                       temp=np.ones_like(v_fitline) * temp_i)
        ax[0].plot(p_fitline, v_fitline, c=a_i)
    # error range here does not make a lot sense, so not supported
    # if (p_err is not None) and (v_err is not None):
    ax[0].errorbar(p_data, v_data, xerr=p_err, yerr=v_err, fmt=' ', c='k',
                   capsize=0, elinewidth=0.5, label='Data', zorder=0)
    points = ax[0].scatter(p_data, v_data, marker='o', s=ms_data, c=temp_data,
                           cmap=c_map, vmin=300., vmax=temp_data.max(),
                           zorder=1)
    if v_residual is None:
        ax[1].scatter(p_data, p_data - p_datafit, marker='o', s=ms_data,
                      c=temp_data, cmap=c_map, vmin=300.,
                      vmax=temp_data.max(), zorder=1)
        ax[1].errorbar(p_data, p_data - p_datafit, yerr=p_err, fmt=' ', c='k',
                       capsize=0, elinewidth=0.5, label='Data', zorder=0)
        ax[1].set_ylabel('$P_{obs} - P_{fit}$')
    else:
        ax[1].scatter(p_data, v_residual, marker='o', s=ms_data, c=temp_data,
                      cmap=c_map, vmin=300., vmax=temp_data.max(), zorder=1)
        ax[1].errorbar(p_data, v_residual, yerr=p_err, fmt=' ', c='k',
                       capsize=0, elinewidth=0.5, label='Data', zorder=0)
        ax[1].set_ylabel('$V_{obs} - V_{fit}$')
    # ax[0].legend()
    position = f.add_axes(cbar_loc)
    f.colorbar(points, orientation="vertical", cax=position,
               ticks=temp_fitline)
    ax[1].axhline(0, c='k', ls='--')
    ax[1].set_xlabel('Pressure (GPa)')
    ax[0].set_ylabel(v_label)
    ax[0].set_title(title)
    plt.tight_layout()
    if pdf_filen is not None:
        f.savefig(pdf_filen)
