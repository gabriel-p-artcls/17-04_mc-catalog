
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
from scipy import stats
from ra_dec_map import ra_dec_plots
from kde_map import kde_2d, kde_1d


def ccc(l1, l2):
    '''
    Concordance correlation coefficient.
    See: https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    '''
    return 2 * np.cov(l1, l2)[0, 1] / (np.var(l1) + np.var(l2) +
                                       (np.mean(l1) - np.mean(l2)) ** 2)


def as_vs_lit_plots(pl_params):
    '''
    Generate ASteCA vs literature values plots.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, xarr, xsigma, yarr, \
        ysigma, zarr, v_min_mp, v_max_mp, par_mean_std, gal_name = pl_params

    # Different limits for \delta plots.
    if i in [2, 5, 8, 11, 14]:
        ax = plt.subplot(gs[i], aspect='auto')
        # 0 line
        plt.axhline(y=par_mean_std[0], xmin=0, xmax=1, color='k', ls='--')
        # Shaded one sigma region.
        if par_mean_std[0] != par_mean_std[1]:
            plt.axhspan(par_mean_std[0] - par_mean_std[1],
                        par_mean_std[0] + par_mean_std[1], facecolor='grey',
                        alpha=0.5, zorder=1)
    else:
        ax = plt.subplot(gs[i], aspect='equal')
        # 1:1 line
        plt.plot([xmin, xmax], [ymin, ymax], 'k', ls='--')

    xy_font_s = 21
    plt.xticks(fontsize=xy_font_s - 6)
    plt.yticks(fontsize=xy_font_s - 6)
    cm = plt.cm.get_cmap('RdYlBu_r')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    if i in [1, 4, 7, 10, 13]:
        ax.set_yticklabels([])

    # Introduce random scatter.
    if i in [0, 1, 2, 9, 10, 11]:
        # 2% of [Fe/H], distance modulus axis ranges.
        ax_ext = (xmax - xmin) * 0.02
    else:
        # No scatter.
        ax_ext = 0.
    # Add random scatter.
    rs_x = xarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))
    rs_y = yarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))

    # Plot all clusters in dictionary.
    SC = plt.scatter(rs_x, rs_y, marker='o', c=zarr, s=75, lw=0.25, cmap=cm,
                     vmin=v_min_mp, vmax=v_max_mp, zorder=3)
    # Plot error bars.
    for j, xy in enumerate(zip(*[rs_x, rs_y])):
        # Only plot error bar if it has a value assigned in the literature.
        if ysigma:
            if ysigma[j] > 0. and xsigma[j] > 0.:
                plt.errorbar(xy[0], xy[1], xerr=xsigma[j], yerr=ysigma[j],
                             ls='none', color='k', elinewidth=0.5, zorder=1)
            elif xsigma[j] > 0. and ysigma[j] < 0.:
                plt.errorbar(xy[0], xy[1], xerr=xsigma[j],
                             ls='none', color='k', elinewidth=0.5, zorder=1)
            elif ysigma[j] > 0. and xsigma[j] <= 0.:
                plt.errorbar(xy[0], xy[1], yerr=ysigma[j], ls='none',
                             color='k', elinewidth=0.5, zorder=1)
    if i in [0, 1, 3, 4, 6, 7, 9, 10]:
        # Text box.
        ob = offsetbox.AnchoredText(gal_name, loc=4,
                                    prop=dict(size=xy_font_s - 4))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    if i in [2, 5, 8, 11]:
        # Text box.
        text1 = r'$\bar{{\Delta}}={:g}$'.format(round(par_mean_std[0], 2))
        text2 = r'$\pm{:g}$'.format(round(par_mean_std[1], 2))
        text = text1 + text2
        ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=xy_font_s - 4))
        ob.patch.set(alpha=0.5)
        ax.add_artist(ob)
        # Position colorbar.
        the_divider = make_axes_locatable(ax)
        color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
        # Colorbar.
        cbar = plt.colorbar(SC, cax=color_axis)
        cbar.set_label(z_lab, fontsize=xy_font_s, labelpad=7)


def make_as_vs_lit_plot(in_params):
    '''
    Prepare parameters and call function to generate ASteca vs literature
    SMC and LMC plots.
    '''

    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, rarr =\
        [in_params[_] for _ in
         ['zarr', 'zsigma', 'aarr', 'asigma', 'earr', 'esigma', 'darr',
          'dsigma', 'rarr']]

    z_all, age_all, ext_all, dm_all = [], [], [], []
    z_delta, age_delta, ext_delta, dm_delta = [], [], [], []
    z_delta_e, age_delta_e, ext_delta_e, dm_delta_e = [], [], [], []
    # SMC/LMC
    for k in [0, 1]:
        # \delta z as ASteCA - literature values.
        z_all += zarr[k][0]
        z_delta += list(np.array(zarr[k][0]) - np.array(zarr[k][1]))
        z_delta_e += list(np.sqrt(np.array(zsigma[k][0])**2 +
                          np.array(zsigma[k][1])**2))
        # \delta log(age) as ASteCA - literature values.
        age_all += aarr[k][0]
        age_delta += list(np.array(aarr[k][0]) - np.array(aarr[k][1]))
        age_delta_e += list(np.sqrt(np.array(asigma[k][0])**2 +
                            np.array(asigma[k][1])**2))
        # \delta E(B-V) as ASteCA - literature values.
        ext_all += earr[k][0]
        ext_delta += list(np.array(earr[k][0]) - np.array(earr[k][1]))
        ext_delta_e += list(np.sqrt(np.array(esigma[k][0])**2 +
                            np.array(esigma[k][1])**2))
        # \delta dm as ASteCA - literature values.
        dm_all += darr[k][0]
        dm_delta += list(np.array(darr[k][0]) - np.array(darr[k][1]))
        dm_delta_e += list(np.sqrt(np.array(dsigma[k][0])**2 +
                           np.array(dsigma[k][1])**2))

    # Replace huge error values for clusters with no extinctions in the
    # literature, with small values. Else they are plotted in the 1:1 plot.
    ext_delta_e = [min(99., abs(_)) for _ in ext_delta_e]

    # Gal Mean & StandDev
    par_mean_std = []
    for span in [z_delta, age_delta, ext_delta, dm_delta]:
        # Filter out -9999999999.9 values added in get_params to replace
        # missing values in .ods file.
        span_filter = []
        for _ in span:
            if abs(_) < 30000:
                span_filter.append(_)
        if span_filter:
            p_mean, p_stdev = np.mean(span_filter), np.std(span_filter)
            par_mean_std.append([p_mean, p_stdev])
            # print p_mean, p_stdev
        else:
            par_mean_std.append([0., 0.])

    # Generate ASteca vs literature plots.
    fig = plt.figure(figsize=(21, 31.25))  # create the top-level container
    # gs = gridspec.GridSpec(2, 4, width_ratios=[1, 0.35, 1, 0.35])
    gs = gridspec.GridSpec(5, 3)

    z_min, z_max = -1.3, 0.
    a_min, a_max = 6.6, 9.8
    ext_min, ext_max = 0., 0.25
    # dm_min, dm_max = 18.62, 19.21
    # dm_min, dm_max = 18.21, 18.79
    dm_min, dm_max = 18.21, 19.19
    dm_span = (dm_max - dm_min) / 2.

    # For old runs where the dist mod range was large.
    # dm_min, dm_max = 17.8, 20.2

    as_lit_pl_lst = [
        # Metallicity LMC/SMC
        [gs, 0, -2.4, 0.45, -2.4, 0.45, '$[Fe/H]_{\mathtt{ASteCA}}$',
         '$[Fe/H]_{lit}$', '', zarr[1][0], zsigma[1][0], zarr[1][1],
            zsigma[1][1], darr[1][0], dm_min, dm_max, [], 'LMC'],
        [gs, 1, -2.4, 0.45, -2.4, 0.45, '$[Fe/H]_{\mathtt{ASteCA}}$', '',
            '', zarr[0][0], zsigma[0][0], zarr[0][1],
            zsigma[0][1], darr[0][0], dm_min, dm_max, [], 'SMC'],
        # Asteca z vs \delta z with lit values.
        [gs, 2, -2.4, 0.45, -1.83, 1.43, '$[Fe/H]_{\mathtt{ASteCA}}$',
            '$\Delta [Fe/H]$', '$\mu_{0; \mathtt{ASteCA}}$', z_all,
            [0.]*len(z_all), z_delta, z_delta_e, dm_all, dm_min, dm_max,
            par_mean_std[0], ''],

        # Age LMC/SMC
        [gs, 3, 5.8, 10.6, 5.8, 10.6, '$\log(age/yr)_{\mathtt{ASteCA}}$',
            '$\log(age/yr)_{lit}$', '', aarr[1][0],
            asigma[1][0], aarr[1][1], asigma[1][1], zarr[1][0], z_min,
            z_max, [], 'LMC'],
        [gs, 4, 5.8, 10.6, 5.8, 10.6, '$\log(age/yr)_{\mathtt{ASteCA}}$',
            '', '', aarr[0][0], asigma[0][0], aarr[0][1], asigma[0][1],
            zarr[0][0], z_min, z_max, [], 'SMC'],
        # Asteca log(age) vs \delta log(age) with lit values.
        [gs, 5, 5.8, 10.6, -2.4, 2.4, '$\log(age/yr)_{\mathtt{ASteCA}}$',
            '$\Delta \log(age/yr)$', '$[Fe/H]_{\mathtt{ASteCA}}$', age_all,
            [0.]*len(age_all), age_delta, age_delta_e, z_all, z_min,
            z_max, par_mean_std[1], ''],

        # Ext LMC/SMC
        [gs, 6, -0.04, 0.29, -0.04, 0.29, '$E(B-V)_{\mathtt{ASteCA}}$',
            '$E(B-V)_{lit}$', '', earr[1][0], esigma[1][0], earr[1][1],
            esigma[1][1], aarr[1][0], a_min, a_max, [], 'LMC'],
        [gs, 7, -0.04, 0.29, -0.04, 0.29, '$E(B-V)_{\mathtt{ASteCA}}$',
            '', '', earr[0][0], esigma[0][0], earr[0][1], esigma[0][1],
            aarr[0][0], a_min, a_max, [], 'SMC'],
        # Asteca E(B-V) vs \delta E(B-V) with lit values.
        [gs, 8, -0.04, 0.29, -0.21, 0.21, '$E(B-V)_{\mathtt{ASteCA}}$',
            '$\Delta E(B-V)$', '$\log(age/yr)_{\mathtt{ASteCA}}$', ext_all,
            [0.]*len(ext_all), ext_delta, ext_delta_e, age_all, a_min, a_max,
            par_mean_std[2], ''],

        # Dits mod LMC/SMC
        [gs, 9, dm_min, dm_max, dm_min, dm_max, '$\mu_{0;\,\mathtt{ASteCA}}$',
            '$\mu_{0;\,lit}$', '', darr[1][0], dsigma[1][0], darr[1][1],
            dsigma[1][1], earr[1][0], ext_min, ext_max, [], 'LMC'],
        [gs, 10, dm_min, dm_max, dm_min, dm_max, '$\mu_{0;\,\mathtt{ASteCA}}$',
            '', '', darr[0][0], dsigma[0][0], darr[0][1], dsigma[0][1],
            earr[0][0], ext_min, ext_max, [], 'SMC'],
        # Asteca dist_mod vs \delta dist_mod with lit values.
        [gs, 11, dm_min, dm_max, -1. * dm_span, dm_span,
            '$\mu_{0;\,\mathtt{ASteCA}}$', '$\Delta \mu_{0}$',
            '$E(B-V)_{\mathtt{ASteCA}}$', dm_all, [0.]*len(dm_all),
            dm_delta, dm_delta_e, ext_all, ext_min, ext_max, par_mean_std[3],
            ''],
    ]

    for pl_params in as_lit_pl_lst:
        as_vs_lit_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_vs_lit_S-LMC.png', dpi=300, bbox_inches='tight')


def as_vs_lit_mass_plots(pl_params):
    '''
    Generate ASteCA vs literature mass values plot.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, xarr, xsigma, yarr, \
        ysigma, carr, v_min_mp, v_max_mp, par_mean_std, gal_name = pl_params

    xy_font_s = 21
    ax = plt.subplot(gs[i], aspect='auto')
    # Different line plots.
    if i == 0:
        ax.set_xticklabels([])
        # 1:1 line
        plt.plot([xmin, xmax], [ymin, ymax], 'k', ls='--')
        # Text box.
        ob = offsetbox.AnchoredText(gal_name, loc=1,
                                    prop=dict(size=xy_font_s - 4))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    else:
        # ax = plt.subplot(gs[i], aspect='auto')
        # 0 line
        plt.axhline(y=par_mean_std[0], xmin=0, xmax=1, color='k', ls='--')
        # Shaded one sigma region.
        if par_mean_std[0] != par_mean_std[1]:
            plt.axhspan(par_mean_std[0] - par_mean_std[1],
                        par_mean_std[0] + par_mean_std[1], facecolor='grey',
                        alpha=0.5, zorder=1)

    plt.xticks(fontsize=xy_font_s - 6)
    plt.yticks(fontsize=xy_font_s - 6)
    cm = plt.cm.get_cmap('RdYlBu_r')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()

    # Introduce random scatter. 1% of axis ranges.
    ax_ext = (xmax - xmin) * 0.01
    # Add random scatter.
    rs_x = xarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))
    rs_y = yarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))

    # Plot all clusters in dictionary.
    SC = plt.scatter(rs_x, rs_y, marker='o', c=carr, s=110, lw=0.25,
                     cmap=cm, vmin=v_min_mp, vmax=v_max_mp, zorder=3)
    # Plot error bars.
    for j, xy in enumerate(zip(*[rs_x, rs_y])):
        if i == 0:
            plt.errorbar(xy[0], xy[1], xerr=xsigma[j], yerr=ysigma[j],
                         ls='none', color='k', elinewidth=0.5, zorder=1)
        else:
            plt.errorbar(xy[0], xy[1], yerr=ysigma[j], ls='none',
                         color='k', elinewidth=0.5, zorder=1)
    if i != 0:
        # Text box.
        pres = [2, 2] if x_lab != '$M_{\mathtt{ASteCA}}\,[M_{\odot}]$' else\
            [0, 0]
        text1 = r'$\bar{{\Delta}}={:g}$'.format(round(par_mean_std[0],
                                                pres[0]))
        text2 = r'$\pm{:g}$'.format(round(par_mean_std[1], pres[1]))
        text = text1 + text2
        ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=xy_font_s - 4))
        ob.patch.set(alpha=0.5)
        ax.add_artist(ob)
    # Position colorbar.
    # the_divider = make_axes_locatable(ax)
    # color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # # Colorbar.
    # cbar = plt.colorbar(SC, cax=color_axis)
    # cbar.set_label(z_lab, fontsize=xy_font_s, labelpad=7)
    # cbar.set_ticks([0., 0.2, 0.4, 0.6, 0.8, 1., 1.2])
    if i == 0:
        # Position colorbar.
        axColor = plt.axes([0.25, 0.83, 0.1, 0.005])
        cbar = plt.colorbar(SC, cax=axColor, orientation="horizontal")
        cbar.set_label(z_lab, fontsize=xy_font_s - 2, labelpad=-45)
        cbar.set_ticks([0.6, 0.8, 1., 1.2])
        cbar.ax.tick_params(labelsize=xy_font_s - 10)


def make_as_vs_lit_mass_plot(in_params):
    '''
    Prepare parameters and call function to generate ASteca vs literature
    SMC and LMC plots.
    '''

    aarr, marr, msigma, zarr, int_colors, cont_ind, phot_disp = [
        in_params[_] for _ in [
            'aarr', 'marr', 'msigma', 'zarr', 'int_colors', 'cont_ind',
            'phot_disp']]

    # SMC/LMC
    age_all, ci_all, ma_all, ma_delta, ma_delta_err = [], [], [], [], []
    for k in [0, 1]:
        for i, a in enumerate(aarr[k][0]):
            # Filter out -9999999999.9 values added in get_params to replace
            # missing values in .ods file.
            if abs(msigma[k][1][i]) < 30000.:
                age_all.append(a)
                ci_all.append(cont_ind[k][i])
                # \delta mass as ASteCA - literature values.
                ma_all.append(marr[k][0][i])
                ma_delta.append(marr[k][0][i] - marr[k][1][i])
                ma_delta_err.append(np.sqrt(msigma[k][0][i]**2 +
                                            msigma[k][1][i]**2))

    # Mean & StandDev
    par_mean_std = [np.mean(ma_delta), np.std(ma_delta)]

    # Replace huge mass error values for clusters with no masses in the
    # literature, with small values. Else they are plotted a the bottom
    # of the 1:1 plot.
    msigma_lit = [min(999, abs(_)) for _ in msigma[0][1]]

    # Generate ASteca vs literature plots.
    fig = plt.figure(figsize=(18.5, 31.25))  # create the top-level container
    gs = gridspec.GridSpec(5, 3)

    cbar_min, cbar_max = 0.6, 1.3

    as_lit_pl_lst = [
        # ASteCA vs literature masses.
        [gs, 0, 10., 3998., 10., 3998., '', '$M_{lit}\,[M_{\odot}]$',
         '$CI$', marr[0][0], msigma[0][0], marr[0][1], msigma_lit,
         cont_ind[0], cbar_min, cbar_max, [], 'SMC (M13)'],
        [gs, 3, 10., 3998., -3500., 1900.,
         '$M_{\mathtt{ASteCA}}\,[M_{\odot}]$',
         '$\Delta M\,[M_{\odot}]$', '$CI$', ma_all, [0.]*len(ma_all),
         ma_delta, ma_delta_err, ci_all, cbar_min, cbar_max, par_mean_std, '']
    ]

    for pl_params in as_lit_pl_lst:
        as_vs_lit_mass_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_vs_lit_mass.png', dpi=300, bbox_inches='tight')


# import matplotlib.colors as mcolors
# def make_colormap(seq):
#     """Return a LinearSegmentedColormap
#     seq: a sequence of floats and RGB-tuples. The floats should be increasing
#     and in the interval (0,1).

#     Source: http://stackoverflow.com/a/16836182/1391441
#     Color names: http://stackoverflow.com/a/29676907/1391441
#     """
#     seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
#     cdict = {'red': [], 'green': [], 'blue': []}
#     for i, item in enumerate(seq):
#         if isinstance(item, float):
#             r1, g1, b1 = seq[i - 1]
#             r2, g2, b2 = seq[i + 1]
#             cdict['red'].append([item, r1, r2])
#             cdict['green'].append([item, g1, g2])
#             cdict['blue'].append([item, b1, b2])
#     return mcolors.LinearSegmentedColormap('CustomMap', cdict)


def kde_plots(pl_params):
    '''
    Generate KDE plots.
    '''
    gs, i, gs_pos, x_lab, y_lab, xarr, xsigma, yarr, ysigma, x_rang, y_rang,\
        size = pl_params

    # Make plot.
    a, b, c, d = gs_pos[i]
    ax = plt.subplot(gs[a:b, c:d])
    xy_font_s = 18
    plt.xlabel(x_lab, fontsize=xy_font_s)
    if i == 6:
        # ax.xaxis.set_label_position('top')
        ax.xaxis.set_label_coords(0.5, 1.06)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.minorticks_on()

    if i not in [4, 5]:
        ax.set_xticklabels([])
    if i not in [2, 4]:
        ax.set_yticklabels([])

    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=-1)
    # Grid density for the KDE evaluation.
    grid_dens = 100

    if i in [0, 1]:
        # Generate map.
        col, lab = ['r', 'b'], ['SMC', 'LMC']
        x, y = kde_1d(np.array(xarr), np.array(xsigma), x_rang, grid_dens)
        # print 'int age:', integrate.simps(y, x)
        # # Confidence intervals.
        # lcb, ucf = lf_cb.monte_carlo_conf_bands(xarr, xsigma, x_rang,
        #                                         grid_dens)
        # plt.fill_between(x, lcb, ucf, alpha=0.3, facecolor=col[i])
        ax.plot(x, y, color=col[i], label=lab[i])

        # Obtain Gaussian fit for 1D KDE.
        # from scipy.optimize import curve_fit
        # def gaus(x, x0, sigma):
        #     return (1/(sigma*np.sqrt(2*np.pi))) *\
        #         np.exp(-0.5*((x - x0)/sigma)**2)

        # mean, sigma = np.mean(x), np.std(x)
        # popt, pcov = curve_fit(gaus, x, y, p0=[mean, sigma])
        # print 'x +- std:', popt[0], popt[1]

    elif i in [6, 7]:
        # Generate map.
        x, y = kde_1d(np.array(xarr), np.array(xsigma), y_rang, grid_dens)
        ax.plot(y, x, color='r', label='SMC')
        # print 'int met:', integrate.simps(y, x)
        # # Confidence intervals.
        # lcb, ucf = lf_cb.monte_carlo_conf_bands(xarr, xsigma, y_rang,
        #                                         grid_dens)
        # plt.fill_betweenx(x, lcb, ucf, alpha=0.3, facecolor='r')
        max_y = max(y)
        #
        x, y = kde_1d(np.array(yarr), np.array(ysigma), y_rang, grid_dens)
        ax.plot(y, x, color='b', label='LMC')
        # print 'int: mass', integrate.simps(y, x)
        # # Confidence intervals.
        # lcb, ucf = lf_cb.monte_carlo_conf_bands(yarr, ysigma, y_rang,
        #                                         grid_dens)
        # plt.fill_betweenx(x, lcb, ucf, alpha=0.3, facecolor='b')
        max_y = max(max(y), max_y) + 0.1*max(max(y), max_y)
    else:
        col = 'r' if i % 2 == 0 else 'b'
        # Generate map.
        ext = [x_rang[0], x_rang[1], y_rang[0], y_rang[1]]
        z = kde_2d(np.array(xarr), np.array(xsigma), np.array(yarr),
                   np.array(ysigma), ext, grid_dens)
        cm = plt.cm.gist_earth_r
        # c = mcolors.ColorConverter().to_rgb
        # cm = make_colormap(
        #     [c('white'), 0.01, c('mintcream'), c('palegoldenrod'), 0.3,
        #      c('wheat'), 0.4, c('khaki'), 0.5, c('darkkhaki'), 0.6,
        #      c('yellowgreen'), 0.7, c('mediumseagreen'), 0.8, c('seagreen'),
        #      0.9, c('green'), 0.95, c('darkgreen'), 0.97, c('black')])
        ax.imshow(z, cmap=cm, extent=ext)
        ax.set_aspect('auto')
        # Error bars.
        # plt.errorbar(xarr, yarr, xerr=xsigma, yerr=ysigma, fmt='none',
        #              elinewidth=0.4, color='k')
        # Define x% of axis ranges.
        xax_ext = (ext[1] - ext[0]) * 0.02
        yax_ext = (ext[3] - ext[2]) * 0.02
        # Random scatter.
        rs_x = np.random.uniform(0., xax_ext, len(xarr))
        rs_y = np.random.uniform(0., yax_ext, len(xarr))
        # Clusters.
        # color='#6b6868'
        plt.scatter(xarr + rs_x, yarr + rs_y, marker='*', color=col,
                    s=11.*np.array(size), lw=0.45, facecolors='none')
    if i in [0, 1]:
        leg = plt.legend(loc='upper left', markerscale=1.,
                         scatterpoints=1, fontsize=xy_font_s - 4)
        leg.get_frame().set_alpha(0.5)
    # Text box
    if i in [2, 3, 4, 5]:
        text = 'SMC' if i % 2 == 0 else 'LMC'
        l = 3 if (i-4) < 0 else 2
        ob = offsetbox.AnchoredText(text, loc=l, prop=dict(size=xy_font_s - 4))
        ob.patch.set(alpha=0.5)
        ax.add_artist(ob)
    # Limits.
    if i in [6, 7]:
        ax.set_xlim(x_rang[0], max_y)
    else:
        ax.set_xlim(x_rang[0], x_rang[1])
    ax.set_ylim(y_rang[0], y_rang[1])


def make_kde_plots(in_params):
    '''
    Prepare parameters and call function to generate SMC and LMC KDE plots.
    '''
    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, marr, msigma,\
        rad_pc = \
        [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
                                'esigma', 'darr', 'dsigma', 'marr', 'msigma',
                                'rad_pc']]

    fig = plt.figure(figsize=(17, 17))  # create the top-level container
    gs = gridspec.GridSpec(6, 6)  # create a GridSpec object

    # Logarithmic mass
    log_mass_smc, e_log_mass_smc = np.log10(marr[0][0]),\
        np.array(msigma[0][0])/np.array(marr[0][0])
    log_mass_lmc, e_log_mass_lmc = np.log10(marr[1][0]),\
        np.array(msigma[1][0])/np.array(marr[1][0])

    # Define extension for each parameter range.
    age_rang, fe_h_rang, log_mass_rang, dist_rang, ext_rang = [6.51, 10.1],\
        [-2.4, 0.25], [1.2, 4.9], [[18.81, 19.14], [18.36, 18.67]],\
        [-0.015, 0.31]
    age_kde_rang, feh_kde_rang, log_m_kde_rang, dist_kde_rang, ext_kde_rang =\
        [0., 1.27], [0., 2.], [0., 2.], [0., 5.1], [0., 2.]

    gs_pos = [[1, 2, 0, 2], [1, 2, 2, 4], [2, 4, 0, 2], [2, 4, 2, 4],
              [4, 6, 0, 2], [4, 6, 2, 4], [2, 4, 4, 5], [4, 6, 4, 5]]

    kde_pl_lst = [
        # SMC
        [gs, 0, gs_pos, '', r'$KDE_{\,\log(age/yr)}$', aarr[0][0],
         asigma[0][0], [], [], age_rang, age_kde_rang, []],
        # LMC
        [gs, 1, gs_pos, '', '', aarr[1][0], asigma[1][0], [], [], age_rang,
         age_kde_rang, []],
        #
        [gs, 2, gs_pos, '', '$[Fe/H]_{\mathtt{ASteCA}}$', aarr[0][0],
         asigma[0][0], zarr[0][0], zsigma[0][0], age_rang, fe_h_rang,
         rad_pc[0]],
        [gs, 3, gs_pos, '', '', aarr[1][0], asigma[1][0],
         zarr[1][0], zsigma[1][0], age_rang, fe_h_rang, rad_pc[1]],
        #
        [gs, 4, gs_pos, '$\log(age/yr)_{\mathtt{ASteCA}}$',
         '$\log(M/M_{\odot})_{\mathtt{ASteCA}}$', aarr[0][0], asigma[0][0],
         log_mass_smc, e_log_mass_smc, age_rang, log_mass_rang, rad_pc[0]],
        [gs, 5, gs_pos, '$\log(age/yr)_{\mathtt{ASteCA}}$', '', aarr[1][0],
         asigma[1][0], log_mass_lmc, e_log_mass_lmc, age_rang, log_mass_rang,
         rad_pc[1]],
        #
        # S/LMC [Fe/H] & log(Mass) KDEs
        [gs, 6, gs_pos, r'$KDE_{\,[Fe/H]}$', '', zarr[0][0], zsigma[0][0],
         zarr[1][0], zsigma[1][0], feh_kde_rang, fe_h_rang, []],
        [gs, 7, gs_pos, r'$KDE_{\,\log(M/M_{\odot})}$', '', log_mass_smc,
         e_log_mass_smc, log_mass_lmc, e_log_mass_lmc, log_m_kde_rang,
         log_mass_rang, []]
    ]
    #
    for pl_params in kde_pl_lst:
        kde_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_kde_maps0.png', dpi=300, bbox_inches='tight')
    fig.clf()

    kde_pl_lst = [
        # SMC
        [gs, 0, gs_pos, '', r'$KDE_{\,\mu_0}$', darr[0][0],
         dsigma[0][0], [], [], dist_rang[0], dist_kde_rang, []],
        # LMC
        [gs, 1, gs_pos, '', '', darr[1][0], dsigma[1][0], [], [], dist_rang[1],
         dist_kde_rang, []],
        #
        [gs, 4, gs_pos, '$\mu_{0;\mathtt{ASteCA}}$',
         '$E(B-V)_{\mathtt{ASteCA}}$', darr[0][0],
         dsigma[0][0], earr[0][0], esigma[0][0], dist_rang[0], ext_rang,
         rad_pc[0]],
        [gs, 5, gs_pos, '$\mu_{0;\mathtt{ASteCA}}$', '', darr[1][0],
         dsigma[1][0], earr[1][0], esigma[1][0], dist_rang[1], ext_rang,
         rad_pc[1]],
        #
        # S/LMC E(B-V) KDEs
        [gs, 7, gs_pos, r'$KDE_{\,E(B-V)}$', '', earr[0][0], esigma[0][0],
         earr[1][0], esigma[1][0], ext_kde_rang, ext_rang, []]
    ]
    #
    for pl_params in kde_pl_lst:
        kde_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_kde_maps1.png', dpi=300, bbox_inches='tight')


def make_ra_dec_plots(in_params, bica_coords):
    '''
    Prepare parameters and call function to generate RA vs DEC positional
    plots for the SMC and LMC.
    '''

    ra, dec, zarr, aarr, earr, darr, marr, rad_pc = [
        in_params[_] for _ in ['ra', 'dec', 'zarr', 'aarr', 'earr', 'darr',
                               'marr', 'rad_pc']]

    # Check plot.
    # cm = plt.cm.get_cmap('RdYlBu_r')
    # bb_ra, bb_dec = zip(*bica_coords)
    # plt.scatter(-1. * np.asarray(bb_ra), bb_dec, marker='.', s=10)
    # plt.scatter(-1. * np.asarray(ra[0]), dec[0], c=darr[0][0], cmap=cm,
    #             marker='o', s=50, lw=0.1, vmin=18.82, vmax=19.08)
    # plt.show()

    # Put both SMC and LMC clusters into a single list.
    ra = ra[0] + ra[1]
    dec = dec[0] + dec[1]
    zarr = zarr[0][0] + zarr[1][0]
    aarr = aarr[0][0] + aarr[1][0]
    earr = earr[0][0] + earr[1][0]
    darr = darr[0][0] + darr[1][0]
    marr = marr[0][0] + marr[1][0]
    rad_pc = rad_pc[0] + rad_pc[1]

    # Sort according to radius value so that larger clusters will be plotted
    # first.
    rad_pc, ra, dec, zarr, aarr, earr, darr, marr = \
        map(list, zip(*sorted(zip(rad_pc, ra, dec, zarr, aarr, earr, darr,
                                  marr), reverse=True)))

    # Bica coords.
    bb_ra, bb_dec = zip(*bica_coords)

    fig = plt.figure(figsize=(15, 20))
    fig.clf()

    ra_dec_pl_lst = [
        [fig, 321, ra, dec, bb_ra, bb_dec, zarr, -2.1, 0., rad_pc, '$[Fe/H]$'],
        [fig, 322, ra, dec, bb_ra, bb_dec, aarr, 6.6, 10., rad_pc,
         '$log(age/yr)$'],
        [fig, 323, ra, dec, bb_ra, bb_dec, earr, 0., 0.3, rad_pc,
         '$E_{(B-V)}$'],
        [fig, 324, ra, dec, bb_ra, bb_dec, darr, 18.4, 18.6, rad_pc,
         '$(m-M)_{\circ}$'],
        [fig, 325, ra, dec, bb_ra, bb_dec, darr, 18.82, 19.08, rad_pc,
         '$(m-M)_{\circ}$'],
        [fig, 326, ra, dec, bb_ra, bb_dec, marr, 100, 30000, rad_pc,
         '$M\,[M_{\odot}]$']
        # [fig, 326, ra, dec, bb_ra, bb_dec, rad_pc, rad_pc,
        # '$r_{clust}\,[pc]$']
    ]

    for pl_params in ra_dec_pl_lst:
        ra_dec_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_RA_DEC.png', dpi=300, bbox_inches='tight')


def lit_ext_plots(pl_params):
    '''
    Generate ASteCA vs literature values plots.
    '''
    gs, i, xmin, xmax, x_lab, y_lab, z_lab, xarr, xsigma, yarr, \
        ysigma, zarr, v_min_mp, v_max_mp, gal_name = pl_params

    xy_font_s = 18
    cm = plt.cm.get_cmap('RdYlBu_r')

    # Different limits for \delta plots.
    ax = plt.subplot(gs[i], aspect='equal')
    # 1:1 line
    plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')

    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()

    # # Introduce random scatter.
    # 5% of axis ranges.
    ax_ext = (xmax - xmin) * 0.05
    # # Add randoms scatter.
    rs_x = xarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))
    rs_y = yarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))

    # Plot all clusters in dictionary.
    SC = plt.scatter(rs_x, rs_y, marker='o', c=zarr, s=70, lw=0.25, cmap=cm,
                     vmin=v_min_mp, vmax=v_max_mp, zorder=3)
    # Plot error bars.
    for j, xy in enumerate(zip(*[rs_x, rs_y])):
        # Only plot error bar if it has a value assigned in the literature.
        if ysigma:
            if ysigma[j] > 0. and xsigma[j] > 0.:
                plt.errorbar(xy[0], xy[1], xerr=xsigma[j], yerr=ysigma[j],
                             ls='none', color='k', elinewidth=0.5, zorder=1)
            elif xsigma[j] > 0. and ysigma[j] < 0.:
                plt.errorbar(xy[0], xy[1], xerr=xsigma[j],
                             ls='none', color='k', elinewidth=0.5, zorder=1)
            elif ysigma[j] > 0. and xsigma[j] < 0.:
                plt.errorbar(xy[0], xy[1], yerr=ysigma[j], ls='none',
                             color='k', elinewidth=0.5, zorder=1)
    if gal_name != '':
        # Text box.
        ob = offsetbox.AnchoredText(gal_name, loc=4, prop=dict(size=xy_font_s))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    if i in [1, 3]:
        # Position colorbar.
        the_divider = make_axes_locatable(ax)
        color_axis = the_divider.append_axes("right", size="5%", pad=0.1)
        # Colorbar.
        cbar = plt.colorbar(SC, cax=color_axis)
        zpad = 10 if z_lab == '$E_{(B-V)}$' else 5
        cbar.set_label(z_lab, fontsize=xy_font_s, labelpad=zpad)


def make_lit_ext_plot(in_params):
    '''
    ASteca vs MCEV vs SandF extinction plot done.
    '''

    aarr, earr, esigma, ext_sf, ext_mcev = \
        [in_params[_] for _ in ['aarr', 'earr', 'esigma', 'ext_sf',
                                'ext_mcev']]

    # Order lists to put max distance values on top.
    # SMC
    ord_mcev_dist_smc, ord_earr_smc_ast, ord_esig_smc_ast, ord_mcev_smc, \
        ord_e_mcev_smc =\
        map(list, zip(*sorted(zip(ext_mcev[0][3], earr[0][0], esigma[0][0],
            ext_mcev[0][0], ext_mcev[0][2]), reverse=False)))
    # LMC
    ord_mcev_dist_lmc, ord_earr_lmc_ast, ord_esig_lmc_ast, ord_mcev_lmc, \
        ord_e_mcev_lmc =\
        map(list, zip(*sorted(zip(ext_mcev[1][3], earr[1][0], esigma[1][0],
            ext_mcev[1][0], ext_mcev[1][2]), reverse=False)))

    # Define values to pass.
    xmin, xmax = -0.014, [0.15, 0.4]
    vmin_sf_SMC, vmax_sf_SMC = min(ext_sf[0][0]), max(ext_sf[0][0])
    vmin_sf_LMC, vmax_sf_LMC = min(ext_sf[1][0]), max(ext_sf[1][0])
    x_lab = '$E(B-V)_{\mathtt{ASteCA}}$'
    y_lab = ['$E(B-V)_{MCEV,\,closer}$', '$E(B-V)_{MCEV,\,max}$']
    z_lab = ['$log(age/yr)_{\mathtt{ASteCA}}$', '$E(B-V)_{SF}$',
             '$dist\,(deg)$']

    fig = plt.figure(figsize=(16, 25))
    gs = gridspec.GridSpec(4, 2)

    ext_pl_lst = [
        # SMC
        [gs, 0, xmin, xmax[0], x_lab, y_lab[0], z_lab[2],
            ord_earr_smc_ast, ord_esig_smc_ast, ord_mcev_smc, ord_e_mcev_smc,
            ord_mcev_dist_smc, vmin_sf_SMC, vmax_sf_SMC, 'SMC'],
        [gs, 1, xmin, xmax[0], x_lab, y_lab[1], z_lab[1],
            earr[0][0], esigma[0][0], ext_mcev[0][1], ext_mcev[0][2],
            ext_sf[0][0], vmin_sf_SMC, vmax_sf_SMC, ''],
        # LMC
        [gs, 2, xmin, xmax[1], x_lab, y_lab[0], z_lab[2],
            ord_earr_lmc_ast, ord_esig_lmc_ast, ord_mcev_lmc, ord_e_mcev_lmc,
            ord_mcev_dist_lmc, vmin_sf_LMC, vmax_sf_LMC, 'LMC'],
        [gs, 3, xmin, xmax[1], x_lab, y_lab[1], z_lab[1],
            earr[1][0], esigma[1][0], ext_mcev[1][1], ext_mcev[1][2],
            ext_sf[1][0], vmin_sf_LMC, vmax_sf_LMC, '']
    ]

    for pl_params in ext_pl_lst:
        lit_ext_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_vs_lit_extin.png', dpi=300, bbox_inches='tight')


def wide_plots(pl_params):
    '''
    Generate plots for integrated colors, concentration parameter, and radius
    (in parsec) vs several parameters.
    '''
    gs, i, xlim, ylim, x_lab, y_lab, z_lab, xarr, xsigma, yarr, ysigma, zarr,\
        rad, gal_name = pl_params
    siz = np.asarray(rad) * 5

    xy_font_s = 16
    cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i])
    # ax.set_aspect('auto')
    xmin, xmax = xlim
    plt.xlim(xmin, xmax)
    if ylim:
        ymin, ymax = ylim
        plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    # Plot all clusters in dictionary.
    SC = plt.scatter(xarr, yarr, marker='o', c=zarr, s=siz, lw=0.25, cmap=cm,
                     zorder=3)
    # Plot x error bar.
    plt.errorbar(xarr, yarr, xerr=xsigma, ls='none', color='k',
                 elinewidth=0.4, zorder=1)
    # Plot y error bar if it is passed.
    if ysigma:
        plt.errorbar(xarr, yarr, yerr=ysigma, ls='none', color='k',
                     elinewidth=0.4, zorder=1)
    # Text box.
    ob = offsetbox.AnchoredText(gal_name, loc=2, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    zpad = 10 if z_lab == '$E_{(B-V)}$' else 5
    cbar.set_label(z_lab, fontsize=xy_font_s, labelpad=zpad)


def make_int_cols_plot(in_params):
    '''
    Prepare parameters and call function to generate integrated color vs Age
    (colored by mass) plots for the SMC and LMC.
    '''

    aarr, asigma, marr, int_colors, rad_pc = [
        in_params[_] for _ in ['aarr', 'asigma', 'marr', 'int_colors',
                               'rad_pc']]

    # Define values to pass.
    xmin, xmax = 6.5, 9.95
    x_lab, y_lab, z_lab = '$log(age/yr)_{\mathtt{ASteCA}}$', \
        '$(C-T_{1})_{0;\,\mathtt{ASteCA}}$', '$M\,[M_{\odot}]$'

    fig = plt.figure(figsize=(16, 25))  # create the top-level container
    gs = gridspec.GridSpec(4, 1)       # create a GridSpec object

    ext_pl_lst = [
        # SMC
        [gs, 0, [xmin, xmax], [], x_lab, y_lab, z_lab, aarr[0][0],
            asigma[0][0], int_colors[0], [], marr[0][0], rad_pc[0], 'SMC'],
        # LMC
        [gs, 1, [xmin, xmax], [], x_lab, y_lab, z_lab, aarr[1][0],
            asigma[1][0], int_colors[1], [], marr[1][0], rad_pc[1], 'LMC']
    ]

    for pl_params in ext_pl_lst:
        wide_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_integ_colors.png', dpi=300, bbox_inches='tight')


def make_concent_plot(in_params):
    '''
    Generate ASteCA concentration parameter (cp) plots, where:
    cp = n_memb / (r_pc **2)
    '''

    zarr, zsigma, aarr, asigma, marr, rad_pc, n_memb, rad_pc = [
        in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'marr',
                               'rad_pc', 'n_memb', 'rad_pc']]

    # Calculate the 'concentration parameter' as the approximate number of
    # (structural) members divided by the area of the cluster in parsecs.
    conc_p = [[], []]
    for j in [0, 1]:
        conc_p[j] = np.asarray(n_memb[j]) / \
            (np.pi * np.asarray(rad_pc[j]) ** 2)

    # Define values to pass.
    xmin, xmax = [6.5, -2.3], [10.4, 0.2]
    x_lab, y_lab, z_lab = ['$log(age/yr)_{\mathtt{ASteCA}}$',
                           '$[Fe/H]_{\mathtt{ASteCA}}$'], \
        '$Concentration\,(N_{memb}/pc^{2})$', '$M\,[M_{\odot}]$'

    fig = plt.figure(figsize=(16, 25))  # create the top-level container
    gs = gridspec.GridSpec(4, 1)       # create a GridSpec object

    conc_pl_lst = [
        # SMC
        [gs, 0, [xmin[0], xmax[0]], [], x_lab[0], y_lab, z_lab, aarr[0][0],
            asigma[0][0], conc_p[0], [], marr[0][0], rad_pc[0], 'SMC'],
        [gs, 1, [xmin[1], xmax[1]], [], x_lab[1], y_lab, z_lab, zarr[0][0],
            zsigma[0][0], conc_p[0], [], marr[0][0], rad_pc[0], 'SMC'],
        # LMC
        [gs, 2, [xmin[0], xmax[0]], [], x_lab[0], y_lab, z_lab, aarr[1][0],
            asigma[1][0], conc_p[1], [], marr[1][0], rad_pc[1], 'LMC'],
        [gs, 3, [xmin[1], xmax[1]], [], x_lab[1], y_lab, z_lab, zarr[1][0],
            zsigma[1][0], conc_p[1], [], marr[1][0], rad_pc[1], 'LMC']
    ]

    for pl_params in conc_pl_lst:
        wide_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/concent_param.png', dpi=300, bbox_inches='tight')


def make_radius_plot(in_params):
    '''
    Plot radius (in pc) versus several parameters.
    '''

    zarr, zsigma, aarr, asigma, marr, msigma, rad_pc, erad_pc, r_core_pc,\
        e_r_core = [in_params[_] for _ in
                    ['zarr', 'zsigma', 'aarr', 'asigma', 'marr', 'msigma',
                    'rad_pc', 'erad_pc', 'r_core_pc', 'e_r_core']]

    # Define values to pass.
    xmin, xmax = 0., 40.
    ymin, ymax = 6., 10.
    x_lab = ['$R_{cl;\,\mathtt{ASteCA}}\,(pc)$',
             '$R_{core;\,\mathtt{ASteCA}}\,(pc)$']
    y_lab = ['$log(age/yr)_{\mathtt{ASteCA}}$', '$[Fe/H]_{\mathtt{ASteCA}}$',
             '$M\,[M_{\odot}]$']
    z_lab = ['$M\,[M_{\odot}]$', '$log(age/yr)_{\mathtt{ASteCA}}$']

    for i, gal_name in enumerate(['SMC', 'LMC']):

        max_r_core = 24 if i == 1 else 42

        fig = plt.figure(figsize=(16, 25))
        gs = gridspec.GridSpec(4, 1)

        rad_pl_lst = [
            [gs, 0, [xmin, xmax], [ymin, ymax], x_lab[0], y_lab[0], z_lab[0],
                rad_pc[i], erad_pc[i], aarr[i][0], asigma[i][0], marr[i][0],
                rad_pc[i], gal_name],
            [gs, 1, [xmin, xmax], [-2.5, 0.5], x_lab[0], y_lab[1], z_lab[0],
                rad_pc[i], erad_pc[i], zarr[i][0], zsigma[i][0], marr[i][0],
                rad_pc[i], gal_name],
            [gs, 2, [xmin, xmax], [-200, 30000], x_lab[0], y_lab[2], z_lab[1],
                rad_pc[i], erad_pc[i], marr[i][0], msigma[i][0], aarr[i][0],
                rad_pc[i], gal_name],
            [gs, 3, [-0.4, max_r_core], [ymin, ymax], x_lab[1], y_lab[0],
                z_lab[0], r_core_pc[i], e_r_core[i], aarr[i][0], asigma[i][0],
                marr[i][0], rad_pc[i], gal_name]
        ]

        for pl_params in rad_pl_lst:
            wide_plots(pl_params)

        # Output png file.
        fig.tight_layout()
        plt.savefig('figures/as_rad_vs_params_' + gal_name + '.png', dpi=300,
                    bbox_inches='tight')


def prob_vs_CI_plot(pl_params):
    '''
    Generate plots for KDE probabilities versus contamination indexes.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, xarr, yarr, zarr, \
        rad, gal_name = pl_params
    siz = np.asarray(rad) * 6

    xy_font_s = 16
    cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i])
    # ax.set_aspect('auto')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    # Plot all clusters in dictionary.
    SC = plt.scatter(xarr, yarr, marker='o', c=zarr, s=siz, lw=0.25, cmap=cm,
                     zorder=3)
    if gal_name != '':
        # Text box.
        ob = offsetbox.AnchoredText(gal_name, loc=2, prop=dict(size=xy_font_s))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    zpad = 10 if z_lab == '$E_{(B-V)}$' else 5
    cbar.set_label(z_lab, fontsize=xy_font_s, labelpad=zpad)


def make_probs_CI_plot(in_params):
    '''
    Plot cluster's ASteCA probabilities versus contamination indexes.
    '''

    zarr, zsigma, aarr, asigma, marr, msigma, rad_pc, kde_prob, cont_ind,\
        n_memb, gal_names = \
        [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'marr',
                                'msigma', 'rad_pc', 'kde_prob', 'cont_ind',
                                'n_memb', 'gal_names']]

    print '* Fraction of clusters with probability < 0.5:'
    for i, gal in enumerate(['SMC', 'LMC']):
        print '  ', gal, float(sum(_ < 0.5 for _ in kde_prob[i])) / \
            float(len(kde_prob[i])), '({})'.format(sum(_ < 0.5 for _ in
                                                       kde_prob[i]))
    print '* Fraction of clusters with probability < 0.25:'
    for i, gal in enumerate(['SMC', 'LMC']):
        print '  ', gal, float(sum(_ < 0.25 for _ in kde_prob[i])) / \
            float(len(kde_prob[i])), '({})'.format(sum(_ < 0.25 for _ in
                                                       kde_prob[i]))

    print '\n* Clusters with n_memb >= 50 & prob =< 0.25'
    for k, gal in enumerate(['SMC', 'LMC']):
        for i, n_m in enumerate(n_memb[k]):
            if n_m >= 50 and kde_prob[k][i] <= 0.25:
                print '', gal, gal_names[k][i], n_m, kde_prob[k][i]

    # Define names of arrays being plotted.
    x_lab, y_lab, z_lab = '$CI_{\mathtt{ASteCA}}$',\
        '$prob_{\mathtt{ASteCA}}$',\
        ['$log(age/yr)_{\mathtt{ASteCA}}$', '$[Fe/H]_{\mathtt{ASteCA}}$',
         '$M\,[M_{\odot}]$', '$M\,[M_{\odot}]$',
         '$log(age/yr)_{\mathtt{ASteCA}}$']
    xmin, xmax, ymin, ymax = -0.01, 1.02, -0.01, 1.02

    fig = plt.figure(figsize=(16, 25))
    gs = gridspec.GridSpec(4, 2)

    prob_CI_pl_lst = [
        # SMC
        [gs, 0, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab[0], cont_ind[0],
            kde_prob[0], aarr[0][0], rad_pc[0], 'SMC'],
        [gs, 1, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab[1], cont_ind[0],
            kde_prob[0], zarr[0][0], rad_pc[0], ''],
        # LMC
        [gs, 2, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab[0], cont_ind[1],
            kde_prob[1], aarr[1][0], rad_pc[1], 'LMC'],
        [gs, 3, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab[1], cont_ind[1],
            kde_prob[1], zarr[1][0], rad_pc[1], ''],
        # Memb number plots
        [gs, 4, -9, 1000, ymin, ymax, '$N_{\mathtt{ASteCA}}$', y_lab, z_lab[0],
         n_memb[0], kde_prob[0], aarr[0][0], rad_pc[0], 'SMC'],
        [gs, 5, -9, 1000, ymin, ymax, '$N_{\mathtt{ASteCA}}$', y_lab, z_lab[0],
         n_memb[1], kde_prob[1], aarr[1][0], rad_pc[1], 'LMC']
    ]

    for pl_params in prob_CI_pl_lst:
        prob_vs_CI_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_prob_vs_CI.png', dpi=300, bbox_inches='tight')


def cross_match_ip_plot_age(pl_params):
    '''
    Generate plots for the cross-matched age values.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, indexes, labels, \
        mark, cols, databases = pl_params

    xy_font_s = 21
    ax = plt.subplot(gs[i])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    ax.tick_params(labelsize=15)

    a, e_a, b, e_b = indexes

    # Plot all clusters for each DB.
    for j, DB in enumerate(databases):
        xarr, yarr = DB[a], DB[b]
        xsigma, ysigma = DB[e_a], DB[e_b]

        db_lab = labels[j] + '$\;(N={})$'.format(len(xarr))
        # Star marker is too small compared to the rest.
        siz = 90. if mark[j] != '*' else 120.
        plt.scatter(xarr, yarr, marker=mark[j], c=cols[j], s=siz,
                    lw=0.3, edgecolor='w', label=db_lab, zorder=3)
        plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')  # 1:1 line
        leg = plt.legend(loc='upper left', markerscale=1.,
                         scatterpoints=1, fontsize=xy_font_s - 7)
        leg.get_frame().set_alpha(0.5)
        # Plot error bars.
        if xsigma:
            for k, xy in enumerate(zip(*[xarr, yarr])):
                if xsigma[k] > 0.:
                    plt.errorbar(xy[0], xy[1], xerr=xsigma[k], ls='none',
                                 color='k', elinewidth=0.2, zorder=1)
        if ysigma:
            for k, xy in enumerate(zip(*[xarr, yarr])):
                if ysigma[k] > 0:
                    plt.errorbar(xy[0], xy[1], yerr=ysigma[k], ls='none',
                                 color='k', elinewidth=0.2, zorder=1)


def make_cross_match_ip_age(cross_match):
    '''
    Plot ASteCA ages versus the values found in the integrated
    photometry databases.
    '''
    # Unpack databases.
    h03, r05, p12 = cross_match[2], cross_match[3], cross_match[6]

    # ASteCA - DB ages
    print '\nDelta (ASteCA - DB) for ages: mean +- std'
    db_name = ['P99', 'P00', 'H03', 'R05', 'C06', 'G10', 'P12']
    for i, db in enumerate(cross_match):
        diff_mean = np.mean(np.array(db[4]) - np.array(db[2]))
        diff_std = np.std(np.array(db[4]) - np.array(db[2]))
        print '{} Delta diffs ages: {:.2f} +- {:.2f}'.format(
            db_name[i], diff_mean, diff_std)

    # First set is for the ages, second for the masses.
    indexes = [4, 5, 2, 3]

    # Labels for each defined plot.
    labels_age = ['H03', 'R05', 'P12']
    mark_age = ['v', '*', 'o']
    cols_age = ['m', 'k', 'b']

    # Define names of arrays being plotted.
    x_lab = '$\log(age/yr)_{\mathtt{ASteCA}}$'
    y_lab = '$\log(age/yr)_{DB}$'

    # Arbitrary size so plots are actually squared.
    fig = plt.figure(figsize=(17.3, 6.3))
    gs = gridspec.GridSpec(1, 3)

    cross_match_lst = [
        # Age cross-match, integrated photometry.
        [gs, 0, 5.8, 10.6, 5.8, 10.6, x_lab, y_lab,
            indexes, labels_age, mark_age, cols_age, [h03, r05, p12]]
    ]

    for pl_params in cross_match_lst:
        cross_match_ip_plot_age(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/cross_match_ip_ages.png', dpi=300,
                bbox_inches='tight')


def cross_match_ip_plot_mass(pl_params):
    '''
    Generate plots for the cross-matched mass values.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, \
        mark, cols, text_box, databases, comb_delta = pl_params

    xy_font_s = 21
    cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    ax.tick_params(labelsize=15)

    a, e_a, b, e_b, s_i, ba_i = 0, 1, 2, 3, 4, 5

    # Plot all clusters for each DB.
    for j, DB in enumerate(databases):
        xarr, yarr = DB[a], DB[b]
        xsigma, ysigma = DB[e_a], DB[e_b]

        # 0 line
        plt.axhline(y=comb_delta[0], xmin=0, xmax=1, color='k', ls='--',
                    lw=0.85)
        # Shaded one sigma region.
        plt.axhspan(comb_delta[0] - comb_delta[1], comb_delta[0] +
                    comb_delta[1], facecolor='grey', alpha=0.2, zorder=1)

        siz = np.array(DB[s_i])*15.
        SC = plt.scatter(xarr, yarr, marker=mark[j], c=DB[ba_i],
                         s=siz, cmap=cm, vmin=-1., vmax=1.5, lw=0.3,
                         edgecolor='k', zorder=3)
        # Plot error bars.
        # x axis error
        for k, xy in enumerate(zip(*[xarr, yarr])):
            if xsigma[k] > 0.:
                plt.errorbar(xy[0], xy[1], xerr=xsigma[k], ls='none',
                             color='k', elinewidth=0.2, zorder=1)
        # y axis error
        for k, xy in enumerate(zip(*[xarr, yarr])):
            if ysigma[k] > 0:
                plt.errorbar(xy[0], xy[1], yerr=ysigma[k], ls='none',
                             color='k', elinewidth=0.2, zorder=1)
    if text_box:
        # Text box.
        ob = offsetbox.AnchoredText(text_box, loc=1,
                                    prop=dict(size=xy_font_s - 5))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    if i == 2:
        # Text box.
        t = r'$\;\;NGC\;419 \,\to\,$' + '\n' + r'$\sim(39, {-}0.87)$'
        ob = offsetbox.AnchoredText(t, loc=4, prop=dict(size=xy_font_s - 7))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
        # Position colorbar.
        axColor = plt.axes([0.885, 0.75, 0.1, 0.023])
        cbar = plt.colorbar(SC, cax=axColor, orientation="horizontal")
        cbar.set_label(r'$\Delta \log(age/yr)$', fontsize=xy_font_s - 3,
                       labelpad=-55)
        cbar.set_ticks([-1., 0., 1.])
        cbar.ax.tick_params(labelsize=xy_font_s - 10)


def make_cross_match_ip_mass(cross_match):
    '''
    Plot ASteCA masses versus the values found in the integrated
    photometry databases.
    '''
    # Unpack databases.
    h03, p12 = cross_match[2], cross_match[6]

    # Calculate correlation coefficients.
    h03_p12_comb_mass = np.array(h03[8] + p12[8])
    ast_h03_p12_comb_mass = np.array(h03[10] + p12[10])
    corr_as = stats.pearsonr(
        ast_h03_p12_comb_mass, ast_h03_p12_comb_mass - h03_p12_comb_mass)
    corr_h03_p12 = stats.pearsonr(
        h03_p12_comb_mass, ast_h03_p12_comb_mass - h03_p12_comb_mass)
    print 'Correlation ASteCA vs Delta mass: {:0.2f}'.format(corr_as[0])
    print 'Correlation H03+P12 vs Delta mass: {:0.2f}'.format(corr_h03_p12[0])

    # Separate clusters with mass < m_limit
    m_low, m_med = 5000., 20000.
    h03_l_mass, p12_l_mass, h03_m_mass, p12_m_mass, h03_h_mass, p12_h_mass =\
        [], [], [], [], [], []
    for cl in zip(*h03):
        # Filter small masses in DBs.
        if cl[8] <= m_low:
            h03_l_mass.append(cl)
        elif m_low < cl[8] <= m_med:
            h03_m_mass.append(cl)
        elif m_med < cl[8]:
            h03_h_mass.append(cl)
            print 'Large mass H03 OC:', cl
    h03_l_mass = zip(*h03_l_mass)
    h03_m_mass = zip(*h03_m_mass)
    h03_h_mass = zip(*h03_h_mass)
    for cl in zip(*p12):
        if cl[8] <= m_low:
            p12_l_mass.append(cl)
        elif m_low < cl[8] <= m_med:
            p12_m_mass.append(cl)
        elif m_med < cl[8]:
            p12_h_mass.append(cl)
            print 'Large mass P12 OC:', cl
    p12_l_mass = zip(*p12_l_mass)
    p12_m_mass = zip(*p12_m_mass)
    p12_h_mass = zip(*p12_h_mass)

    # Low ASteCA - DB masses
    print '\nDelta (ASteCA - DB) for mass_DB<5000: mean +- std'
    db_name = ['H03', 'P12']
    for i, low_m_db in enumerate([h03_l_mass, p12_l_mass]):
        diff_mean = np.mean(np.array(low_m_db[10]) - np.array(low_m_db[8]))
        diff_std = np.std(np.array(low_m_db[10]) - np.array(low_m_db[8]))
        print '{} Delta diffs small mass: {:.0f} +- {:.0f}'.format(
            db_name[i], diff_mean, diff_std)
    comb_delta_l = (np.array(h03_l_mass[10] + p12_l_mass[10]) -
                    np.array(h03_l_mass[8] + p12_l_mass[8]))
    print ("Combined Delta(ASteCA-DB) mean +- std; mass_DB<5000:"
           " {:.2f} +- {:.2f}".format(np.mean(comb_delta_l),
                                      np.std(comb_delta_l)))
    # Relative difference for combined H03-P12 sample.
    comb_l_mass = (np.array(h03_l_mass[10] + p12_l_mass[10]) -
                   np.array(h03_l_mass[8] + p12_l_mass[8])) /\
        (np.array(h03_l_mass[10] + p12_l_mass[10]) +
         np.array(h03_l_mass[8] + p12_l_mass[8]))
    comb_mean, comb_std = np.mean(comb_l_mass), np.std(comb_l_mass)
    comb_rel_delta_l = [comb_mean, comb_std]
    print ("Combined relative Delta mean +- std; mass_DB<5000:"
           " {:.2f} +- {:.2f}".format(comb_mean, comb_std))

    # Medium ASteCA - DB masses
    print '\nDelta (ASteCA - DB) for 5000<mass_DB<20000: mean +- std'
    for i, med_m_db in enumerate([h03_m_mass, p12_m_mass]):
        diff_mean = np.mean(np.array(med_m_db[10]) - np.array(med_m_db[8]))
        diff_std = np.std(np.array(med_m_db[10]) - np.array(med_m_db[8]))
        print '{} Delta diffs medium mass: {:.0f} +- {:.0f}'.format(
            db_name[i], diff_mean, diff_std)
    comb_delta_m = (np.array(h03_m_mass[10] + p12_m_mass[10]) -
                    np.array(h03_m_mass[8] + p12_m_mass[8]))
    print ("Combined Delta(ASteCA-DB) mean +- std; 5000<mass_DB<20000:"
           " {:.2f} +- {:.2f}".format(np.mean(comb_delta_m),
                                      np.std(comb_delta_m)))
    # Relative difference for combined H03-P12 sample.
    comb_m_mass = (np.array(h03_m_mass[10] + p12_m_mass[10]) -
                   np.array(h03_m_mass[8] + p12_m_mass[8])) /\
        (np.array(h03_m_mass[10] + p12_m_mass[10]) +
         np.array(h03_m_mass[8] + p12_m_mass[8]))
    comb_mean, comb_std = np.mean(comb_m_mass), np.std(comb_m_mass)
    comb_rel_delta_m = [comb_mean, comb_std]
    print ("Combined relative Delta mean +- std; 5000<mass_DB<20000:"
           " {:.2f} +- {:.2f}".format(comb_mean, comb_std))

    # Large ASteCA - DB masses
    print '\nDelta (ASteCA - DB) for mass_DB>20000: mean +- std'
    print 'H03 OCs:', len(h03_h_mass[0]), 'P12 OCs:', len(p12_h_mass[0])
    for i, h_m_db in enumerate([h03_h_mass, p12_h_mass]):
        diff_mean = np.mean(np.array(h_m_db[10]) - np.array(h_m_db[8]))
        diff_std = np.std(np.array(h_m_db[10]) - np.array(h_m_db[8]))
        print '{} Delta diffs large mass: {:.0f} +- {:.0f}'.format(
            db_name[i], diff_mean, diff_std)
    comb_delta_h = (np.array(h03_h_mass[10] + p12_h_mass[10]) -
                    np.array(h03_h_mass[8] + p12_h_mass[8]))
    print ("Combined Delta(ASteCA-DB) mean +- std; mass_DB>20000:"
           " {:.2f} +- {:.2f}".format(np.mean(comb_delta_h),
                                      np.std(comb_delta_h)))
    # Relative difference for combined H03-P12 sample.
    comb_h_mass = (np.array(h03_h_mass[10] + p12_h_mass[10]) -
                   np.array(h03_h_mass[8] + p12_h_mass[8])) /\
        (np.array(h03_h_mass[10] + p12_h_mass[10]) +
         np.array(h03_h_mass[8] + p12_h_mass[8]))
    comb_mean, comb_std = np.mean(comb_h_mass), np.std(comb_h_mass)
    comb_rel_delta_h = [comb_mean, comb_std]
    print ("Combined mean +- std for Delta(ASteCA-DB) mass_DB>20000:"
           " {:.2f} +- {:.2f}".format(comb_mean, comb_std))

    # Low masses relative differences ASteCA  - DBs
    h03_mass_diff_l = (np.array(h03_l_mass[10]) - np.array(h03_l_mass[8])) /\
        (np.array(h03_l_mass[10]) + np.array(h03_l_mass[8]))
    p12_mass_diff_l = (np.array(p12_l_mass[10]) - np.array(p12_l_mass[8])) /\
        (np.array(p12_l_mass[10]) + np.array(p12_l_mass[8]))
    # Errors.
    A, s_A, B, s_B = np.array(h03_l_mass[10]), np.array(h03_l_mass[11]),\
        np.array(h03_l_mass[8]), np.array(h03_l_mass[9])
    h03_delta_err_l = list((2./(A+B)**2) * np.sqrt((B*s_A)**2 + (A*s_B)**2))
    A, s_A, B, s_B = np.array(p12_l_mass[10]), np.array(p12_l_mass[11]),\
        np.array(p12_l_mass[8]), np.array(p12_l_mass[9])
    p12_delta_err_l = list((2./(A+B)**2) * np.sqrt((B*s_A)**2 + (A*s_B)**2))
    # Medium masses relative differences ASteCA  - DBs
    h03_mass_diff_m = (np.array(h03_m_mass[10]) - np.array(h03_m_mass[8])) /\
        (np.array(h03_m_mass[10]) + np.array(h03_m_mass[8]))
    p12_mass_diff_m = (np.array(p12_m_mass[10]) - np.array(p12_m_mass[8])) /\
        (np.array(p12_m_mass[10]) + np.array(p12_m_mass[8]))
    # Errors.
    A, s_A, B, s_B = np.array(h03_m_mass[10]), np.array(h03_m_mass[11]),\
        np.array(h03_m_mass[8]), np.array(h03_m_mass[9])
    h03_delta_err_m = list((2./(A+B)**2) * np.sqrt((B*s_A)**2 + (A*s_B)**2))
    A, s_A, B, s_B = np.array(p12_m_mass[10]), np.array(p12_m_mass[11]),\
        np.array(p12_m_mass[8]), np.array(p12_m_mass[9])
    p12_delta_err_m = list((2./(A+B)**2) * np.sqrt((B*s_A)**2 + (A*s_B)**2))
    # Large masses relative differences ASteCA  - DBs
    h03_mass_diff_h = (np.array(h03_h_mass[10]) - np.array(h03_h_mass[8])) /\
        (np.array(h03_h_mass[10]) + np.array(h03_h_mass[8]))
    p12_mass_diff_h = (np.array(p12_h_mass[10]) - np.array(p12_h_mass[8])) /\
        (np.array(p12_h_mass[10]) + np.array(p12_h_mass[8]))
    # Errors
    A, s_A, B, s_B = np.array(h03_h_mass[10]), np.array(h03_h_mass[11]),\
        np.array(h03_h_mass[8]), np.array(h03_h_mass[9])
    h03_delta_err_h = list((2./(A+B)**2) * np.sqrt((B*s_A)**2 + (A*s_B)**2))
    A, s_A, B, s_B = np.array(p12_h_mass[10]), np.array(p12_h_mass[11]),\
        np.array(p12_h_mass[8]), np.array(p12_h_mass[9])
    p12_delta_err_h = list((2./(A+B)**2) * np.sqrt((B*s_A)**2 + (A*s_B)**2))

    # Low mass ASteCA-DBs. Use age deltas for Colors.
    colors_h03_l, colors_p12_l = np.array(h03_l_mass[4]) -\
        np.array(h03_l_mass[2]), np.array(p12_l_mass[4]) -\
        np.array(p12_l_mass[2])
    # Use this for averaged values in the x axis. Replace these lists in the 0th
    # positions of both sub-lists, in the delta_DBs_X lists.
    h03_avr_l_mass = (np.array(h03_l_mass[10]) + np.array(h03_l_mass[8])) / 2.
    p12_avr_l_mass = (np.array(p12_l_mass[10]) + np.array(p12_l_mass[8])) / 2.
    delta_DBs_l = [[h03_l_mass[8], h03_l_mass[9], h03_mass_diff_l,
                    h03_delta_err_l, h03_l_mass[19], colors_h03_l],
                   [p12_l_mass[8], p12_l_mass[9], p12_mass_diff_l,
                    p12_delta_err_l, p12_l_mass[19], colors_p12_l]]
    # Medium mass ASteCA-DBs.
    colors_h03_m, colors_p12_m = np.array(h03_m_mass[4]) -\
        np.array(h03_m_mass[2]), np.array(p12_m_mass[4]) -\
        np.array(p12_m_mass[2])
    scale = 10.**4
    # Averaged masses in x axis.
    # h03_avr_m_mass = (np.array(h03_m_mass[10]) + np.array(h03_m_mass[8])) / 2.
    # p12_avr_m_mass = (np.array(p12_m_mass[10]) + np.array(p12_m_mass[8])) / 2.
    delta_DBs_m = [[np.array(h03_m_mass[8]) / scale,
                    np.array(h03_m_mass[9]) / scale, h03_mass_diff_m,
                    h03_delta_err_m, h03_m_mass[19], colors_h03_m],
                   [np.array(p12_m_mass[8]) / scale,
                    np.array(p12_m_mass[9]) / scale, p12_mass_diff_m,
                    p12_delta_err_m, p12_m_mass[19], colors_p12_m]]
    # Large mass ASteCA-DBs.
    colors_h03_h, colors_p12_h = np.array(h03_h_mass[4]) -\
        np.array(h03_h_mass[2]), np.array(p12_h_mass[4]) -\
        np.array(p12_h_mass[2])
    # Averaged masses in x axis.
    # h03_avr_h_mass = (np.array(h03_h_mass[10]) + np.array(h03_h_mass[8])) / 2.
    # p12_avr_h_mass = (np.array(p12_h_mass[10]) + np.array(p12_h_mass[8])) / 2.
    delta_DBs_h = [[np.array(h03_h_mass[8]) / scale,
                    np.array(h03_h_mass[9]) / scale, h03_mass_diff_h,
                    h03_delta_err_h, h03_h_mass[19], colors_h03_h],
                   [np.array(p12_h_mass[8]) / scale,
                    np.array(p12_h_mass[9]) / scale, p12_mass_diff_h,
                    p12_delta_err_h, p12_h_mass[19], colors_p12_h]]

    # Define data to pass.
    databases = [delta_DBs_l, delta_DBs_m, delta_DBs_h]

    # Labels for each defined plot.
    mark_mass = ['v', 'o']
    cols_mass = ['m', 'b']

    # Define names of arrays being plotted.
    x_lab = ['$M_{DBs}\,[M_{\odot}]$', '$M_{DBs}\,[10^{-4}M_{\odot}]$']
    y_lab = [r'$\overline{\Delta M_r}\;\;(\mathtt{ASteCA}-DBs)$', '']
    # Limits when using DBs masses in x axis.
    l_mass_lims = [-50., 4990., -1.19, 1.19]
    m_mass_lims = [0.5, 2.05, -1.19, 0.45]
    h_mass_lims = [2.05, 10.6, -1.02, -0.51]
    # Limits when using averaged masses in x axis.
    # l_mass_lims = [-50., 5990., -1.19, 1.19]
    # m_mass_lims = [0.3, 3.05, -1.19, 0.45]
    # h_mass_lims = [1.05, 10.6, -1.02, -0.51]

    # Arbitrary size so plots are actually squared.
    fig = plt.figure(figsize=(19.3, 6.3))
    gs = gridspec.GridSpec(1, 3)

    cross_match_lst = [
        # Mass cross_match (low mass)
        [gs, 0, l_mass_lims[0], l_mass_lims[1], l_mass_lims[2], l_mass_lims[3],
         x_lab[0], y_lab[0], mark_mass, cols_mass,
         '$M_{DBs}\leq 5000\,[M_{\odot}]$', databases[0], comb_rel_delta_l],
        # Mass cross_match (medium masses)
        [gs, 1, m_mass_lims[0], m_mass_lims[1], m_mass_lims[2], m_mass_lims[3],
         x_lab[1], y_lab[1], mark_mass, cols_mass,
         '$5000<M_{DBs}\leq 20000\,[M_{\odot}]$', databases[1],
         comb_rel_delta_m],
        # Mass cross_match (large masses)
        [gs, 2, h_mass_lims[0], h_mass_lims[1], h_mass_lims[2], h_mass_lims[3],
         x_lab[1], y_lab[1], mark_mass, cols_mass,
         '$M_{DBs}>20000\,[M_{\odot}]$', databases[2], comb_rel_delta_h]
    ]

    for pl_params in cross_match_lst:
        cross_match_ip_plot_mass(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/cross_match_ip_mass.png', dpi=300,
                bbox_inches='tight')


def cross_match_if_plot(pl_params):
    '''
    Generate plots for the cross-matched isochrone fitted OCs.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, data, labels, mark, cols, \
        kde_cont = pl_params

    xy_font_s = 21
    ax = plt.subplot(gs[i])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    ax.tick_params(labelsize=15)
    if i in [1, 2]:
        # Origin lines.
        plt.plot([-10, 10], [0., 0.], 'k', ls='--')
        plt.plot([0., 0.], [-10, 10], 'k', ls='--')
    else:
        # 1:1 line
        plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')
    # Plot all clusters for each DB.
    for j, DB in enumerate(data):
        xarr, yarr = DB[0], DB[2]
        xsigma, ysigma = DB[1], DB[3]
        siz = 60. if mark[j] != '*' else 80.

        if i in [0, 2]:
            db_lab = labels[j] + '$\;(N={})$'.format(len(xarr))
        else:
            db_lab = labels[j]
        plt.scatter(xarr, yarr, marker=mark[j], c=cols[j], s=siz,
                    lw=0.25, edgecolor='w', label=db_lab, zorder=3)
        # Plot error bars.
        if xsigma:
            for k, xy in enumerate(zip(*[xarr, yarr])):
                x_err = xsigma[k] if 0. < xsigma[k] < 5. else 0.
                plt.errorbar(xy[0], xy[1], xerr=x_err,
                             ls='none', color='k', elinewidth=0.2, zorder=1)
        if ysigma:
            for k, xy in enumerate(zip(*[xarr, yarr])):
                y_err = ysigma[k] if 0. < ysigma[k] < 5. else 0.
                plt.errorbar(xy[0], xy[1], yerr=y_err,
                             ls='none', color='k', elinewidth=0.2, zorder=1)
        # Plot KDE.
        if kde_cont:
            x, y, kde = kde_cont
            # plt.imshow(np.rot90(kde), cmap=plt.cm.YlOrBr, extent=ext_range)
            plt.contour(x, y, kde, 5, colors='k', linewidths=0.6)

    # Legend.
    leg = plt.legend(loc='upper left', markerscale=1., scatterpoints=1,
                     fontsize=xy_font_s - 7)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.5)
    ax.set_aspect('auto')


def make_cross_match_if(cross_match, in_params):
    '''
    Plot the differences between extinction and age for ASteCA values versus
    Washington values (ie: Piatti et al. values) and ASteCA values versus
    the databases where the isochrone fitting method was used.
    '''
    # unpack databases.
    p99, p00, h03, r05, c06, g10, p12 = cross_match

    aarr, asigma, earr, esigma = \
        [in_params[_] for _ in ['aarr', 'asigma', 'earr', 'esigma']]

    # Define lists of ASteCA minus literature values.
    # SMC ASteCA minus literature diffs.
    diffs_lit_ages_smc = np.array(aarr[0][0]) - np.array(aarr[0][1])
    diffs_lit_exts_smc = np.array(earr[0][0]) - np.array(earr[0][1])
    # LMC ASteCA minus literature diffs.
    diffs_lit_ages_lmc, diffs_lit_exts_lmc = [], []
    for i, lit_ext in enumerate(earr[1][1]):
        # Remove 99.9 values from 'M' reference that contains no extinction
        # estimates.
        if abs(lit_ext) < 5:
            diffs_lit_ages_lmc.append(aarr[1][0][i] - aarr[1][1][i])
            diffs_lit_exts_lmc.append(earr[1][0][i] - earr[1][1][i])

    # Define lists of difference between ages and extinctions.
    # Age indexes (DB, ASteCA, lit) -> 2, 4, 6
    # Ext indexes (DB, ASteCA, lit) -> 14, 15, 17

    # P99 ASteCA minus database diffs.
    diffs_db_ages_p99 = np.array(p99[4]) - np.array(p99[2])
    # Same for extinctions.
    diffs_db_exts_p99 = np.array(p99[15]) - np.array(p99[14])
    # P00 ASteCA minus database diffs (no extinction information).
    diffs_db_ages_p00 = np.array(p00[4]) - np.array(p00[2])
    diffs_db_exts_p00 = np.array(p00[15]) - np.array(p00[14])
    # C06
    diffs_db_ages_c06 = np.array(c06[4]) - np.array(c06[2])
    diffs_db_exts_c06 = np.array(c06[15]) - np.array(c06[14])
    # G10
    diffs_db_ages_g10 = np.array(g10[4]) - np.array(g10[2])
    diffs_db_exts_g10 = np.array(g10[15]) - np.array(g10[14])

    # Calculate std, means and medians for the age differences.
    txt = ['SMC', 'LMC', 'P99', 'P00', 'C06', 'G10']
    dbs = [diffs_lit_ages_smc, diffs_lit_ages_lmc, diffs_db_ages_p99,
           diffs_db_ages_p00, diffs_db_ages_c06, diffs_db_ages_g10]
    for i, db in enumerate(dbs):
        print "{}, diff ages = {:.3f} +- {:.3f}".format(
                txt[i], np.mean(db), np.std(db))

    # Obtain a Gaussian KDE for each plot.
    # Define x,y grid.
    gd_c = complex(0, 100)
    kde_cont = []
    for xarr, yarr in [
            [list(diffs_db_ages_p99) + list(diffs_db_ages_p00) +
             list(diffs_db_ages_c06) + list(diffs_db_ages_g10),
             list(diffs_db_exts_p99) + list(diffs_db_exts_p00) +
             list(diffs_db_exts_c06) + list(diffs_db_exts_g10)],
            [list(diffs_lit_ages_smc) + list(diffs_lit_ages_lmc),
             list(diffs_lit_exts_smc) + list(diffs_lit_exts_lmc)]
    ]:
        values = np.vstack([xarr, yarr])
        kernel = stats.gaussian_kde(values)
        xmin, xmax, ymin, ymax = -2., 2., -1., 1.
        x, y = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
        positions = np.vstack([x.ravel(), y.ravel()])
        # Evaluate kernel in grid positions.
        k_pos = kernel(positions)
        kde = np.reshape(k_pos.T, x.shape)
        kde_cont.append([x, y, kde])

    # Order data to plot.
    # Extinction vs ages differences.
    lit_data = [[diffs_lit_ages_smc, [], diffs_lit_exts_smc, []],
                [diffs_lit_ages_lmc, [], diffs_lit_exts_lmc, []]]
    db_data = [[diffs_db_ages_p99, [], diffs_db_exts_p99, []],
               [diffs_db_ages_p00, [], diffs_db_exts_p00, []],
               [diffs_db_ages_c06, [], diffs_db_exts_c06, []],
               [diffs_db_ages_g10, [], diffs_db_exts_g10, []]]
    age_ast_DB_data = [[p99[4], p99[5], p99[2], p99[3]],
                       [p00[4], p00[5], p00[2], p00[3]],
                       [c06[4], c06[5], c06[2], c06[3]],
                       [g10[4], g10[5], g10[2], g10[3]]]

    labels = [['P99', 'P00', 'C06', 'G10'], ['SMC', 'LMC']]
    mark = [['>', '^', 'v', '<'], ['*', '*']]
    cols = [['chocolate', 'r', 'c', 'g'], ['m', 'b']]

    # Define names of arrays being plotted.
    x_lab = ['$\Delta \log(age/yr)_{\mathtt{ASteCA}-DB}$',
             '$\Delta \log(age/yr)_{\mathtt{ASteCA}-lit}$',
             '$\log(age/yr)_{\mathtt{ASteCA}}$']
    y_lab = ['$\Delta E(B-V)_{\mathtt{ASteCA}-DB}$',
             '$\Delta E(B-V)_{\mathtt{ASteCA}-lit}$',
             '$\log(age/yr)_{DB}$']
    xmm, ymm = [-1.5, 1.5, 5.8], [-0.19, 0.19, 10.6]

    # Arbitrary size so plots are actually squared.
    fig = plt.figure(figsize=(20, 6.3))
    gs = gridspec.GridSpec(1, 3)

    cross_match_lst = [
        # Age 1:1, isoch fit lit vs ASteCA.
        [gs, 0, xmm[2], ymm[2], xmm[2], ymm[2], x_lab[2], y_lab[2],
            age_ast_DB_data, labels[0], mark[0], cols[0], []],
        # Age vs ext diff for ASteCA vs databases.
        [gs, 1, xmm[0], xmm[1], ymm[0], ymm[1], x_lab[0], y_lab[0],
            db_data, labels[0], mark[0], cols[0], kde_cont[0]],
        # Age vs ext diff for ASteCA vs literature.
        [gs, 2, xmm[0], xmm[1], ymm[0], ymm[1], x_lab[1], y_lab[1],
            lit_data, labels[1], mark[1], cols[1], kde_cont[1]]
    ]

    for pl_params in cross_match_lst:
        cross_match_if_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/cross_match_if.png', dpi=300,
                bbox_inches='tight')


def pl_DBs_ASteCA_CMDs(pl_params):
    '''
    Star's membership probabilities on cluster's photom diagram.
    '''
    gs, i, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax, cl, db,\
        gal, cl_reg_fit, cl_reg_no_fit, synth_stars, lit_isoch, asteca_isoch,\
        db_z, db_a, db_e, db_d, as_z, as_a, as_e, as_d, as_m = pl_params

    letter = ['a', '', 'b', '', '', '', 'c', '', 'd', '', '', '', 'e', '',
              'f', '', '', '', 'g', '', 'h', '', '', '', 'i', '', 'j', '']

    # DB/outlier isoch fit.
    ax0 = plt.subplot(gs[i])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Add text box.
    ob0 = offsetbox.AnchoredText(letter[i], loc=2, prop=dict(size=12))
    ob0.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax0.add_artist(ob0)
    if db == 'outliers':
        t = 'Lit'
    elif db == 'largemass':
        t = 'Synthetic'
    else:
        t = db
    text = gal + '-' + cl + ' ({})'.format(t)
    ob1 = offsetbox.AnchoredText(text, loc=1, prop=dict(size=11))
    ob1.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax0.add_artist(ob1)
    text1 = r'$z={}$'.format(db_z)
    text2 = '\n' + r'$log(age/yr)={}$'.format(float(db_a))
    text3 = '\n' + r'$E_{{(B-V)}}={}$'.format(db_e)
    text4 = '\n' + r'$dm={}$'.format(db_d)
    text5 = '\n' + r'$M={}\,M_{{\odot}}$'.format(int(float(as_m)))
    if db != 'largemass':
        text = text1 + text2 + text3 + text4
        locat = 3
    else:
        text = text1 + text2 + text3 + text4 + text5
        locat = 4 if letter[i] in ['b', 'd'] else 3
    ob2 = offsetbox.AnchoredText(text, loc=locat, prop=dict(size=11))
    ob2.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax0.add_artist(ob2)
    # Set minor ticks
    ax0.minorticks_on()
    ax0.xaxis.set_major_locator(MultipleLocator(1.0))
    # Plot grid.
    ax0.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    col_select_fit, c_iso = '#4682b4', 'r'
    if db != 'largemass':
        # Plot stars used in the best fit process.
        cl_reg_x = cl_reg_fit[0] + cl_reg_no_fit[0]
        cl_reg_y = cl_reg_fit[1] + cl_reg_no_fit[1]
    else:
        cl_reg_x, cl_reg_y = synth_stars
    plt.scatter(cl_reg_x, cl_reg_y, marker='o',
                c=col_select_fit, s=40, cmap=cm, lw=0.5, zorder=4)
    # Plot isochrone.
    plt.plot(lit_isoch[0], lit_isoch[1], c=c_iso, lw=1.2, zorder=5)

    # ASteCA isoch fit.
    ax1 = plt.subplot(gs[i + 1])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('')
    ax1.axes.yaxis.set_ticklabels([])
    # Add text box.
    ob0 = offsetbox.AnchoredText(letter[i], loc=2, prop=dict(size=12))
    ob0.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax1.add_artist(ob0)
    text = gal + '-' + cl + r'$\,(\mathtt{ASteCA})$'
    ob1 = offsetbox.AnchoredText(text, loc=1, prop=dict(size=11))
    ob1.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax1.add_artist(ob1)
    if db != 'largemass':
        text1 = r'$z={}$'.format(as_z)
        text2 = '\n' + r'$log(age/yr)={}$'.format(as_a)
        text3 = '\n' + r'$E_{{(B-V)}}={}$'.format(as_e)
        text4 = '\n' + r'$dm={}$'.format(as_d)
        text5 = '\n' + r'$M={}\,M_{{\odot}}$'.format(int(float(as_m)))
        text = text1 + text2 + text3 + text4 + text5
        ob = offsetbox.AnchoredText(text, loc=3, prop=dict(size=11))
        ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
        ax1.add_artist(ob)
    # Set minor ticks
    ax1.minorticks_on()
    ax1.xaxis.set_major_locator(MultipleLocator(1.0))
    # Plot grid.
    ax1.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Get extreme values for colorbar.
    lst_comb = cl_reg_fit[2] + cl_reg_no_fit[2]
    v_min_mp, v_max_mp = round(min(lst_comb), 2), round(max(lst_comb), 2)
    col_select_fit, col_select_no_fit, c_iso = cl_reg_fit[2], \
        cl_reg_no_fit[2], 'g'
    siz = 40 if db != 'largemass' else 20
    # Plot stars *not* used in the best fit process.
    plt.scatter(cl_reg_no_fit[0], cl_reg_no_fit[1], marker='o',
                c=col_select_no_fit, s=siz-5, cmap=cm, lw=0.5, alpha=0.5,
                vmin=v_min_mp, vmax=v_max_mp, zorder=2)
    # Plot stars used in the best fit process.
    plt.scatter(cl_reg_fit[0], cl_reg_fit[1], marker='o',
                c=col_select_fit, s=siz, cmap=cm, lw=0.5, vmin=v_min_mp,
                vmax=v_max_mp, zorder=4)
    # Plot isochrone.
    plt.plot(asteca_isoch[0], asteca_isoch[1], c=c_iso, lw=1.2, zorder=5)


def make_DB_ASteCA_CMDs(db, db_cls):
    '''
    '''
    for k, cl_lst in enumerate(db_cls):

        fig = plt.figure(figsize=(30, 25))
        gs = gridspec.GridSpec(5, 6)

        i, j, db_sat_cmd_lst = 0, 1, []
        for cl_data in cl_lst:

            x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, cl, db, gal, \
                cl_reg_fit, cl_reg_no_fit, synth_stars, lit_isoch,\
                asteca_isoch, db_z, db_a, db_e, db_d, as_z, as_a, as_e,\
                as_d, as_m = cl_data

            db_sat_cmd_lst.append(
                [gs, i, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, '(C-T_1)',
                    'T_1', cl, db, gal, cl_reg_fit, cl_reg_no_fit, synth_stars,
                    lit_isoch, asteca_isoch, db_z, db_a, db_e, db_d, as_z,
                    as_a, as_e, as_d, as_m])

            # Plotting positions.
            if (j % 2 == 0):  # even
                i += 4
            else:  # odd
                i += 2
            j += 1

        for pl_params in db_sat_cmd_lst:
            pl_DBs_ASteCA_CMDs(pl_params)

        # Output png file.
        fig.tight_layout()
        if db in ['outliers', 'largemass']:
            r_path = 'figures/'
        else:
            r_path = 'figures/DB_fit/'
        fig_name = r_path + db + '_VS_asteca_' + str(k) + '.png'
        plt.savefig(fig_name, dpi=150, bbox_inches='tight')


def pl_errors(pl_params):
    '''
    '''
    gs, i, xmin, xmax, ymin, ymax, x, y, z, rad, x_lab, y_lab =\
        pl_params
    siz = np.asarray(rad) * 1.

    xy_font_s = 16
    ax = plt.subplot(gs[i, 0:4])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    cm = plt.cm.get_cmap('RdYlBu_r')
    # Introduce random scatter.
    if i == 0:
        # 1% of axis ranges.
        ax_ext = (xmax - xmin) * 0.02
    elif i == 4:
        ax_ext = (xmax - xmin) * 0.01
    else:
        ax_ext = (xmax - xmin) * 0.025
    # Add random scatter.
    r_x = x + np.random.uniform(-ax_ext, ax_ext, len(x))
    SC = plt.scatter(r_x, y, marker='o', c=z, edgecolor='k', s=siz,
                     cmap=cm, lw=0.25, zorder=4)  # vmin=0., vmax=1.,
    # Horizontal histogram
    axHisty = plt.subplot(gs[i, 4:5])
    axHisty.minorticks_on()
    axHisty.grid(b=True, which='major', color='white', linestyle='-', lw=0.5,
                 zorder=1)
    axHisty.hist(y, bins=20, orientation='horizontal', normed=1,
                 histtype='stepfilled', color='#22596a', zorder=3)
    axHisty.set_ylim(ax.get_ylim())
    axHisty.set_xticklabels([])
    axHisty.set_yticklabels([])
    axHisty.set_axis_bgcolor('#ced6d8')
    axHisty.axhline(y=np.mean(y), color='r', ls='--', lw=2, zorder=5)
    # Text box.
    if i != 4:
        txt = r'$\bar{{e}}={:.2f}$'.format(np.mean(y))
    else:
        txt = r'$\bar{{e}}={:.0f}$'.format(np.mean(y))
    ob = offsetbox.AnchoredText(txt, loc=1, prop=dict(size=xy_font_s-1))
    ob.patch.set(alpha=0.85)
    axHisty.add_artist(ob)
    print 'Mean {}: {}'.format(y_lab, np.mean(y))

    if i == 0:
        # Position colorbar.
        axColor = plt.axes([0.1, 0.85, 0.3, 0.005])
        cbar = plt.colorbar(SC, cax=axColor, orientation="horizontal")
        cbar.set_label(r'$CI$', fontsize=xy_font_s - 2, labelpad=-37)
        cbar.set_ticks([0.2, 0.4, 0.6, 0.8, 1., 1.2])
        cbar.ax.tick_params(labelsize=xy_font_s - 8)


def make_errors_plots(in_params):
    '''
    '''
    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, marr, msigma,\
        rarr, cont_ind = [
            in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
                                   'esigma', 'darr', 'dsigma', 'marr',
                                   'msigma', 'rarr', 'cont_ind']]
                                   # , 'kde_prob', 'phot_disp']]

    ci = cont_ind[0] + cont_ind[1]
    r_arr = rarr[0][0] + rarr[1][0]
    z_arr = zarr[0][0] + zarr[1][0]
    z_sigma = zsigma[0][0] + zsigma[1][0]
    # # Transform [Fe/H] values to z
    # z_arr = list((10**np.array(z_arr)) * 0.0152)
    # z_sigma = list(np.array(z_sigma) * np.array(z_arr) * np.log(10.))
    #
    a_arr = aarr[0][0] + aarr[1][0]
    a_sigma = asigma[0][0] + asigma[1][0]
    e_arr = earr[0][0] + earr[1][0]
    e_sigma = esigma[0][0] + esigma[1][0]
    d_arr = darr[0][0] + darr[1][0]
    d_sigma = dsigma[0][0] + dsigma[1][0]
    m_arr = marr[0][0] + marr[1][0]
    m_sigma = zsigma[0][0] + msigma[1][0]
    # p_disp = phot_disp[0] + phot_disp[1]
    # probs = kde_prob[0] + kde_prob[1]

    # Order lists to put min rad values on top.
    ord_r, ord_z, ord_zs, ord_a, ord_as, ord_e, ord_es, ord_d, ord_ds, ord_m,\
        ord_ms, ord_ci = map(list, zip(*sorted(zip(
            r_arr, z_arr, z_sigma, a_arr, a_sigma, e_arr, e_sigma, d_arr,
            d_sigma, m_arr, m_sigma, ci), reverse=True)))
    # ord_prob, ord_p_disp

    # Select colorbar parameter.
    # ord_X = np.array(ord_prob) / np.array(ord_p_disp)
    ord_X = ord_ci

    fig = plt.figure(figsize=(10, 23))
    gs = gridspec.GridSpec(5, 5)

    errors_lst = [
        [gs, 0, -2.4, 0.11, -0.03, 2.1, ord_z, ord_zs, ord_X, ord_r,
            '$[Fe/H]$', '$e_{[Fe/H]}$'],
        # [gs, 0, 0., 0.016, 0., 0.01, ord_z, ord_zs, ord_X, ord_r,
        #     '$[Fe/H]$', '$e_{[Fe/H]}$'],
        [gs, 1, 6.51, 10.1, -0.03, 1.1, ord_a, ord_as, ord_X, ord_r,
            r'$\log(aye/yr)$', '$e_{\log(aye/yr)}$'],
        [gs, 2, -0.02, 0.32, -0.01, 0.11, ord_e, ord_es, ord_X, ord_r,
            '$E_{B-V}$', '$e_{E_{B-V}}$'],
        [gs, 3, 18.28, 19.19, 0.007, 0.083, ord_d, ord_ds, ord_X, ord_r,
            '$\mu_{0}$', '$e_{\mu_{0}}$'],
        [gs, 4, -210, 30000, -210, 4450, ord_m, ord_ms, ord_X, ord_r,
            '$M\,[M_{\odot}]$', '$e_{M}$']
    ]

    for pl_params in errors_lst:
        pl_errors(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/errors_asteca.png', dpi=300, bbox_inches='tight')


def pl_amr(pl_params):
    '''
    Plot AMRs.
    '''

    gs, i, age_vals, met_weighted, age_gyr, amr_lit, feh, rad_pc, x_lab,\
        y_lab, ast_lit = pl_params

    xy_font_s = 16
    ax = plt.subplot(gs[i])
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    plt.xlim(-0.02, 8.4)
    l = r'$\,\mathtt{ASteCA}$' if ast_lit == 0 else 'Literature'

    if i == 0:
        plt.ylim(-2.45, 0.4)
        ax.set_xticklabels([])
        col, leg = ['r', 'b'], ['SMC', 'LMC']
        for k in [0, 1]:
            for x_ed in amr_lit[k]:
                # vertical lines
                vl_lim = [-5., -2.2] if k == 0 else [0.2, 1.]
                ax.vlines(x=x_ed, ymin=vl_lim[0], ymax=vl_lim[1],
                          linestyle='-', color=col[k], lw=1.5, zorder=1)
            # Introduce random scatter.
            ax_x, ax_y = 0.15, 0.03
            # Add random scatter.
            rs_x = age_gyr[k][0] + np.random.uniform(-ax_x, ax_x,
                                                     len(age_gyr[k][0]))
            rs_y = feh[k] + np.random.uniform(-ax_y, ax_y, len(feh[k]))
            siz = np.array(rad_pc[k]) * 5.
            plt.scatter(rs_x, rs_y, marker='*', s=siz, edgecolors=col[k],
                        facecolor='none', lw=0.4, label=leg[k] + ' OCs',
                        zorder=4)
            # 1 sigma error regions.
            y_err_min = np.array(met_weighted[k][0]) -\
                np.array(met_weighted[k][1])
            y_err_max = np.array(met_weighted[k][0]) +\
                np.array(met_weighted[k][1])
            plt.fill_between(age_vals[k], y_err_min, y_err_max, alpha=0.1,
                             color=col[k], zorder=5)
            # ASteCA/Literature AMR.
            plt.plot(age_vals[k], met_weighted[k][0], c=col[k], lw=1.2,
                     label=leg[k] + '(' + l + ')', zorder=8)
        # Legend.
        leg0 = plt.legend(loc='upper right',  # bbox_to_anchor=(1., 0.12),
                          handlelength=2.5, scatterpoints=1,
                          fontsize=xy_font_s - 8)
        leg0.legendHandles[0]._sizes = [10]
        leg0.legendHandles[1]._sizes = [10]
        leg0.get_frame().set_alpha(0.85)

    # Literature values.
    elif i == 1:
        gal, k, c_as, leg_cut, l_thick = 'LMC', 1, 'b', 6, 6
        ax.set_xticklabels([])
        ymin, ymax = -1.21, -0.1
        col = ['#ff661a', '#8080ff', 'y', 'g', 'c', 'k', 'm', '#33cc33',
               '#b22222']
        c_dash = [[8, 4], [2, 2], [8, 4], [2, 2], [8, 4, 2, 4, 2, 4],
                  [8, 4, 2, 4], [1000, 1], [8, 4, 2, 4, 2, 4], [8, 4, 2, 4]]
        amr_lab = ['PT98', 'G98', 'C08a', 'HZ09', 'R12', 'PG13', 'M14-0',
                   'M14-1', 'M14-2']
    elif i == 2:
        gal, k, c_as, leg_cut, l_thick = 'SMC', 0, 'r', 7, 5
        ymin, ymax = -1.34, -0.29
        plt.xlabel(x_lab, fontsize=xy_font_s)
        col = ['#ff661a', '#8080ff', 'y', 'c', '#33cc33', 'm', 'g', '#b22222',
               '#b22222', 'k']
        c_dash = [[8, 4], [8, 4, 2, 4], [8, 4], [2, 2], [8, 4, 2, 4, 2, 4],
                  [1000, 1], [8, 4, 2, 4, 2, 4], [8, 4], [2, 2],
                  [8, 4, 2, 4]]
        amr_lab = ['PT98', 'HZ04', 'C08b', 'N09', 'TB09-1', 'TB09-2', 'TB09-3',
                   'C13-B', 'C13-C', 'PG13']
    if i in [1, 2]:
        plt.ylim(ymin, ymax)
        ax.set_title(gal, x=0.5, y=0.92, fontsize=xy_font_s - 4,
                     bbox=dict(facecolor=(1, 1, 1, 0.5),
                               edgecolor=(0, 0, 0, 1)))
        hand1, hand2 = [], []  # Store handles.
        for j, amr in enumerate(amr_lit):
            # The continuous line for the AMR should be thinner.
            l_w = 1.5 if j != l_thick else 0.85
            pl, = plt.plot(amr[0], amr[1], color=col[j], label=amr_lab[j],
                           dashes=c_dash[j], lw=l_w, zorder=3)
            # Using two legends.
            if j < leg_cut:
                hand1.append(pl)
            else:
                hand2.append(pl)
        # ASteCA/Literature values.
        pl, = plt.plot(age_vals[k], met_weighted[k][0], c=c_as, lw=1.2,
                       label=l, zorder=5)
        hand2.append(pl)
        # Legend.
        leg1 = plt.legend(handles=hand1, loc='upper right', handlelength=3.5,
                          scatterpoints=1, fontsize=xy_font_s - 8)
        leg2 = plt.legend(handles=hand2, loc='lower left', handlelength=3.5,
                          scatterpoints=1, fontsize=xy_font_s - 8)
        leg1.get_frame().set_alpha(0.85)
        leg2.get_frame().set_alpha(0.85)
        plt.gca().add_artist(leg1)


def make_amr_plot(in_params, amr_lit, amr_asteca):
    '''
    Make age-metallicity relation plot for both galaxies.
    '''

    rad_pc = in_params['rad_pc']

    fig = plt.figure(figsize=(5.25, 13.5))
    gs = gridspec.GridSpec(3, 1)

    amr_lit_smc, amr_lit_lmc = amr_lit
    age_vals, met_weighted, age_gyr, feh_f, age_rang_MCs, ast_lit = amr_asteca

    amr_lst = [
        [gs, 0, age_vals, met_weighted, age_gyr, age_rang_MCs, feh_f, rad_pc,
         '', '$[Fe/H]$', ast_lit],
        [gs, 1, age_vals, met_weighted, age_gyr, amr_lit_lmc, [], [],
         '', '$[Fe/H]$', ast_lit],
        [gs, 2, age_vals, met_weighted, age_gyr, amr_lit_smc, [], [],
         '$Age\,[Gyr]$', '$[Fe/H]$', ast_lit]
    ]

    for pl_params in amr_lst:
        pl_amr(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/AMR_asteca.png', dpi=300, bbox_inches='tight')


def h03_p12_mass_plots(pl_params):
    '''
    Generate ASteCA vs literature mass values plot.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, xarr, yarr, \
        carr, v_min_mp, v_max_mp, par_mean_std, m_limit = pl_params

    xy_font_s = 21
    ax = plt.subplot(gs[i], aspect='auto')
    # Text box.
    ob = offsetbox.AnchoredText(m_limit, loc=1, prop=dict(size=xy_font_s - 4))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # 0 line
    plt.axhline(y=par_mean_std[0], xmin=0, xmax=1, color='k', ls='--')
    # Shaded one sigma region.
    if par_mean_std[0] != par_mean_std[1]:
        plt.axhspan(par_mean_std[0] - par_mean_std[1],
                    par_mean_std[0] + par_mean_std[1], facecolor='grey',
                    alpha=0.3, zorder=1)

    plt.xticks(fontsize=xy_font_s - 6)
    plt.yticks(fontsize=xy_font_s - 6)
    cm = plt.cm.get_cmap('RdYlBu_r')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()

    # Plot all clusters in dictionary.
    SC = plt.scatter(xarr, yarr, marker='o', c=carr, s=80, lw=0.25,
                     cmap=cm, vmin=v_min_mp, vmax=v_max_mp, zorder=3)
    # # Text box.
    # text1 = r'$\overline{{\Delta M}}={:.1f}$'.format(par_mean_std[0])
    # r = 2 if i == 0 else 0
    # text2 = r'$\pm{:g}\,[M_{{\odot}}]$'.format(
    #     round(par_mean_std[1], r))
    # text = text1 + text2
    # ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=xy_font_s - 4))
    # ob.patch.set(alpha=0.5)
    # ax.add_artist(ob)
    # Position colorbar.
    if i == 2:
        axColor = plt.axes([0.885, 0.2, 0.1, 0.023])
        cbar = plt.colorbar(SC, cax=axColor, orientation="horizontal")
        cbar.set_label(z_lab, fontsize=xy_font_s - 3, labelpad=-45)
        cbar.set_ticks([-2., -1., 0., 1.])
        cbar.ax.tick_params(labelsize=xy_font_s - 10)


def make_cross_match_h03_p12(cross_match_h03_p12):
    """
    Plot H03 versus P12 cross matched clusters.
    """
    a_h03, a_p12, m_h03, m_p12 = cross_match_h03_p12

    a_l, a_m, a_g, m_l, m_m, m_g = [[], []], [[], []], [[], []], [[], []],\
        [[], []], [[], []]
    m_low, m_med = 5000., 20000.
    for i, a_h in enumerate(a_h03):
        a_p, m_h, m_p = a_p12[i], m_h03[i], m_p12[i]
        avr_mass = .5 * (m_h + m_p)
        # Separate by average mass limit.
        if avr_mass <= m_low:
            a_l[0].append(a_p)
            a_l[1].append(a_h)
            m_l[0].append(m_p)
            m_l[1].append(m_h)
        elif m_low < avr_mass <= m_med:
            a_m[0].append(a_p)
            a_m[1].append(a_h)
            m_m[0].append(m_p)
            m_m[1].append(m_h)
        elif m_med < avr_mass:
            a_g[0].append(a_p)
            a_g[1].append(a_h)
            m_g[0].append(m_p)
            m_g[1].append(m_h)
            if avr_mass > 100000.:
                print 'Large average DB mass OCs:', m_h, m_p, avr_mass

    # Mean & StandDev. P12-H03
    # Low mass region.
    m_l_delta = (np.array(m_l[0]) - np.array(m_l[1])) /\
        (np.array(m_l[0]) + np.array(m_l[1]))
    m_l_mean_std = [np.mean(m_l_delta), np.std(m_l_delta)]
    print 'Low mass Delta M_r mean +- std:', m_l_mean_std
    a_l_delta = np.array(a_l[0]) - np.array(a_l[1])
    # Medium mass region.
    m_m_delta = (np.array(m_m[0]) - np.array(m_m[1])) /\
        (np.array(m_m[0]) + np.array(m_m[1]))
    m_m_mean_std = [np.mean(m_m_delta), np.std(m_m_delta)]
    print 'Medium mass Delta M_r mean +- std:', m_m_mean_std
    a_m_delta = np.array(a_m[0]) - np.array(a_m[1])
    # Large mass region.
    m_g_delta = (np.array(m_g[0]) - np.array(m_g[1])) /\
        (np.array(m_g[0]) + np.array(m_g[1]))
    m_g_mean_std = [np.mean(m_g_delta), np.std(m_g_delta)]
    print 'Large mass Delta M_r mean +- std:', m_g_mean_std
    a_g_delta = np.array(a_g[0]) - np.array(a_g[1])

    scale = 10**4
    mass_l_avrg = ((np.array(m_l[0]) + np.array(m_l[1])) * 0.5)
    mass_m_avrg = ((np.array(m_m[0]) + np.array(m_m[1])) * 0.5) / scale
    mass_g_avrg = ((np.array(m_g[0]) + np.array(m_g[1])) * 0.5) / scale

    # Generate ASteca vs literature plots.
    fig = plt.figure(figsize=(19.3, 6.3))
    gs = gridspec.GridSpec(1, 3)

    xmin, xmax = [-50., 0.5, 2.05], [5000., 2.05, 10.5]
    ymin, ymax = [-1.19, -1.19, -1.29], [1.3, 1.19, 0.7]

    cbar_min = min(a_l_delta.min(), a_m_delta.min(), a_g_delta.min())
    cbar_max = max(a_l_delta.max(), a_m_delta.max(), a_g_delta.max())
    print 'Delta age min, max:', cbar_min, cbar_max

    as_lit_pl_lst = [
        # Low H03 vs P12 masses.
        [gs, 0, xmin[0], xmax[0], ymin[0], ymax[0],
         r'$\overline{M}_{DBs}\,[M_{\odot}]$',
         r'$\overline{\Delta M_r}\;\;(P12-H03)$', '',
         mass_l_avrg, m_l_delta, a_l_delta, cbar_min, cbar_max,
         m_l_mean_std, r'$\overline{M}_{DBs}\leq 5000\,[M_{\odot}]$'],
        # Medium H03 vs P12 masses.
        [gs, 1, xmin[1], xmax[1], ymin[1], ymax[1],
         r'$\overline{M}_{DBs}\,[10^{-4} M_{\odot}]$', '', '',
         mass_m_avrg, m_m_delta, a_m_delta, cbar_min, cbar_max,
         m_m_mean_std, r'$5000<\overline{M}_{DBs}\leq 20000\,[M_{\odot}]$'],
        # Large H03 vs P12 masses.
        [gs, 2, xmin[2], xmax[2], ymin[2], ymax[2],
         r'$\overline{M}_{DBs}\,[10^{-4} M_{\odot}]$', '',
         r'$\Delta \log(age/yr)$',
         mass_g_avrg, m_g_delta, a_g_delta, cbar_min, cbar_max,
         m_g_mean_std, r'$\overline{M}_{DBs}>20000\,[M_{\odot}]$']
    ]

    for pl_params in as_lit_pl_lst:
        h03_p12_mass_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/H03_P12_mass.png', dpi=300, bbox_inches='tight')


def age_mass_corr_plot(pl_params):
    '''
    Generate plots for the cross-matched isochrone fitted OCs.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, data, labels, mark, cols, \
        kde_cont = pl_params

    xy_font_s = 21
    ax = plt.subplot(gs[i])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    ax.tick_params(labelsize=13)
    if i == 1:
        ax.set_yticklabels([])

    # Origin lines.
    plt.plot([-10, 10], [0., 0.], 'k', ls='--')
    plt.plot([0., 0.], [-10, 10], 'k', ls='--')
    # Plot all clusters.
    xarr, yarr = data[0], data[1]
    plt.scatter(xarr, yarr, marker=mark, c=cols, s=50,
                lw=0.25, edgecolor='w', label=labels, zorder=0)
    # Plot KDE.
    x, y, kde = kde_cont
    plt.contour(x, y, kde, 5, colors='k', linewidths=1.5, zorder=5)

    # Legend.
    leg = plt.legend(loc='upper left', markerscale=1., scatterpoints=1,
                     fontsize=xy_font_s - 7)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.5)
    ax.set_aspect('auto')


def make_age_mass_corr(cross_match, cross_match_h03_p12):
    '''
    Plot the differences between extinction and age for ASteCA values versus
    Washington values (ie: Piatti et al. values) and ASteCA values versus
    the databases where the isochrone fitting method was used.
    '''
    # unpack databases.
    h03, p12 = cross_match[2], cross_match[6]
    a_h03, a_p12, m_h03, m_p12 = cross_match_h03_p12

    # Define lists of difference between ages and extinctions.
    # Age indexes (DB, ASteCA) -> 2, 4
    # Mass indexes (DB, ASteCA) -> 8, 10

    # ASteCA minus database log(age) diffs.
    diffs_ages_h03 = np.array(h03[4]) - np.array(h03[2])
    diffs_ages_p12 = np.array(p12[4]) - np.array(p12[2])
    # Combine into single list
    as_h03_p12_ages = list(diffs_ages_h03) + list(diffs_ages_p12)

    # ASteCA minus database log(mass) diffs.
    diffs_mass_h03 = np.log(np.array(h03[10])) - np.log(np.array(h03[8]))
    diffs_mass_p12 = np.log(np.array(p12[10])) - np.log(np.array(p12[8]))
    # Combine into single list
    as_h03_p12_mass = list(diffs_mass_h03) + list(diffs_mass_p12)

    # P12-H03 ages and masses.
    p12_h03_ages = np.array(a_p12) - np.array(a_h03)
    p12_h03_mass = np.log(np.array(m_p12)) - np.log(np.array(m_h03))

    # Order data to plot.
    age_mass_data = [[as_h03_p12_ages, as_h03_p12_mass],
                     [p12_h03_ages, p12_h03_mass]]

    xmm = [-1.95, 1.95]
    ymm = [-3.95, 3.95]
    # Obtain a Gaussian KDE for each plot.
    # Define x,y grid.
    gd_c = complex(0, 100)
    kde_cont = []
    for xarr, yarr in age_mass_data:
        values = np.vstack([xarr, yarr])
        kernel = stats.gaussian_kde(values)
        xmin, xmax, ymin, ymax = xmm[0], xmm[1], ymm[0], ymm[1]
        x, y = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
        positions = np.vstack([x.ravel(), y.ravel()])
        # Evaluate kernel in grid positions.
        k_pos = kernel(positions)
        kde = np.reshape(k_pos.T, x.shape)
        kde_cont.append([x, y, kde])

    # Define names of arrays being plotted.
    x_lab = ['$\Delta \log(age/yr)_{\mathtt{ASteCA}-DBs}$',
             '$\Delta \log(age/yr)_{P12-H03}$']
    y_lab = ['$\Delta \log(M/M_{\odot})$', '']
    labels = ['$\mathtt{ASteCA}-DBs$', 'P12-H03']
    mark = ['>', '^']
    cols = ['c', 'm']

    # Arbitrary size so plots are actually squared.
    fig = plt.figure(figsize=(17, 6))
    gs = gridspec.GridSpec(1, 3)

    cross_match_lst = [
        # Age vs mass delta plots for ASteCA-DBs
        [gs, 0, xmm[0], xmm[1], ymm[0], ymm[1], x_lab[0], y_lab[0],
            age_mass_data[0], labels[0], mark[0], cols[0], kde_cont[0]],
        # Age vs mass delta plots for P12-H03
        [gs, 1, xmm[0], xmm[1], ymm[0], ymm[1], x_lab[1], y_lab[1],
            age_mass_data[1], labels[1], mark[1], cols[1], kde_cont[1]]
    ]

    for pl_params in cross_match_lst:
        age_mass_corr_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/age_mass_corr.png', dpi=300, bbox_inches='tight')
