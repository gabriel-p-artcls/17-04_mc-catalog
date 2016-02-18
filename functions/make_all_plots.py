
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Ellipse
# import statsmodels.api as sm
from scipy import stats
from scipy.stats import ks_2samp

from ra_dec_map import ra_dec_plots
from kde_2d import kde_map
from amr_kde import age_met_rel
import lin_fit_conf_bands as lf_cb


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
    if x_lab == '$[Fe/H]_{ASteCA}$':
        # 1% of axis ranges.
        ax_ext = (xmax - xmin) * 0.01
    elif x_lab == '$(m-M)_{0;\,ASteCA}$':
        # 5% of axis ranges.
        ax_ext = (xmax - xmin) * 0.05
    else:
        # No scatter.
        ax_ext = 0.
    # Add randoms scatter.
    rs_x = xarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))
    rs_y = yarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))

    # Plot all clusters in dictionary.
    SC = plt.scatter(rs_x, rs_y, marker='o', c=zarr, s=110, lw=0.25, cmap=cm,
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
    if i in [0, 1]:
        # Text box.
        ob = offsetbox.AnchoredText(gal_name, loc=4,
                                    prop=dict(size=xy_font_s - 4))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    if i in [2, 5, 8, 11, 14]:
        # Text box.
        pres = [2, 2] if i in [2, 5, 8, 11] else [0, 0]
        text1 = r'$\bar{{y}}={:g}$'.format(round(par_mean_std[0], pres[0]))
        text2 = r'$\sigma={:g}$'.format(round(par_mean_std[1], pres[1]))
        text = text1 + '\n' + text2
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

    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, marr, msigma, \
        rarr = [in_params[_] for _ in
                ['zarr', 'zsigma', 'aarr', 'asigma', 'earr', 'esigma', 'darr',
                'dsigma', 'marr', 'msigma', 'rarr']]

    # SMC/LMC
    z_all, age_all, ext_all, dm_all, ma_all = [], [], [], [], []
    z_delta, age_delta, ext_delta, dm_delta, ma_delta = [], [], [], [], []
    for k in [0, 1]:

        # \delta z as ASteCA - literature values.
        z_all += zarr[k][0]
        z_delta += list(np.array(zarr[k][0]) - np.array(zarr[k][1]))
        # \delta log(age) as ASteCA - literature values.
        age_all += aarr[k][0]
        age_delta += list(np.array(aarr[k][0]) - np.array(aarr[k][1]))
        # \delta E(B-V) as ASteCA - literature values.
        ext_all += earr[k][0]
        ext_delta += list(np.array(earr[k][0]) - np.array(earr[k][1]))
        # \delta dm as ASteCA - literature values.
        dm_all += darr[k][0]
        dm_delta += list(np.array(darr[k][0]) - np.array(darr[k][1]))
        # \delta mass as ASteCA - literature values.
        ma_all += marr[k][0]
        ma_delta += list(np.array(marr[k][0]) - np.array(marr[k][1]))

    # Shaded area that contains 9X% of the clusters.
    # par_9x_span = []
    # idx_9x = int(68 * len(zarr[k][0]) / 100)
    # for span in [z_delta, age_delta, ext_delta, dm_delta]:
    #     abs_v = sorted([abs(_) for _ in span])
    #     par_9x_span.append(abs_v[idx_9x])

    # K_S test.
    # Null hypothesis: that 2 independent samples are drawn from the same
    # continuous distribution (sample sizes can be different)
    #
    # If the K-S statistic is small or the p-value is high, then we cannot
    # reject the hypothesis that the distributions of the two samples are the
    # same.
    # For two identical distributions the KS value will be small and
    # the p-value high.

    print 'SMC'
    # Mean only for those clusters with ASteCA age values closer than 0.5
    # to literature values.
    age_smc_f = []
    for a in list(np.array(aarr[0][0]) - np.array(aarr[0][1])):
        if abs(a) <= 0.5:
            age_smc_f.append(a)
    print 'Age mean for Delta log(age)<0.5:', np.mean(age_smc_f)
    print 'Met vals mean/std, AS:', np.mean(zarr[0][0]), \
        np.std(zarr[0][0])
    print 'Met vals mean/std, Lit:', np.mean(zarr[0][1]), \
        np.std(zarr[0][1])
    print 'Met vals CCC:', ccc(zarr[0][0], zarr[0][1])
    print 'Met vals PCC', np.corrcoef(zarr[0][0], zarr[0][1])[0, 1]
    ks, pval = ks_2samp(zarr[0][0], zarr[0][1])
    print 'Met vals K-S:', ks, pval
    print 'Age vals CCC:', ccc(aarr[0][0], aarr[0][1])
    print 'Age vals PCC:', np.corrcoef(aarr[0][0], aarr[0][1])[0, 1]
    ks, pval = ks_2samp(aarr[0][0], aarr[0][1])
    print 'Age vals K-S:', ks, pval, '\n'

    print 'LMC'
    # Mean only for those clusters with ASteCA age values closer than 0.5
    # to literature values.
    age_lmc_f = []
    for a in list(np.array(aarr[1][0]) - np.array(aarr[1][1])):
        if abs(a) <= 0.5:
            age_lmc_f.append(a)
    print 'Age mean for Delta log(age)<0.5:', np.mean(age_lmc_f)
    print 'Met vals mean/std, AS:', np.mean(zarr[1][0]), \
        np.std(zarr[1][0])
    # Filter out clusters with no metal values in the literature (.ods file)
    z_lmc_lit_f, z_lmc_ast_f = [], []
    for z_ast, z_lit in zip(*[zarr[1][0], zarr[1][1]]):
        if abs(z_lit) < 10000:
            z_lmc_lit_f.append(z_lit)
            z_lmc_ast_f.append(z_ast)
    print 'Met vals mean/std, Lit:', np.mean(z_lmc_lit_f), \
        np.std(z_lmc_lit_f)
    print 'Met vals CCC:', ccc(z_lmc_ast_f, z_lmc_lit_f)
    print 'Met vals PCC:', np.corrcoef(z_lmc_ast_f, z_lmc_lit_f)[0, 1]
    ks, pval = ks_2samp(z_lmc_ast_f, z_lmc_lit_f)
    print 'Met vals K-S:', ks, pval
    print 'Age vals CCC:', ccc(aarr[1][0], aarr[1][1])
    print 'Age vals PCC:', np.corrcoef(aarr[1][0], aarr[1][1])[0, 1]
    ks, pval = ks_2samp(aarr[1][0], aarr[1][1])
    print 'Age vals K-S:', ks, pval

    # print 'Gal  Mean  StandDev'
    par_mean_std = []
    for span in [z_delta, age_delta, ext_delta, dm_delta, ma_delta]:
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

    ext_min, ext_max = 0., 0.3
    # dm_min, dm_max = 18.62, 19.21
    # dm_min, dm_max = 18.21, 18.79
    dm_min, dm_max = 18.21, 19.21
    dm_span = (dm_max - dm_min) / 2.

    # For old runs where the dist mod range was large.
    # dm_min, dm_max = 17.8, 20.2

    as_lit_pl_lst = [
        # Metallicity LMC/SMC
        [gs, 0, -2.4, 0.45, -2.4, 0.45, '$[Fe/H]_{ASteCA}$', '$[Fe/H]_{lit}$',
            '', zarr[1][0], zsigma[1][0], zarr[1][1],
            zsigma[1][1], aarr[1][0], 6.6, 9.8, [], 'LMC'],
        [gs, 1, -2.4, 0.45, -2.4, 0.45, '$[Fe/H]_{ASteCA}$', '',
            '', zarr[0][0], zsigma[0][0], zarr[0][1],
            zsigma[0][1], aarr[0][0], 6.6, 9.8, [], 'SMC'],
        # Asteca z vs \delta z with lit values.
        [gs, 2, -2.4, 0.45, -1.43, 1.43, '$[Fe/H]_{ASteCA}$',
            '$\Delta [Fe/H]$', '$log(age/yr)_{ASteCA}$', z_all,
            [], z_delta, [], age_all, 6.6, 9.8, par_mean_std[0],
            ''],

        # Age LMC/SMC
        [gs, 3, 5.8, 10.6, 5.8, 10.6, '$log(age/yr)_{ASteCA}$',
            '$log(age/yr)_{lit}$', '', aarr[1][0],
            asigma[1][0], aarr[1][1], asigma[1][1], earr[1][0], ext_min,
            ext_max, [], ''],
        [gs, 4, 5.8, 10.6, 5.8, 10.6, '$log(age/yr)_{ASteCA}$',
            '', '', aarr[0][0], asigma[0][0], aarr[0][1], asigma[0][1],
            earr[0][0], ext_min, ext_max, [], ''],
        # Asteca log(age) vs \delta log(age) with lit values.
        [gs, 5, 5.8, 10.6, -2.4, 2.4, '$log(age/yr)_{ASteCA}$',
            '$\Delta log(age/yr)$', '$E(B-V)_{ASteCA}$', age_all,
            [], age_delta, [], ext_all, ext_min, ext_max, par_mean_std[1], ''],

        # Ext LMC/SMC
        [gs, 6, -0.04, 0.29, -0.04, 0.29, '$E(B-V)_{ASteCA}$',
            '$E(B-V)_{lit}$', '', earr[1][0], esigma[1][0], earr[1][1],
            esigma[1][1], aarr[1][0], 6.6, 9.8, [], ''],
        [gs, 7, -0.04, 0.29, -0.04, 0.29, '$E(B-V)_{ASteCA}$',
            '', '', earr[0][0], esigma[0][0], earr[0][1], esigma[0][1],
            aarr[0][0], 6.6, 9.8, [], ''],
        # Asteca E(B-V) vs \delta E(B-V) with lit values.
        [gs, 8, -0.04, 0.29, -0.21, 0.21, '$E(B-V)_{ASteCA}$',
            '$\Delta E(B-V)$', '$log(age/yr)_{ASteCA}$', ext_all, [],
            ext_delta, [], age_all, 6.6, 9.8, par_mean_std[2], ''],

        # Dits mod LMC/SMC
        [gs, 9, dm_min, dm_max, dm_min, dm_max, '$(m-M)_{0;\,ASteCA}$',
            '$(m-M)_{0;\,lit}$', '', darr[1][0], dsigma[1][0], darr[1][1],
            dsigma[1][1], aarr[1][0], 6.6, 9.8, [], ''],
        [gs, 10, dm_min, dm_max, dm_min, dm_max, '$(m-M)_{0;\,ASteCA}$',
            '', '', darr[0][0], dsigma[0][0], darr[0][1], dsigma[0][1],
            aarr[0][0], 6.6, 9.8, [], ''],
        # Asteca dist_mod vs \delta dist_mod with lit values.
        [gs, 11, dm_min, dm_max, -1. * dm_span, dm_span,
            '$(m-M)_{0;\,ASteCA}$', '$\Delta (m-M)_{0}$',
            '$log(age/yr)_{ASteCA}$', dm_all, [],
            dm_delta, [], age_all, 6.6, 9.8, par_mean_std[3], ''],

        # # ASteCA vs literature masses.
        # [gs, 8, 10., 5000., 10., 5000., '$M_{ASteCA}\,(M_{\odot})$',
        #     '$M_{lit}\,(M_{\odot})$', '$log(age/yr)_{ASteCA}$', marr[k][0],
        #     msigma[k][0], marr[k][1], msigma[k][1], aarr[k][0], 6.6, 9.8, [],
        #     galax],
        # [gs, 9, 10., 4000., -2000., 2000., '$M_{ASteCA}\,(M_{\odot})$',
        #     '$\Delta M_{\odot}$', '$log(age/yr)_{ASteCA}$', marr[k][0],
        #     msigma[k][0], ma_delta, [], aarr[k][0], 6.6, 9.8, par_mean_std[4],
        #     galax]
    ]
    #
    for pl_params in as_lit_pl_lst:
        as_vs_lit_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_vs_lit_S-LMC.png', dpi=300, bbox_inches='tight')


def kde_plots(pl_params):
    '''
    Generate KDE plots.
    '''
    gs, i, x_lab, y_lab, xarr, xsigma, yarr, ysigma, x_rang, y_rang = pl_params

    ext = [x_rang[0], x_rang[1], y_rang[0], y_rang[1]]

    # Generate maps.
    z = kde_map(np.array(xarr), np.array(xsigma), np.array(yarr),
                np.array(ysigma), ext, 100)

    # Make plot.
    ax = plt.subplot(gs[i])
    xy_font_s = 18
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)

    # cm = plt.cm.get_cmap('RdYlBu_r')
    cm = plt.cm.gist_earth_r
    ax.imshow(z, cmap=cm, extent=ext)
    ax.set_aspect('auto')
    # Error bars.
    # plt.errorbar(xarr, yarr, xerr=xsigma, yerr=ysigma, fmt='none',
    #              elinewidth=0.4, color='k')
    # Define 1% of axis ranges.
    xax_ext = (ext[1] - ext[0]) * 0.001
    yax_ext = (ext[3] - ext[2]) * 0.001
    # Random scatter.
    rs_x = np.random.uniform(0., xax_ext, len(xarr))
    rs_y = np.random.uniform(0., yax_ext, len(xarr))
    # Clusters.
    # color='#6b6868'
    plt.scatter(xarr + rs_x, yarr + rs_y, marker='*', color='r', s=40,
                lw=0.5, facecolors='none')
    ax.set_xlim(ext[0], ext[1])
    ax.set_ylim(ext[2], ext[3])


def make_kde_plots(galax, k, in_params):
    '''
    Prepare parameters and call function to generate SMC and LMC KDE plots.
    '''
    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, marr, msigma = \
        [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
                                'esigma', 'darr', 'dsigma', 'marr', 'msigma']]

    fig = plt.figure(figsize=(14, 25))  # create the top-level container
    gs = gridspec.GridSpec(4, 2)  # create a GridSpec object

    # Define extension for each parameter range.
    age_rang, fe_h_rang, mass_rang = [6.4, 10.1], [-2.4, 0.15], [-100., 30500.]
    if galax == 'SMC':
        E_bv_rang, dist_mod_rang = [-0.01, 0.15], [18.75, 19.25]
    else:
        E_bv_rang, dist_mod_rang = [-0.01, 0.3], [18.25, 18.75]

    # Age in Gyrs.
    age_gyr = [10 ** (np.asarray(aarr[k][0]) - 9),
               np.asarray(asigma[k][0]) * np.asarray(aarr[k][0]) *
               np.log(10) / 5.]
    age_gyr_rang = [0., 6.6]

    kde_pl_lst = [
        [gs, 0, '$log(age/yr)_{ASteCA}$', '$[Fe/H]_{ASteCA}$', aarr[k][0],
            asigma[k][0], zarr[k][0], zsigma[k][0], age_rang, fe_h_rang],
        [gs, 1, '$log(age/yr)_{ASteCA}$', '$M_{ASteCA}\,(M_{\odot})$',
            aarr[k][0], asigma[k][0], marr[k][0], msigma[k][0], age_rang,
            mass_rang],
        [gs, 2, '$(m-M)_{\circ;\,ASteCA}$', '$E(B-V)_{ASteCA}$', darr[k][0],
            dsigma[k][0], earr[k][0], esigma[k][0], dist_mod_rang, E_bv_rang],
        [gs, 3, '$M_{ASteCA}\,(M_{\odot})$', '$[Fe/H]_{ASteCA}$', marr[k][0],
            msigma[k][0], zarr[k][0], zsigma[k][0], mass_rang, fe_h_rang],
        [gs, 4, '$Age_{ASteCA}\,(Gyr)$', '$[Fe/H]_{ASteCA}$', age_gyr[0],
            age_gyr[1], zarr[k][0], zsigma[k][0], age_gyr_rang, fe_h_rang]
    ]
    #
    for pl_params in kde_pl_lst:
        kde_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_kde_maps_' + galax + '.png', dpi=300,
                bbox_inches='tight')


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
         '$M\,(M_{\odot})$']
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
    x_lab = '$E(B-V)_{ASteCA}$'
    y_lab = ['$E(B-V)_{MCEV,\,closer}$', '$E(B-V)_{MCEV,\,max}$']
    z_lab = ['$log(age/yr)_{ASteCA}$', '$E(B-V)_{SF}$', '$dist\,(deg)$']

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
    x_lab, y_lab, z_lab = '$log(age/yr)_{ASteCA}$', \
        '$(C-T_{1})_{0;\,ASteCA}$', '$M\,(M_{\odot})$'

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
    x_lab, y_lab, z_lab = ['$log(age/yr)_{ASteCA}$', '$[Fe/H]_{ASteCA}$'], \
        '$Concentration\,(N_{memb}/pc^{2})$', '$M\,(M_{\odot})$'

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
    x_lab = ['$R_{cl;\,ASteCA}\,(pc)$', '$R_{core;\,ASteCA}\,(pc)$']
    y_lab = ['$log(age/yr)_{ASteCA}$', '$[Fe/H]_{ASteCA}$', '$M\,(M_{\odot})$']
    z_lab = ['$M\,(M_{\odot})$', '$log(age/yr)_{ASteCA}$']

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

    print '\n* Clusters with n_memb > 50 & prob < 0.5'
    for k, gal in enumerate(['SMC', 'LMC']):
        for i, n_m in enumerate(n_memb[k]):
            if n_m > 50 and kde_prob[k][i] < 0.5:
                print '', gal, gal_names[k][i], n_m, kde_prob[k][i]

    # Define names of arrays being plotted.
    x_lab, y_lab, z_lab = '$CI_{ASteCA}$', '$prob_{ASteCA}$', \
        ['$log(age/yr)_{ASteCA}$', '$[Fe/H]_{ASteCA}$', '$M\,(M_{\odot})$',
            '$M\,(M_{\odot})$', '$log(age/yr)_{ASteCA}$']
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
        [gs, 4, -9, 1000, ymin, ymax, '$N_{ASteCA}$', y_lab, z_lab[0],
         n_memb[0], kde_prob[0], aarr[0][0], rad_pc[0], 'SMC'],
        [gs, 5, -9, 1000, ymin, ymax, '$N_{ASteCA}$', y_lab, z_lab[0],
         n_memb[1], kde_prob[1], aarr[1][0], rad_pc[1], 'LMC']
    ]

    for pl_params in prob_CI_pl_lst:
        prob_vs_CI_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_prob_vs_CI.png', dpi=300, bbox_inches='tight')


def plot_dist_2_cent(pl_params):
    '''
    Generate plots for KDE probabilities versus contamination indexes.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, xarr, yarr,\
        zarr, ysigma, v_min, v_max, rad, gal_name = pl_params
    siz = np.asarray(rad) * 6

    xy_font_s = 16
    cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i])
    # ax.set_aspect('auto')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    if i in [0, 2, 4, 6]:
        plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    # Plot all clusters in dictionary.
    # Introduce random scatter.
    # X% of axis ranges.
    ax_ext = (xmax - xmin) * 0.02
    # Add randoms scatter.
    if i in [2, 3]:
        rs_x = xarr + np.random.uniform(-ax_ext, ax_ext, len(xarr))
    else:
        rs_x = xarr
    SC = plt.scatter(rs_x, yarr, marker='o', c=zarr, s=siz, lw=0.25, cmap=cm,
                     vmin=v_min, vmax=v_max, zorder=3)
    # Plot y error bar if it is passed.
    if ysigma:
        plt.errorbar(xarr, yarr, yerr=ysigma, ls='none', color='k',
                     elinewidth=0.4, zorder=1)

    if i in [2, 3]:
        # Linear regression of metallicity, NOT weighted by errors.
        fit_nw = lf_cb.non_weigth_linear_fit(xarr, yarr)
        plt.plot(xarr, fit_nw(xarr), '-k', lw=0.8)
        # Linear regression and confidence bands or metallicity gradient,
        # weighted by their errors in [Fe/H].
        a, b, sa, sb, rchi2, dof = lf_cb.weight_linear_fit(
            np.array(xarr), np.array(yarr), np.array(ysigma))
        # Confidence bands.
        lcb, ucb, x = lf_cb.confband(np.array(xarr), np.array(yarr), a, b,
                                     conf=0.95)
        fit_w = np.poly1d([a, b])
        plt.plot(x, fit_w(x), '--g')
        plt.fill_between(x, lcb, ucb, alpha=0.3, facecolor='gray')
        # Text box.
        text = r'$[Fe/H]_{{rg}}={:.2f}\pm{:.2f}\,dex\,kpc^{{-1}}$'.format(
            a * 1000., sa * 1000.)
        ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=xy_font_s - 3))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)

    if gal_name != '':
        # Text box.
        ob = offsetbox.AnchoredText(gal_name, loc=1,
                                    prop=dict(size=xy_font_s - 2))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
    # if i in [1, 3, 5, 7]:
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    zpad = 10 if z_lab == '$E_{(B-V)}$' else 5
    cbar.set_label(z_lab, fontsize=xy_font_s, labelpad=zpad)


def make_dist_2_cents(in_params):
    '''
    Plot ASteCA distances to center of either MC galaxy.
    '''

    zarr, zsigma, aarr, asigma, earr, esigma, marr, msigma, rad_pc, cont_ind,\
        dist_cent, gal_names, ra, dec = \
        [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
                                'esigma', 'marr', 'msigma', 'rad_pc',
                                'cont_ind', 'dist_cent', 'gal_names', 'ra',
                                'dec']]

    # # Print info to screen.
    # for j, gal in enumerate(['SMC', 'LMC']):
    #     for i, cl in enumerate(gal_names[j]):
    #         if dist_cent[j][i] > 4000 and aarr[j][0][i] < 8.5:
    #             print gal, cl, ra[j][i], dec[j][i], dist_cent[j][i],\
    #                 '{:.5f}'.format(zarr[j][0][i]), aarr[j][0][i]
    #
    # print 'SMC, ASteCA:', np.mean(zarr[0][0]), np.std(zarr[0][0])
    # print 'LMC, ASteCA:', np.mean(zarr[1][0]), np.std(zarr[1][0])

    # Define names of arrays being plotted.
    x_lab, yz_lab = '$R_{GC}\,[pc]$', \
        ['$log(age/yr)_{ASteCA}$', '$[Fe/H]_{ASteCA}$', '$M\,(M_{\odot})$',
            '$E(B-V)_{ASteCA}$']
    xmin, xmax = 0, 7500
    vmin_met, vmax_met = -2.1, 0.29
    # vmin_mas, vmax_mas = 1000, 28000
    vmin_ext, vmax_ext = 0., 0.3
    vmin_age, vmax_age = 6.7, 9.7

    fig = plt.figure(figsize=(14, 25))
    gs = gridspec.GridSpec(4, 2)

    dist_2_cent_pl_lst = [
        # SMC
        [gs, 0, xmin, xmax, 6.6, 10.1, x_lab, yz_lab[0], yz_lab[1],
            dist_cent[0], aarr[0][0], zarr[0][0], asigma[0][0], vmin_met,
            vmax_met, rad_pc[0], 'SMC'],
        [gs, 2, xmin, xmax, -2.4, 0.4, x_lab, yz_lab[1], yz_lab[0],
            dist_cent[0], zarr[0][0], aarr[0][0], zsigma[0][0], vmin_age,
            vmax_age, rad_pc[0], 'SMC'],
        [gs, 4, xmin, xmax, 0., 30000, x_lab, yz_lab[2], yz_lab[3],
            dist_cent[0], marr[0][0], earr[0][0], msigma[0][0], vmin_ext,
            vmax_ext, rad_pc[0], ''],
        [gs, 6, xmin, xmax, -0.01, 0.11, x_lab, yz_lab[3], yz_lab[0],
            dist_cent[0], earr[0][0], aarr[0][0], esigma[0][0], vmin_age,
            vmax_age, rad_pc[0], ''],
        # LMC
        [gs, 1, xmin, xmax, 6.6, 10.1, x_lab, yz_lab[0], yz_lab[1],
            dist_cent[1], aarr[1][0], zarr[1][0], asigma[1][0], vmin_met,
            vmax_met, rad_pc[1], 'LMC'],
        [gs, 3, xmin, xmax, -2.4, 0.4, x_lab, yz_lab[1], yz_lab[0],
            dist_cent[1], zarr[1][0], aarr[1][0], zsigma[1][0], vmin_age,
            vmax_age, rad_pc[1], 'LMC'],
        [gs, 5, xmin, xmax, 0., 30000, x_lab, yz_lab[2], yz_lab[3],
            dist_cent[1], marr[1][0], earr[1][0], msigma[1][0], vmin_ext,
            vmax_ext, rad_pc[1], ''],
        [gs, 7, xmin, xmax, -0.01, 0.31, x_lab, yz_lab[3], yz_lab[0],
            dist_cent[1], earr[1][0], aarr[1][0], esigma[1][0], vmin_age,
            vmax_age, rad_pc[1], '']
    ]

    for pl_params in dist_2_cent_pl_lst:
        plot_dist_2_cent(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_dist_2_cent.png', dpi=300, bbox_inches='tight')


def cross_match_plot(pl_params):
    '''
    Generate plots for the cross-matched age and mass values.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, indexes, labels, \
        mark, cols, text_box, databases = pl_params

    a, e_a, b, e_b = indexes

    xy_font_s = 21
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i])
    ax.set_aspect('equal')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    # Plot all clusters for each DB.
    for j, DB in enumerate(databases):
        if DB:
            xarr, yarr = DB[a], DB[b]
            xsigma, ysigma = DB[e_a], DB[e_b]

            # # Fit y = s*x + i line.
            # from scipy import stats
            # slope, intrcpt, r_v, p_v, std_err = stats.linregress(xarr, yarr)
            # print '\n', labels[j]
            # print 'y=a*x+b fit:', slope, intrcpt, r_v ** 2, std_err, '\n'

            # # Fit y = s*x line to data, ie: x=0 --> y=0 (intercept=0).
            # # x needs to be a column vector instead of a 1D vector for this.
            # x = np.asarray(xarr)[:, np.newaxis]
            # lin_fit = np.linalg.lstsq(x, yarr)
            # print lin_fit
            # slope = lin_fit[0][0]

            # # Fit y = s*x + i line.
            # model = sm.OLS(yarr, xarr)
            # results = model.fit()
            # # print results.summary()
            # slope, std_err = results.params[0], results.bse[0]
            # db_lab = labels[j] + '$\;(N={},\,s={:.2f}),\,SE={:.3f}$'.format(
            #     len(xarr), slope, std_err)

            db_lab = labels[j] + '$\;(N={})$'.format(len(xarr))
            # Star marker is too small compared to the rest.
            siz = 60. if mark[j] != '*' else 90.
            plt.scatter(xarr, yarr, marker=mark[j], c=cols[j], s=siz,
                        lw=0.25, edgecolor='w', label=db_lab, zorder=3)
            # Plot error bars.
            if xsigma:
                for k, xy in enumerate(zip(*[xarr, yarr])):
                    x_err = xsigma[k] if 0. < xsigma[k] < 5. else 0.
                    plt.errorbar(xy[0], xy[1], xerr=x_err, ls='none',
                                 color='k', elinewidth=0.2, zorder=1)
            if ysigma:
                for k, xy in enumerate(zip(*[xarr, yarr])):
                    y_err = ysigma[k] if 0. < ysigma[k] < 5. else 0.
                    plt.errorbar(xy[0], xy[1], yerr=y_err, ls='none',
                                 color='k', elinewidth=0.2, zorder=1)
            # Legend.
            if y_lab != '$mass_{DB}\,[M_{\odot}]$':
                leg = plt.legend(loc='upper left', markerscale=1.,
                                 scatterpoints=1, fontsize=xy_font_s - 7)
                leg.get_frame().set_alpha(0.5)
    plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')  # 1:1 line
    if text_box:
        # Text box.
        ob = offsetbox.AnchoredText(text_box, loc=4,
                                    prop=dict(size=xy_font_s - 5))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)


def make_cross_match(cross_match):
    '''
    Plot ASteCA ages and masses versus the values found in several databases.
    '''
    # unpack databases.
    p99, p00, h03, r05, c06, g10, p12 = cross_match

    # Labels for each defined plot.
    # labels_isoch_analy = ['Pietrzynski & Udalski (1999)',
    #                       'Pietrzynski & Udalski (2000)',
    #                       'Chiosi et al. (2006)', 'Glatt et al. (2010)']
    # labels_integ_photo = ['Hunter et al. (2003)',
    #                       'Rafelski & Zaritsky (2005)',
    #                       'Popescu et al. (2012)']
    # labels_smc = ['Hunter et al. (2003)', 'Rafelski & Zaritsky (2005)',
    #               'Chiosi et al. (2006)', 'Glatt et al. (2010)']
    # labels_lmc = ['Pietrzynski & Udalski (2000)', 'Hunter et al. (2003)',
    #               'Glatt et al. (2010)', 'Popescu et al. (2012)']
    labels_isoch_analy = ['P99', 'P00', 'C06', 'G10']
    labels_integ_photo = ['H03', 'R05', 'P12']
    labels_smc = ['P99', 'H03', 'R05', 'C06', 'G10']
    labels_lmc = ['P00', 'H03', 'G10', 'P12']
    labels_mass = ['H03', 'P12']
    labels = [labels_isoch_analy, labels_integ_photo, labels_smc, labels_lmc,
              labels_mass]

    mark = [['>', '^', 'v', '<'], ['v', '*', 'o'],
            ['>', 'v', '*', 'v', '<'], ['^', 'v', '<', 'o'], ['v', 'o']]
    cols = [['chocolate', 'r', 'c', 'g'], ['m', 'k', 'b'],
            ['chocolate', 'm', 'k', 'c', 'g'], ['r', 'm', 'g', 'b'],
            ['m', 'b']]

    # Text boxes.
    text_box = ['Isochrone fitting', 'Integrated photometry', 'SMC', 'LMC',
                '$M_{\odot}<5000$']

    # Separate SMC from LMC clusters in H03 and G10 databases.
    h03_smc, h03_lmc, g10_smc, g10_lmc = [], [], [], []
    for cl in zip(*h03):
        if cl[0] == 'SMC':
            h03_smc.append(cl)
        else:
            h03_lmc.append(cl)
    h03_smc, h03_lmc = zip(*h03_smc), zip(*h03_lmc)
    for cl in zip(*g10):
        if cl[0] == 'SMC':
            g10_smc.append(cl)
        else:
            g10_lmc.append(cl)
    g10_smc, g10_lmc = zip(*g10_smc), zip(*g10_lmc)

    # Separate clusters with mass < 5000
    h03_low_mass, p12_low_mass = [], []
    for cl in zip(*h03):
        if cl[6] <= 5000.:
            h03_low_mass.append(cl)
    h03_low_mass = zip(*h03_low_mass)
    for cl in zip(*p12):
        if cl[6] <= 5000.:
            p12_low_mass.append(cl)
    p12_low_mass = zip(*p12_low_mass)

    # Define data to pass.
    databases = [[p99, p00, c06, g10], [h03, r05, p12],
                 [p99, h03_smc, r05, c06, g10_smc],
                 [p00, h03_lmc, g10_lmc, p12],
                 [h03_low_mass, p12_low_mass], [h03, p12]]

    # First set is for the ages, second for the masses.
    indexes = [[4, 5, 2, 3], [10, 11, 8, 9]]

    # Concordance correlation coef between AsteCA ages and DBs ages.
    ast_if = list(p99[4]) + list(p00[4]) + list(c06[4]) + list(g10[4])
    dbs_if = list(p99[2]) + list(p00[2]) + list(c06[2]) + list(g10[2])
    ast_ip = list(h03[4]) + list(r05[4]) + list(p12[4])
    dbs_ip = list(h03[2]) + list(r05[2]) + list(p12[2])
    print 'AsteCA vs DBs (isoch fit) age conc corr coef:', ccc(ast_if, dbs_if)
    print 'ASteCA vs DBs (isoch fit) Pearson corr coef:', np.corrcoef(
        ast_if, dbs_if)[0, 1]
    print 'AsteCA vs DBs (integ mag) age conc corr coef:', ccc(ast_ip, dbs_ip)
    print 'ASteCA vs DBs (integ mag) Pearson corr coef:', np.corrcoef(
        ast_ip, dbs_ip)[0, 1]

    # Define names of arrays being plotted.
    x_lab = ['$log(age/yr)_{ASteCA}$', '$mass_{ASteCA}\,[M_{\odot}]$']
    y_lab = ['$log(age/yr)_{DB}$', '$mass_{DB}\,[M_{\odot}]$']
    z_lab = ['$mass_{ASteCA}\,[M_{\odot}]$', '$log(age/yr)_{ASteCA}$']
    xymin, xymax = [5.8, -69.], [10.6, 5000, 30000]

    fig = plt.figure(figsize=(16, 25))
    gs = gridspec.GridSpec(4, 2)

    cross_match_lst = [
        # Age cross-match, isoch fit.
        [gs, 0, xymin[0], xymax[0], xymin[0], xymax[0], x_lab[0], y_lab[0],
            z_lab[0], indexes[0], labels[0], mark[0], cols[0], text_box[0],
            databases[0]],

        # # Extinction G10 vs P99, P00, C06. This block needs the
        # # 'matched_clusters_G10.dat' file to be used as
        # # 'matched_clusters.dat', to produce the correct plot.
        # [gs, 1, -0.01, 0.31, -0.01, 0.321, '$E(B-V)_{G10}$', '$E(B-V)_{DB}$',
        #     z_lab[0], [13, 14, 12, 14], ['P99', 'C06'], ['>', 'v'],
        #     ['chocolate', 'c'], text_box[0], [p99, c06]]

        # Age cross-match, integrated photometry.
        [gs, 1, xymin[0], xymax[0], xymin[0], xymax[0], x_lab[0], y_lab[0],
            z_lab[0], indexes[0], labels[1], mark[1], cols[1], text_box[1],
            databases[1]],
        # Age cross-match, SMC.
        [gs, 2, xymin[0], xymax[0], xymin[0], xymax[0], x_lab[0], y_lab[0],
            z_lab[0], indexes[0], labels[2], mark[2], cols[2], text_box[2],
            databases[2]],
        # Age cross-match, LMC.
        [gs, 3, xymin[0], xymax[0], xymin[0], xymax[0], x_lab[0], y_lab[0],
            z_lab[0], indexes[0], labels[3], mark[3], cols[3], text_box[3],
            databases[3]],
        # Mass cross_match (all)
        [gs, 4, xymin[1], xymax[2], xymin[1], xymax[2], x_lab[1], y_lab[1],
            z_lab[1], indexes[1], labels[4], mark[4], cols[4], [],
            databases[5]],
        # Mass cross_match (low mass)
        [gs, 5, xymin[1], xymax[1], xymin[1], xymax[1], x_lab[1], y_lab[1],
            z_lab[1], indexes[1], labels[4], mark[4], cols[4], text_box[4],
            databases[4]]
    ]

    for pl_params in cross_match_lst:
        cross_match_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/cross_match.png', dpi=300, bbox_inches='tight')


def cross_match_age_ext_plot(pl_params):
    '''
    Generate plots for the cross-matched age and mass values.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, data, labels, mark, cols, \
        kde_cont = pl_params

    xy_font_s = 21
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    ax = plt.subplot(gs[i])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    if i in [3, 5, 7]:
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

        if i == 0:
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


def make_cross_match_age_ext(cross_match, in_params):
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
    # P99 liter minus database diffs.
    diffs_lit_db_ages_p99 = np.array(p99[6]) - np.array(p99[2])
    # Same for extinctions.
    diffs_db_exts_p99 = np.array(p99[15]) - np.array(p99[14])
    diffs_lit_db_exts_p99 = np.array(p99[17]) - np.array(p99[14])

    # P00 ASteCA minus database diffs.
    diffs_db_ages_p00 = np.array(p00[4]) - np.array(p00[2])

    # C06
    diffs_db_ages_c06 = np.array(c06[4]) - np.array(c06[2])
    diffs_lit_db_ages_c06 = np.array(c06[6]) - np.array(c06[2])
    diffs_db_exts_c06 = np.array(c06[15]) - np.array(c06[14])
    diffs_lit_db_exts_c06 = np.array(c06[17]) - np.array(c06[14])

    # G10
    diffs_db_ages_g10 = np.array(g10[4]) - np.array(g10[2])
    diffs_lit_db_ages_g10 = np.array(g10[6]) - np.array(g10[2])
    diffs_db_exts_g10 = np.array(g10[15]) - np.array(g10[14])
    diffs_lit_db_exts_g10 = np.array(g10[17]) - np.array(g10[14])

    # Calculate std, means and medians for the age differences.
    txt = ['SMC', 'LMC', 'P99', 'P00', 'C06', 'G10']
    dbs = [diffs_lit_ages_smc, diffs_lit_ages_lmc, diffs_db_ages_p99,
           diffs_db_ages_p00, diffs_db_ages_c06, diffs_db_ages_g10]
    for i, db in enumerate(dbs):
        print '{}, diff ages std = {:.3f}'.format(txt[i], np.std(db))
    #     print '{}, mean = {:.3f}'.format(txt[i], np.mean(db))
    #     print '{} median: {:.3f}'.format(txt[i], np.median(db))

    # median_db = [
    #     [[np.median(diffs_db_exts_p99), np.median(diffs_db_ages_p99)],
    #      [np.median(diffs_db_exts_c06), np.median(diffs_db_ages_c06)],
    #      [np.median(diffs_db_exts_g10), np.median(diffs_db_ages_g10)]],
    #     [[np.median(diffs_lit_exts_smc), np.median(diffs_lit_ages_smc)],
    #      [np.median(diffs_lit_exts_lmc), np.median(diffs_lit_ages_lmc)]]
    # ]
    print 'DB exts median:', np.median(list(diffs_db_exts_p99) +
                                       list(diffs_db_exts_c06) +
                                       list(diffs_db_exts_g10))
    print 'DB ages median:', np.median(list(diffs_db_ages_p99) +
                                       list(diffs_db_ages_c06) +
                                       list(diffs_db_ages_g10))
    print 'AS exts median:', np.median(list(diffs_lit_exts_smc) +
                                       list(diffs_lit_exts_lmc))
    print 'AS ages median:', np.median(list(diffs_lit_ages_smc) +
                                       list(diffs_lit_ages_lmc))

    # Obtain a Gaussian KDE for each plot.
    # Define x,y grid.
    gd_c = complex(0, 100)
    kde_cont = []
    for xarr, yarr in [
            [list(diffs_db_ages_p99) + list(diffs_db_ages_c06) +
             list(diffs_db_ages_g10), list(diffs_db_exts_p99) +
             list(diffs_db_exts_c06) + list(diffs_db_exts_g10)],
            [list(diffs_lit_ages_smc) + list(diffs_lit_ages_lmc),
             list(diffs_lit_exts_smc) + list(diffs_lit_exts_lmc)],
            [list(diffs_lit_db_ages_p99) + list(diffs_lit_db_ages_c06) +
             list(diffs_lit_db_ages_g10), list(diffs_lit_db_exts_p99) +
             list(diffs_lit_db_exts_c06) + list(diffs_lit_db_exts_g10)]
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
               [diffs_db_ages_c06, [], diffs_db_exts_c06, []],
               [diffs_db_ages_g10, [], diffs_db_exts_g10, []]]
    lit_db_data = [[diffs_lit_db_ages_p99, [], diffs_lit_db_exts_p99, []],
                   [diffs_lit_db_ages_c06, [], diffs_lit_db_exts_c06, []],
                   [diffs_lit_db_ages_g10, [], diffs_lit_db_exts_g10, []]]
    # 1:1 plots.
    ext_lit_data = [[earr[0][0], esigma[0][0], earr[0][1], esigma[0][1]],
                    [earr[1][0], esigma[1][0], earr[1][1], esigma[1][1]]]
    ext_DB_data = [[p99[15], p99[16], p99[14], []],
                   [c06[15], c06[16], c06[14], []],
                   [g10[15], g10[16], g10[14], []]]
    ext_lit_DB_data = [[p99[17], p99[18], p99[14], []],
                       [c06[17], c06[18], c06[14], []],
                       [g10[17], g10[18], g10[14], []]]
    age_lit_DB_data = [[p99[6], p99[7], p99[2], p99[3]],
                       [p00[6], p00[7], p00[2], p00[3]],
                       [c06[6], c06[7], c06[2], c06[3]],
                       [g10[6], g10[7], g10[2], g10[3]]]

    # Concordance correlation coef between Literature ages and DBs ages.
    lit = list(p99[6]) + list(p00[6]) + list(c06[6]) + list(g10[6])
    dbs = list(p99[2]) + list(p00[2]) + list(c06[2]) + list(g10[2])
    print 'Lit vs DBs age conc corr coef:', ccc(lit, dbs)
    print 'Lit vs DBs Pearson corr coef:', np.corrcoef(lit, dbs)[0, 1]

    labels = [['P99', 'C06', 'G10'], ['SMC', 'LMC'],
              ['P99', 'P00', 'C06', 'G10']]
    mark = [['>', 'v', '<'], ['o', '*'], ['>', '^', 'v', '<']]
    cols = [['chocolate', 'c', 'g'], ['m', 'b'], ['chocolate', 'r', 'c', 'g']]

    # Define names of arrays being plotted.
    x_lab = ['$\Delta log(age/yr)_{ASteCA-DB}$',
             '$\Delta log(age/yr)_{ASteCA-lit}$',
             '$E(B-V)_{ASteCA}$', '$E(B-V)_{lit}$',
             '$\Delta log(age/yr)_{lit-DB}$', '$log(age/yr)_{lit}$']
    y_lab = ['$\Delta E(B-V)_{ASteCA-DB}$', '$\Delta E(B-V)_{ASteCA-lit}$',
             '$E(B-V)_{DB}$', '$E(B-V)_{lit}$', '$\Delta E(B-V)_{lit-DB}$',
             '$log(age/yr)_{DB}$']
    xmm, ymm = [-1.5, 1.5, -0.019, 0.31], [-0.19, 0.19]

    # The arbitrary size is set so that the fig sizes are equivalent to that
    # of the cross_match plots.
    fig = plt.figure(figsize=(13.77, 25.2))
    gs = gridspec.GridSpec(4, 2)

    cross_match_lst = [
        # Age 1:1, literature vs databases.
        [gs, 0, 5.8, 10.5, 5.8, 10.5, x_lab[5], y_lab[5],
            age_lit_DB_data, labels[2], mark[2], cols[2], []],

        # Extinction 1:1, literature vs DBs.
        [gs, 2, xmm[2], xmm[3], xmm[2], xmm[3], x_lab[3], y_lab[2],
            ext_lit_DB_data, labels[0], mark[0], cols[0], []],
        # Age vs ext diff for literature vs DBs.
        [gs, 3, xmm[0], xmm[1], ymm[0], ymm[1], x_lab[4], y_lab[4],
            lit_db_data, labels[0], mark[0], cols[0], kde_cont[2]],

        # Extinction 1:1, ASteCA vs databases.
        [gs, 4, xmm[2], xmm[3], xmm[2], xmm[3], x_lab[2], y_lab[2],
            ext_DB_data, labels[0], mark[0], cols[0], []],
        # Age vs ext diff for ASteCA vs databases.
        [gs, 5, xmm[0], xmm[1], ymm[0], ymm[1], x_lab[0], y_lab[0],
            db_data, labels[0], mark[0], cols[0], kde_cont[0]],

        # Extinction 1:1, ASteCA vs literature.
        [gs, 6, xmm[2], xmm[3], xmm[2], xmm[3], x_lab[2], y_lab[3],
            ext_lit_data, labels[1], mark[1], cols[1], []],
        # Age vs ext diff for ASteCA vs literature.
        [gs, 7, xmm[0], xmm[1], ymm[0], ymm[1], x_lab[1], y_lab[1],
            lit_data, labels[1], mark[1], cols[1], kde_cont[1]]
    ]

    for pl_params in cross_match_lst:
        cross_match_age_ext_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/cross_match_age_ext.png', dpi=300,
                bbox_inches='tight')


def pl_DBs_ASteCA_CMDs(pl_params):
    '''
    Star's membership probabilities on cluster's photom diagram.
    '''
    gs, i, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax, cl, db,\
        gal, cl_reg_fit, cl_reg_no_fit, lit_isoch, asteca_isoch, db_z, db_a,\
        db_e, db_d, as_z, as_a, as_e, as_d = pl_params

    # DB isoch fit.
    ax = plt.subplot(gs[i])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Add text box.
    text = '$' + cl + '-' + gal + '\,({})$'.format(db)
    ob1 = offsetbox.AnchoredText(text, loc=1, prop=dict(size=11))
    ob1.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax.add_artist(ob1)
    text1 = r'$z={}$'.format(db_z)
    text2 = '\n' + r'$log(age/yr)={}$'.format(float(db_a))
    text3 = '\n' + r'$E_{{(B-V)}}={}$'.format(db_e)
    text4 = '\n' + r'$dm={}$'.format(db_d)
    text = text1 + text2 + text3 + text4
    ob2 = offsetbox.AnchoredText(text, loc=2, prop=dict(size=11))
    ob2.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax.add_artist(ob2)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Plot grid.
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')
    col_select_fit, c_iso = '#4682b4', 'r'
    # Plot stars used in the best fit process.
    cl_reg_x = cl_reg_fit[0] + cl_reg_no_fit[0]
    cl_reg_y = cl_reg_fit[1] + cl_reg_no_fit[1]
    plt.scatter(cl_reg_x, cl_reg_y, marker='o',
                c=col_select_fit, s=40, cmap=cm, lw=0.5, zorder=4)
    # Plot isochrone.
    plt.plot(lit_isoch[0], lit_isoch[1], c=c_iso, lw=1.2, zorder=5)

    # ASteCA isoch fit.
    ax = plt.subplot(gs[i + 1])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Add text box.
    text = '$' + cl + '-' + gal + '\,(ASteCA)$'
    ob1 = offsetbox.AnchoredText(text, loc=1, prop=dict(size=11))
    ob1.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax.add_artist(ob1)
    text1 = r'$z={}$'.format(as_z)
    text2 = '\n' + r'$log(age/yr)={}$'.format(as_a)
    text3 = '\n' + r'$E_{{(B-V)}}={}$'.format(as_e)
    text4 = '\n' + r'$dm={}$'.format(as_d)
    text = text1 + text2 + text3 + text4
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=11))
    ob.patch.set(boxstyle='square,pad=-0.2', alpha=0.75)
    ax.add_artist(ob)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # Plot grid.
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    # This reversed colormap means higher prob stars will look redder.
    cm = plt.cm.get_cmap('RdYlBu_r')

    # Get extreme values for colorbar.
    lst_comb = cl_reg_fit[2] + cl_reg_no_fit[2]
    v_min_mp, v_max_mp = round(min(lst_comb), 2), round(max(lst_comb), 2)
    col_select_fit, col_select_no_fit, c_iso = cl_reg_fit[2], \
        cl_reg_no_fit[2], 'g'
    # Plot stars *not* used in the best fit process.
    plt.scatter(cl_reg_no_fit[0], cl_reg_no_fit[1], marker='o',
                c=col_select_no_fit, s=35, cmap=cm, lw=0.5, alpha=0.5,
                vmin=v_min_mp, vmax=v_max_mp, zorder=2)
    # Plot stars used in the best fit process.
    plt.scatter(cl_reg_fit[0], cl_reg_fit[1], marker='o',
                c=col_select_fit, s=40, cmap=cm, lw=0.5, vmin=v_min_mp,
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
                cl_reg_fit, cl_reg_no_fit, lit_isoch, asteca_isoch, db_z, \
                db_a, db_e, db_d, as_z, as_a, as_e, as_d = cl_data

            db_sat_cmd_lst.append(
                [gs, i, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, '(C-T_1)',
                    'T_1', cl, db, gal, cl_reg_fit, cl_reg_no_fit, lit_isoch,
                    asteca_isoch, db_z, db_a, db_e, db_d, as_z, as_a, as_e,
                    as_d])

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
        fig_name = 'figures/DB_fit/' + db + '_VS_asteca_' + str(k) + '.png'
        plt.savefig(fig_name, dpi=150, bbox_inches='tight')


def pl_errors(pl_params):
    '''
    '''
    gs, i, xmin, xmax, ymin, ymax, x, y, z, rad, x_lab, y_lab =\
        pl_params
    siz = np.asarray(rad) * 1.

    xy_font_s = 16
    ax = plt.subplot(gs[i])
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
    if x_lab == '$[Fe/H]_{ASteCA}$':
        # 1% of axis ranges.
        ax_ext = (xmax - xmin) * 0.02
    elif x_lab == '$M_{\odot;\,ASteCA}$':
        ax_ext = (xmax - xmin) * 0.01
    else:
        ax_ext = (xmax - xmin) * 0.025
    # Add random scatter.
    r_x = x + np.random.uniform(-ax_ext, ax_ext, len(x))
    SC = plt.scatter(r_x, y, marker='o', c=z, edgecolor='k', s=siz,
                     cmap=cm, lw=0.25, zorder=4) # vmin=0., vmax=1., 
    # # Text box.
    # ob = offsetbox.AnchoredText(gal_name, loc=4, prop=dict(size=xy_font_s))
    # ob.patch.set(alpha=0.85)
    # ax.add_artist(ob)
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    zpad = 10
    cbar.set_label(r'$CI$', fontsize=xy_font_s - 2, labelpad=zpad)
    cbar.ax.tick_params(labelsize=10)


def make_errors_plots(in_params):
    '''
    '''

    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, marr, msigma,\
        rarr, cont_ind, kde_prob, phot_disp = [
            in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
                                   'esigma', 'darr', 'dsigma', 'marr',
                                   'msigma', 'rarr', 'cont_ind', 'kde_prob',
                                   'phot_disp']]

    p_disp = phot_disp[0] + phot_disp[1]
    ci = cont_ind[0] + cont_ind[1]
    probs = kde_prob[0] + kde_prob[1]
    r_arr = rarr[0][0] + rarr[1][0]
    z_arr = zarr[0][0] + zarr[1][0]
    z_sigma = zsigma[0][0] + zsigma[1][0]
    a_arr = aarr[0][0] + aarr[1][0]
    a_sigma = asigma[0][0] + asigma[1][0]
    e_arr = earr[0][0] + earr[1][0]
    e_sigma = esigma[0][0] + esigma[1][0]
    d_arr = darr[0][0] + darr[1][0]
    d_sigma = dsigma[0][0] + dsigma[1][0]
    m_arr = marr[0][0] + marr[1][0]
    m_sigma = zsigma[0][0] + msigma[1][0]

    # Order lists to put min rad values on top.
    ord_r, ord_z, ord_zs, ord_a, ord_as, ord_e, ord_es, ord_d, ord_ds, ord_m,\
        ord_ms, ord_ci, ord_prob, ord_p_disp = map(list, zip(*sorted(zip(
            r_arr, z_arr, z_sigma, a_arr, a_sigma, e_arr, e_sigma, d_arr,
            d_sigma, m_arr, m_sigma, ci, probs, p_disp), reverse=True)))

    # Select colorbar parameter.
    ord_X = np.array(ord_prob) / np.array(ord_p_disp)

    fig = plt.figure(figsize=(10, 20))
    gs = gridspec.GridSpec(5, 1)

    errors_lst = [
        [gs, 0, -2.4, 0.11, -0.03, 2.1, ord_z, ord_zs, ord_X, ord_r,
            '$[Fe/H]_{ASteCA}$', '$e_{[Fe/H]}$'],
        [gs, 1, 6.51, 10.1, -0.03, 1.1, ord_a, ord_as, ord_X, ord_r,
            '$log(aye/yr)_{ASteCA}$', '$e_{log(aye/yr)}$'],
        [gs, 2, -0.02, 0.32, -0.01, 0.11, ord_e, ord_es, ord_X, ord_r,
            '$E(B-V)_{ASteCA}$', '$e_{E(B-V)}$'],
        [gs, 3, 18.28, 19.19, 0.007, 0.083, ord_d, ord_ds, ord_X, ord_r,
            '$(m-M)_{\circ;\,ASteCA}$', '$e_{(m-M)_{\circ}}$'],
        [gs, 4, -210, 30000, -210, 4450, ord_m, ord_ms, ord_X, ord_r,
            '$M_{\odot;\,ASteCA}$', '$e_{M_{\odot}}$']
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

    gs, i, age_vals, met_weighted, age_gyr, amr_lit, zarr, x_lab,\
        y_lab = pl_params

    xy_font_s = 16
    ax = plt.subplot(gs[i])
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()
    plt.xlim(-0.02, 8.4)
    if i == 0:
        plt.ylim(-2.45, 0.4)
        ax.set_xticklabels([])
        col, leg = ['r', 'b'], ['SMC', 'LMC']
        for k in [1, 0]:
            # ASteCA values.
            plt.plot(age_vals[k], met_weighted[k][0], c=col[k],
                     label=leg[k] + ' (ASteCA)', zorder=3)
            # Introduce random scatter in Age (Gyr).
            # 2% of axis ranges.
            ax_ext = max(age_gyr[k][0]) * 0.02
            # Add randoms scatter.
            rs_x = age_gyr[k][0] + np.random.uniform(-ax_ext, ax_ext,
                                                     len(age_gyr[k][0]))
            plt.scatter(rs_x, zarr[k][0], marker='*', s=25, edgecolors=col[k],
                        facecolor='none', lw=0.4, label=leg[k], zorder=1)
            # ASteCA 1 sigma error regions.
            y_err_min = np.array(met_weighted[k][0]) -\
                np.array(met_weighted[k][1])
            y_err_max = np.array(met_weighted[k][0]) +\
                np.array(met_weighted[k][1])
            plt.fill_between(age_vals[k], y_err_min, y_err_max, alpha=0.1,
                             color=col[k])
        # Legend.
        leg0 = plt.legend(loc='lower right', handlelength=2.5, scatterpoints=1,
                          fontsize=xy_font_s - 8)
        leg0.get_frame().set_alpha(0.85)

    # Literature values.
    elif i == 1:
        ax.set_xticklabels([])
        plt.ylim(-1.23, -0.06)
        ax.set_title("LMC", x=0.5, y=0.92, fontsize=xy_font_s - 4)
        col = ['m', 'k', 'g', 'c', 'y']
        c_dash = [[8, 4], [8, 4, 2, 4], [2, 2], [8, 4, 2, 4, 2, 4], [8, 4]]
        amr_lab = ['PT98', 'PG03', 'C08', 'HZ09', 'R12']
        for j, amr in enumerate(amr_lit):
            plt.plot(amr[0], amr[1], color=col[j], label=amr_lab[j],
                     dashes=c_dash[j], lw=1.5, zorder=1)
        # ASteCA values.
        plt.plot(age_vals[1], met_weighted[1][0], c='b', label='ASteCA',
                 zorder=3)
        # Legend.
        leg2 = plt.legend(loc='lower left', handlelength=3.5, scatterpoints=1,
                          fontsize=xy_font_s - 8)
        leg2.get_frame().set_alpha(0.85)
    elif i == 2:
        plt.ylim(-1.39, -0.39)
        plt.xlabel(x_lab, fontsize=xy_font_s)
        ax.set_title("SMC", x=0.5, y=0.92, fontsize=xy_font_s - 4)
        col = ['m', 'k', 'g', 'c', 'y', 'y', '#b22222', '#b22222']
        c_dash = [[8, 4], [8, 4, 2, 4], [2, 2], [8, 4, 2, 4, 2, 4],
                  [8, 4, 2, 4, 2, 4], [8, 4], [2, 2], [8, 4, 2, 4]]
        amr_lab = ['PT98', 'PG03', 'HZ04', 'N09', 'TB09-1', 'TB09-2',
                   'C13-B', 'C13-C']
        for j, amr in enumerate(amr_lit):
            plt.plot(amr[0], amr[1], color=col[j], label=amr_lab[j],
                     dashes=c_dash[j], lw=1.5, zorder=1)
        # ASteCA values.
        plt.plot(age_vals[0], met_weighted[0][0], c='r', label='ASteCA',
                 zorder=3)
        # Legend.
        leg1 = plt.legend(loc='lower left', handlelength=3.5, scatterpoints=1,
                          fontsize=xy_font_s - 8)
        leg1.get_frame().set_alpha(0.85)


def make_amr_plot(in_params, amr_lit):
    '''
    Make age-metallicity relation plot for both galaxies.
    '''

    zarr, zsigma, aarr, asigma = [in_params[_] for _ in ['zarr', 'zsigma',
                                                         'aarr', 'asigma']]

    # First index k indicates the galaxy (0 for SMC, 1 for LMC), the second
    # index 0 indicates ASteCA values.
    # k=0 -> SMC, k=1 ->LMC
    age_gyr, age_vals, met_weighted = [[], []], [[], []], [[], []]
    for k in [0, 1]:
        # Age in Gyrs.
        age_gyr[k] = [10 ** (np.asarray(aarr[k][0]) - 9),
                      np.asarray(asigma[k][0]) * np.asarray(aarr[k][0]) *
                      np.log(10) / 5.]
        # Weighted metallicity values for an array of ages.
        # Max limit on very large met errors.
        zsig = [min(2., _) for _ in zsigma[k][0]]
        age_vals[k], met_weighted[k] = age_met_rel(
            age_gyr[k][0], age_gyr[k][1], zarr[k][0], zsig)

    fig = plt.figure(figsize=(5.25, 13.5))
    gs = gridspec.GridSpec(3, 1)

    amr_lit_smc, amr_lit_lmc = amr_lit

    amr_lst = [
        [gs, 0, age_vals, met_weighted, age_gyr, [], zarr, '', '$[Fe/H]$'],
        [gs, 1, age_vals, met_weighted, age_gyr, amr_lit_lmc, zarr,
         '', '$[Fe/H]$'],
        [gs, 2, age_vals, met_weighted, age_gyr, amr_lit_smc, zarr,
         '$Age\,(Gyr)$', '$[Fe/H]$']
    ]

    for pl_params in amr_lst:
        pl_amr(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/AMR_asteca.png', dpi=300, bbox_inches='tight')


def pl_angles(in_pars):
    '''
    '''
    gs, i, xlab, ylab, ang_pars = in_pars
    xmin, xmax, ymin, ymax, xi, yi, zi, mean_pos, i_pa_std, e_w, e_h, theta = \
        ang_pars

    # Define gridspec ranges.
    x_gdsp = [[0, 2], [2, 4], [4, 6]]
    gal_id = ['SMC', 'LMC']
    if i in [0, 1, 2]:
        y1, y2 = 0, 1
        x1, x2 = x_gdsp[i]
        gal_name = gal_id[0]
    elif i in [3, 4, 5]:
        y1, y2 = 1, 2
        x1, x2 = x_gdsp[i-3]
        gal_name = gal_id[1]

    ax = plt.subplot(gs[y1:y2, x1:x2])
    xy_font_s = 12
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.tick_params(axis='both', which='major', labelsize=xy_font_s - 3)
    # Set axis labels
    plt.xlabel(xlab, fontsize=xy_font_s)
    plt.ylabel(ylab, fontsize=xy_font_s)

    t = ['M1: Match deprojected distances',
         'M2: Minimize distance to plane (fixed)',
         'M3: Minimize distance to plane (free)']
    if i in [0, 1, 2]:
        ax.set_xticklabels([])
        plt.title(t[i], fontsize=11)
    if i in [1, 2, 4, 5]:
        ax.set_yticklabels([])

    sub_idx = ['1', '3', '5', '2', '4', '6']
    # Set minor ticks
    ax.minorticks_on()
    if i in [0, 3]:
        cm = plt.cm.get_cmap('RdBu_r')
    else:
        cm = plt.cm.get_cmap('RdBu')
    # Only draw units on axis (ie: 1, 2, 3)
    ax.xaxis.set_major_locator(MultipleLocator(20.))
    ax.yaxis.set_major_locator(MultipleLocator(20.))
    # Plot density map.
    SC = plt.imshow(np.rot90(zi), origin='upper',
                    extent=[xi.min(), xi.max(), yi.min(), yi.max()],
                    cmap=cm)
    # Plot contour curves.
    curv_num = 150
    plt.contour(np.rot90(zi), curv_num, colors='k', linewidths=0.2,
                origin='upper',
                extent=[xi.min(), xi.max(), yi.min(), yi.max()])
    plt.scatter(mean_pos[0], mean_pos[1], c='w', s=45, zorder=3)
    # Standard dev ellipse.
    ellipse = Ellipse(xy=(mean_pos[0], mean_pos[1]),
                      width=e_w, height=e_h, angle=theta,
                      edgecolor='w', fc='None', lw=0.85, zorder=3)
    ax.add_patch(ellipse)
    # Text box.
    text1 = r'$\Theta_{{{}}}={:.0f}^{{\circ}}\pm{:.0f}$'.format(
        sub_idx[i], mean_pos[1], i_pa_std[1])
    text2 = r'$i_{{{}}}={:.0f}^{{\circ}}\pm{:.0f}$'.format(
        sub_idx[i], mean_pos[0], i_pa_std[0])
    text = text1 + '\n' + text2
    ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Gal name.
    ob = offsetbox.AnchoredText(gal_name, loc=1, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    cb_lab = [r'$ccc$', r'$\overline{{|d_{{p}}}|}$',
              r'$\overline{{|d_{{p}}}|}$', r'$ccc$',
              r'$\overline{{|d_{{p}}}|}$', r'$\overline{{|d_{{p}}}|}$']
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    cbar.set_label(cb_lab[i], fontsize=xy_font_s, labelpad=4, y=0.5)
    cbar.ax.tick_params(labelsize=xy_font_s - 3)
    if i not in [0, 3]:
        cbar.ax.invert_yaxis()
    # Square
    ax.set_aspect(aspect='auto')


def diag_plots(in_pars):
    """
    """
    gs, i, xlab, ylab, ang_pars = in_pars

    x_gdsp = [[[0, 1], [2, 3], [4, 5]], [[1, 2], [3, 4],  [5, 6]]]
    y1, y2 = 2, 3
    x11, x12 = x_gdsp[0][i]
    x21, x22 = x_gdsp[1][i]

    xy_font_s = 12
    cm = plt.cm.get_cmap('RdYlBu_r')

    # Unpack
    xmin, xmax, ymin, ymax, xi, yi, zi, mean_pos, i_pa_std, e_w, e_h, theta = \
        ang_pars[i]
    # Left plot.
    right_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2, width_ratios=[1, 6], subplot_spec=gs[y1:y2, x11:x12],
        wspace=0.)
    axl = plt.subplot(right_gs[1])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlab, fontsize=xy_font_s)
    plt.tick_params(axis='both', which='major', labelsize=xy_font_s - 3)
    axl.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
             zorder=1)
    if i == 0:
        axl.set_ylabel(ylab, fontsize=xy_font_s)
    else:
        axl.set_yticklabels([])
    axl.tick_params(axis='both', which='major', labelsize=9)
    plt.plot([xmin, xmax], [ymin, ymax], 'k', ls='--')
    plt.scatter(xi, yi, marker='o', c=zi, edgecolor='k', s=25,
                cmap=cm, lw=0.25, zorder=4)
    # Gal name.
    text1 = r'$[\Theta_{{{}}},i_{{{}}}]$'.format(i*2+1, i*2+1)
    text = 'SMC\n' + text1
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    axl.add_artist(ob)
    # Text box.
    text2 = r'$ccc={:.2f}$'.format(mean_pos[i][0])
    text3 = r'$\overline{{|d_{{p}}|}}={:.2f}$'.format(mean_pos[i][1])
    text = text2 + '\n' + text3
    ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=xy_font_s - 1))
    ob.patch.set(alpha=0.85)
    axl.add_artist(ob)

    # Unpack
    xmin, xmax, ymin, ymax, xi, yi, zi, mean_pos, i_pa_std, e_w, e_h, theta = \
        ang_pars[i+3]
    # Right pot.
    right_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2, width_ratios=[6, 1], subplot_spec=gs[y1:y2, x21:x22],
        wspace=0.)
    axr = plt.subplot(right_gs[0])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel(xlab, fontsize=xy_font_s)
    plt.tick_params(axis='both', which='major', labelsize=xy_font_s - 3)
    axr.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
             zorder=1)
    axr.set_yticklabels([])
    plt.plot([xmin, xmax], [ymin, ymax], 'k', ls='--')
    SC = plt.scatter(xi, yi, marker='o', c=zi, edgecolor='k', s=25,
                     cmap=cm, lw=0.25, zorder=4)
    # Gal name.
    text1 = r'$[\Theta_{{{}}},i_{{{}}}]$'.format(i*2+2, i*2+2)
    text = 'LMC\n' + text1
    ob = offsetbox.AnchoredText(text, loc=2, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    axr.add_artist(ob)
    # Text box.
    text2 = r'$ccc={:.2f}$'.format(mean_pos[i][0])
    text3 = r'$\overline{{|d_{{p}}|}}={:.2f}$'.format(mean_pos[i][1])
    text = text2 + '\n' + text3
    ob = offsetbox.AnchoredText(text, loc=4, prop=dict(size=xy_font_s - 1))
    ob.patch.set(alpha=0.85)
    axr.add_artist(ob)
    # Position colorbar.
    the_divider = make_axes_locatable(axr)
    color_axis = the_divider.append_axes("right", size="5%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    cbar.set_label(r'$log(aye/yr)_{ASteCA}$', fontsize=xy_font_s, labelpad=4,
                   y=0.5)
    cbar.ax.tick_params(labelsize=xy_font_s - 3)


def make_angles_plot(gal_str_pars):
    '''
    Plot inclination and position angles density maps for both galaxies.
    '''

    fig = plt.figure(figsize=(14.25, 13.5))
    gs = gridspec.GridSpec(3, 6)

    # New gridspecs for bottom rectangular plots.
    # ***Values selected by hand***
    gs2 = gridspec.GridSpec(3, 6)
    gs2.update(wspace=0.0, bottom=0.029, left=0.031, right=0.95)
    gs3 = gridspec.GridSpec(3, 6)
    gs3.update(wspace=0.0, bottom=0.029, left=0.043, right=0.965)
    gs4 = gridspec.GridSpec(3, 6)
    gs4.update(wspace=0.0, bottom=0.029, left=0.05, right=0.976)

    xlab = ['', r'Inclination ($i^{\circ}$)', r'$d_{GC}\,(Kpc)$']
    ylab = ['', r'Position angle ($\Theta^{\circ}$)',
            r'$d_{[\Theta_m,i_m]}\,(Kpc)$']

    str_lst = [
        # SMC dens map
        [gs, 0, xlab[0], ylab[1], gal_str_pars[0][0]],
        [gs, 1, xlab[0], ylab[0], gal_str_pars[0][1]],
        [gs, 2, xlab[0], ylab[0], gal_str_pars[0][2]],
        # LMC dens map
        [gs, 3, xlab[1], ylab[1], gal_str_pars[0][3]],
        [gs, 4, xlab[1], ylab[0], gal_str_pars[0][4]],
        [gs, 5, xlab[1], ylab[0], gal_str_pars[0][5]]
    ]

    for pl_params in str_lst:
        pl_angles(pl_params)

    str_lst = [
        # S/LMC CCC 1:1 diagonal diagrams.
        [gs2, 0, xlab[2], ylab[2], gal_str_pars[1]],
        [gs3, 1, xlab[2], ylab[0], gal_str_pars[1]],
        [gs4, 2, xlab[2], ylab[0], gal_str_pars[1]]
    ]
    for pl_params in str_lst:
        diag_plots(pl_params)

    fig.tight_layout()
    # plt.savefig('figures/MCs_deproj_dist_angles.png', dpi=300,
    #             bbox_inches='tight')
    plt.savefig('MCs_deproj_dist_angles.png', dpi=150, bbox_inches='tight')


def pl_rho_var(in_pars):
    '''
    '''
    # gs, i, xlab, ylab, gal_name, ang_pars = in_pars
    # r_min, N_clust, inc_best, inc_std, pa_best, pa_std, ccc_best =\
    #     zip(*ang_pars)

    # if i in [0, 1]:
    #     y1, y2, y3 = zip(*pa_best)
    #     e_y = pa_std
    # elif i in [2, 3]:
    #     y1, y2, y3 = zip(*inc_best)
    #     e_y = inc_std
    # ccc_cols = zip(*ccc_best)

    gs, i, xlab, ylab, gal_name, v_min, v_max, ang_pars = in_pars
    r_min, N_clust, inc_mean, inc_std, pa_mean, pa_std, ccc_mean =\
        zip(*ang_pars)

    if i in [0, 1]:
        yi, e_yi = pa_mean, pa_std
    elif i in [2, 3]:
        yi, e_yi = inc_mean, inc_std

    right_gs = gridspec.GridSpecFromSubplotSpec(
        1, 2, width_ratios=[40, 1], subplot_spec=gs[i], wspace=0.05)
    ax = plt.subplot(right_gs[0])
    if i in [1, 3]:
        # Only draw subplot containing colorbar for these subplots.
        color_axis = plt.subplot(right_gs[1])

    ax2 = ax.twiny()
    xy_font_s = 12
    ax.set_xlim(-0.5, max(r_min)+0.5)
    if i in [0, 1]:
        plt.ylim(0., 230.)
    else:
        plt.ylim(-89., 89.)
    # Increase font since it's LaTeX and it looks small.
    ax.set_xlabel(xlab, fontsize=xy_font_s+3)
    ax.set_ylabel(ylab, fontsize=xy_font_s)
    # Set minor ticks
    ax.minorticks_on()
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.3,
            zorder=1)
    # cm = plt.cm.get_cmap('RdBu_r')
    # m = ['>', 'v', '^']
    # for j, y in enumerate([y1, y2, y3]):
    #     sx = np.array(r_min) + j*0.1
    #     plt.errorbar(sx, y, yerr=e_y, ls='none', color='k',
    #                  elinewidth=0.5, zorder=3)
    #     SC = plt.scatter(
    #         sx, y, marker=m[j], c=ccc_cols[j], edgecolor='k',
    #         s=80, cmap=cm, lw=0.2, vmin=0., vmax=1., zorder=4)

    plt.errorbar(r_min, yi, yerr=e_yi, ls='none', color='k', elinewidth=0.5,
                 zorder=3)
    cm = plt.cm.get_cmap('RdBu')
    SC = plt.scatter(r_min, yi, marker='o', c=ccc_mean, edgecolor='k', s=80,
                     cmap=cm, lw=0.2, vmin=v_min, vmax=v_max, zorder=4)
    # Set font size for the three axes.
    ax.set_xticklabels(ax.get_xticks(), fontsize=xy_font_s-3)
    ax.set_yticklabels(ax.get_yticks(), fontsize=xy_font_s-3)
    # Set range, ticks, and label for the second x axis.
    # Second x axis.
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(r_min)
    if i in [0, 1]:
        ax.set_xticklabels([])
        ax2.set_xticklabels(N_clust, fontsize=xy_font_s-3)
        ax2.set_xlabel(r"$N_{clusters}$", fontsize=xy_font_s+3, labelpad=10.5)
    else:
        ax2.set_xticklabels([])
    if i in [1, 3]:
        ax.set_yticklabels([])
    # Filling that shows the range of values in the literature.
    f = [[120., 161.6], [114.2, 175.], [38., 83.], [16., 39.7]]
    plt.fill_between([-9., 9.], f[i][0], f[i][1], color='grey', alpha='0.2',
                     zorder=1)
    # Gal name.
    ob = offsetbox.AnchoredText(gal_name, loc=1, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Position colorbar.
    if i in [1, 3]:
        cbar = plt.colorbar(SC, cax=color_axis)
        cb_lab = r'$\overline{{|d_{{p}}}|}$'
        cbar.set_label(cb_lab, fontsize=xy_font_s, labelpad=6, y=0.5)
        cbar.ax.tick_params(labelsize=xy_font_s - 3)
        cbar.ax.invert_yaxis()


def make_rho_min_plot(rho_plot_pars):
    '''
    Plot variation of the inclination and position angles with the selected
    minimum projected angular density value, for both galaxies.
    '''

    fig = plt.figure(figsize=(11.5, 11))
    gs = gridspec.GridSpec(2, 2)

    labels = [r'Inclination ($i^{\circ}$)',
              r'Position angle ($\Theta^{\circ}$)', r'$\rho_{min}$']

    # Extract min and max CCC values.
    ccc_mean = zip(*rho_plot_pars[0])[-1] + zip(*rho_plot_pars[1])[-1]
    v_min, v_max = min(ccc_mean), max(ccc_mean)

    str_lst = [
        # [gs, 0, '', labels[1], 'SMC', rho_plot_pars[0]],
        # [gs, 1, '', '', 'LMC', rho_plot_pars[1]],
        # [gs, 2, labels[2], labels[0], 'SMC', rho_plot_pars[0]],
        # [gs, 3, labels[2], '', 'LMC', rho_plot_pars[1]]
        [gs, 0, '', labels[1], 'SMC', v_min, v_max, rho_plot_pars[0]],
        [gs, 1, '', '', 'LMC', v_min, v_max, rho_plot_pars[1]],
        [gs, 2, labels[2], labels[0], 'SMC', v_min, v_max, rho_plot_pars[0]],
        [gs, 3, labels[2], '', 'LMC', v_min, v_max, rho_plot_pars[1]]
    ]

    for pl_params in str_lst:
        pl_rho_var(pl_params)

    fig.tight_layout()
    # plt.savefig('figures/MCs_angles_var_w_rho.png', dpi=300,
    #             bbox_inches='tight')
    plt.savefig('MCs_angles_var_w_rho.png', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    rho_plot_pars = [[[0.0, 89, -71.135672033726252, 76.714077100490456, 50.861586766820743, 38.32345848712243, 1.0375530985205295], [0.5, 87, 17.43788156252084, 73.536746084355187, 109.93369061718404, 42.541549144111976, 1.0439379418678658], [1.0, 49, 42.428234431070344, 60.368943745258363, 99.615689627560769, 40.121405771967183, 1.1276678381517391], [1.5, 31, 27.557379861729999, 52.150867035281152, 90.767386854974731, 45.102293883944391, 1.252288517832606], [2.0, 24, 8.0952301240493654, 52.61098231435377, 90.615582049773352, 31.678776820572025, 1.3353770118232866], [2.5, 18, 5.4514194107377802, 56.780473636754493, 102.40456209615986, 27.919100072786811, 1.2857535342829622], [3.0, 15, 5.0301541503560214, 59.80224224607845, 104.85266326575095, 26.974447513509798, 1.2289344950057322], [3.5, 10, -44.197611422586554, 50.983468462530332, 71.904368461075123, 41.866828318076585, 1.3100775535144198], [4.0, 7, -23.129139062801297, 45.369018519254752, 53.963591271446397, 33.610281172748685, 1.1826356424601696]], [[0.0, 150, 25.323146024158603, 52.300330951873136, 159.48214375936985, 22.257460522931943, 1.0535850336215742], [0.5, 146, 24.650741770469207, 49.246715966266883, 159.33448682046139, 20.883458107436041, 1.050688651127091], [1.0, 137, 23.384853209322415, 44.265431722849769, 159.18858726876783, 22.51689391895232, 1.0339773090191755], [1.5, 127, 21.796480188728726, 45.710349151428503, 159.86272280425098, 22.071000347706402, 1.0012459139603553], [2.0, 87, 15.209619738317508, 23.679747403043208, 170.03645919238599, 42.334219885888693, 0.93716995629040944], [2.5, 62, 15.121132982345188, 22.229900310451505, 168.16235535529694, 50.165792260206025, 0.9134585017692104], [3.0, 42, 10.972077150955121, 17.620401283815998, 146.03220223934318, 48.249131637487814, 0.93920061922466835], [3.5, 39, 11.866543841572947, 19.405873344557289, 147.82091738601366, 50.253661289425295, 0.95078505321284201], [4.0, 28, 12.562343206761108, 16.859104998275843, 152.21216713712388, 47.121314606879963, 0.88553743394976048]]]
    make_rho_min_plot(rho_plot_pars)
