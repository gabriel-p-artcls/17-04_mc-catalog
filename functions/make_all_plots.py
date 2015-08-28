
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
# import statsmodels.api as sm
from scipy import stats

from functions.ra_dec_map import ra_dec_plots
from functions.kde_2d import kde_map
import functions.CMD_obs_vs_asteca as cmd


def as_vs_lit_plots(pl_params):
    '''
    Generate ASteCA vs literature values plots.
    '''
    gs, i, xmin, xmax, x_lab, y_lab, z_lab, xarr, xsigma, yarr, ysigma, \
        zarr, gal_name = pl_params

    xy_font_s = 18
    cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i], aspect='equal')
    # Different limits for \delta log(age) plot.
    if i != 4:
        plt.ylim(xmin, xmax)
        # 1:1 line
        plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')
    else:
        # ax = plt.subplot(gs[i], aspect='auto')
        plt.ylim(-2.4, 2.4)
        # 0 line
        plt.plot([5, 12], [0, 0], 'k', ls='--')
        plt.axhspan(-0.5, 0.5, facecolor='grey', alpha=0.5, zorder=1)

    plt.xlim(xmin, xmax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
            zorder=1)
    ax.minorticks_on()

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
    SC = plt.scatter(rs_x, rs_y, marker='o', c=zarr, s=70, lw=0.25, cmap=cm,
                     zorder=3)
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

    # # Plot IDs
    # if i == 4 and 
    # for label, x, y in zip(labels, x_data, y_data):
    #     plt.annotate(
    #         label, xy=(x, y),
    #         xytext=(4, 4),
    #         textcoords='offset points',
    #         ha='left', va='bottom',
    #         fontsize=20, color='#1A5233', zorder=4)

    # Text box.
    ob = offsetbox.AnchoredText(gal_name, loc=4, prop=dict(size=xy_font_s))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="5%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    zpad = 10 if z_lab == '$E_{(B-V)}$' else 5
    cbar.set_label(z_lab, fontsize=xy_font_s, labelpad=zpad)


def make_as_vs_lit_plot(galax, k, in_params):
    '''
    Prepare parameters and call function to generate ASteca vs literature
    SMC and LMC plots.
    '''

    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, rarr = \
        [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
                                'esigma', 'darr', 'dsigma', 'rarr']]

    # \delta log(age) as ASteCA - literature values.
    age_delta = np.array(aarr[k][0]) - np.array(aarr[k][1])

    # Generate ASteca vs literature plots.
    fig = plt.figure(figsize=(16, 25))  # create the top-level container
    # gs = gridspec.GridSpec(2, 4, width_ratios=[1, 0.35, 1, 0.35])
    gs = gridspec.GridSpec(4, 2)

    if galax == 'SMC':
        dm_min, dm_max = 18.62, 19.21
    else:
        dm_min, dm_max = 18.21, 18.79

    # For old runs where the dist mod range was large.
    # dm_min, dm_max = 17.8, 20.2

    as_lit_pl_lst = [
        [gs, 0, -2.4, 0.45, '$[Fe/H]_{ASteCA}$', '$[Fe/H]_{lit}$',
            '$log(age/yr)_{ASteCA}$', zarr[k][0], zsigma[k][0], zarr[k][1],
            zsigma[k][1], aarr[k][0], galax],
        [gs, 1, 5.8, 10.6, '$log(age/yr)_{ASteCA}$', '$log(age/yr)_{lit}$',
            '$E(B-V)_{ASteCA}$', aarr[k][0], asigma[k][0], aarr[k][1],
            asigma[k][1], earr[k][0], galax],
        [gs, 2, -0.04, 0.29, '$E(B-V)_{ASteCA}$', '$E(B-V)_{lit}$',
            '$log(age/yr)_{ASteCA}$', earr[k][0], esigma[k][0], earr[k][1],
            esigma[k][1], aarr[k][0], galax],
        [gs, 3, dm_min, dm_max, '$(m-M)_{0;\,ASteCA}$', '$(m-M)_{0;\,lit}$',
            '$log(age/yr)_{ASteCA}$', darr[k][0], dsigma[k][0], darr[k][1],
            dsigma[k][1], aarr[k][0], galax],
        # Asteca log(age) vs \delta log(age) with lit values.
        [gs, 4, 5.8, 10.6, '$log(age/yr)_{ASteCA}$', '$\Delta log(age/yr)$',
            '$E(B-V)_{ASteCA}$', aarr[k][0], asigma[k][0], age_delta,
            [], earr[k][0], galax]
        # [gs, 4, 1., 599., '$rad_{ASteCA} (px)$', '$rad_{lit} (px)$',
        #     '$log(age/yr)_{ASteCA}$', rarr[k][0], [], rarr[k][1], [],
        #     aarr[k][0], galax]
    ]
    #
    for pl_params in as_lit_pl_lst:
        as_vs_lit_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_vs_lit_' + galax + '.png', dpi=300)


def kde_plots(pl_params):
    '''
    Generate KDE plots.
    '''
    gs, i, x_lab, y_lab, xarr, xsigma, yarr, ysigma, x_rang, y_rang = pl_params

    ext = [x_rang[0], x_rang[1], y_rang[0], y_rang[1]]

    # Generate maps.
    z = kde_map(np.array(xarr), np.array(xsigma), np.array(yarr),
                np.array(ysigma), ext)

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
    plt.scatter(xarr + rs_x, yarr + rs_y, marker='*', color='#6b6868', s=40,
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
    age_rang, fe_h_rang, mass_rang = [6., 10.], [-2.4, 0.15], [-100., 30500.]
    if galax == 'SMC':
        E_bv_rang, dist_mod_rang = [-0.01, 0.15], [18.75, 19.25]
    else:
        E_bv_rang, dist_mod_rang = [-0.01, 0.3], [18.25, 18.75]

    kde_pl_lst = [
        [gs, 0, '$log(age/yr)$', '$[Fe/H]$', aarr[k][0], asigma[k][0],
            zarr[k][0], zsigma[k][0], age_rang, fe_h_rang],
        [gs, 1, '$log(age/yr)$', '$M\,(M_{\odot})$', aarr[k][0], asigma[k][0],
            marr[k][0], msigma[k][0], age_rang, mass_rang],
        [gs, 2, '$(m-M)_0$', '$E_{(B-V)}$', darr[k][0], dsigma[k][0],
            earr[k][0], esigma[k][0], dist_mod_rang, E_bv_rang],
        [gs, 3, '$M\,(M_{\odot})$', '$[Fe/H]$', marr[k][0], msigma[k][0],
            zarr[k][0], zsigma[k][0], mass_rang, fe_h_rang]
        # [gs, 4, '$log(age/yr)$', '$M\,(M_{\odot})$', aarr[k][0],
        # asigma[k][0], marr[k][0], msigma[k][0]],
        # [gs, 5, '$log(age/yr)$', '$M\,(M_{\odot})$', aarr[k][0],
        # asigma[k][0], marr[k][0], msigma[k][0]]
    ]
    #
    for pl_params in kde_pl_lst:
        kde_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_kde_maps_' + galax + '.png', dpi=300)


def make_ra_dec_plots(in_params, bica_coords):
    '''
    Prepare parameters and call function to generate RA vs DEC positional
    plots for the SMC and LMC.
    '''

    ra, dec, zarr, aarr, earr, darr, marr, rad_pc = [
        in_params[_] for _ in ['ra', 'dec', 'zarr', 'aarr', 'earr', 'darr',
                               'marr', 'rad_pc']]

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
        [fig, 321, ra, dec, bb_ra, bb_dec, zarr, rad_pc, '$[Fe/H]$'],
        [fig, 322, ra, dec, bb_ra, bb_dec, aarr, rad_pc, '$log(age/yr)$'],
        [fig, 323, ra, dec, bb_ra, bb_dec, earr, rad_pc, '$E_{(B-V)}$'],
        [fig, 324, ra, dec, bb_ra, bb_dec, darr, rad_pc, '$(m-M)_0$'],
        [fig, 325, ra, dec, bb_ra, bb_dec, marr, rad_pc, '$M\,(M_{\odot})$'],
        [fig, 326, ra, dec, bb_ra, bb_dec, rad_pc, rad_pc, '$r_{clust}\,[pc]$']
    ]

    for pl_params in ra_dec_pl_lst:
        ra_dec_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_RA_DEC.png', dpi=300)


def make_lit_ext_plot(in_params):
    '''
    Prepare parameters and call function to generate RA vs DEC positional
    plots for the SMC and LMC.
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
    xmin, xmax = -0.02, [0.15, 0.32, 0.4]
    x_lab = '$E(B-V)_{ASteCA}$'
    y_lab = ['$E(B-V)_{MCEV,\,closer}$', '$E(B-V)_{MCEV,\,max}$']
    z_lab = ['$log(age/yr)_{ASteCA}$', '$E(B-V)_{SF}$', '$dist\,(deg)$']

    fig = plt.figure(figsize=(16, 25))
    gs = gridspec.GridSpec(4, 2)

    ext_pl_lst = [
        # SMC
        [gs, 0, xmin, xmax[0], x_lab, y_lab[0], z_lab[2], ord_earr_smc_ast,
            ord_esig_smc_ast, ord_mcev_smc, ord_e_mcev_smc, ord_mcev_dist_smc,
            'SMC'],
        [gs, 1, xmin, xmax[2], x_lab, y_lab[1], z_lab[1], earr[0][0],
            esigma[0][0], ext_mcev[0][1], ext_mcev[0][2], ext_sf[0][0], 'SMC'],
        # LMC
        [gs, 2, xmin, xmax[1], x_lab, y_lab[0], z_lab[2], ord_earr_lmc_ast,
            ord_esig_lmc_ast, ord_mcev_lmc, ord_e_mcev_lmc, ord_mcev_dist_lmc,
            'LMC'],
        [gs, 3, xmin, xmax[2], x_lab, y_lab[1], z_lab[1], earr[1][0],
            esigma[1][0], ext_mcev[1][1], ext_mcev[1][2], ext_sf[1][0], 'LMC']
    ]

    for pl_params in ext_pl_lst:
        as_vs_lit_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_vs_lit_extin.png', dpi=300)


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
    plt.savefig('figures/as_integ_colors.png', dpi=300)


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
    plt.savefig('figures/concent_param.png', dpi=300)


def make_radius_plot(in_params):
    '''
    Plot radius (in pc) versus several parameters.
    '''

    zarr, zsigma, aarr, asigma, marr, msigma, n_memb, rad_pc, erad_pc, \
        r_core_pc, e_r_core = \
        [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma',
                                'marr', 'msigma', 'n_memb', 'rad_pc',
                                'erad_pc', 'r_core_pc', 'e_r_core']]

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
        plt.savefig('figures/as_rad_vs_params_' + gal_name + '.png', dpi=300)


def prob_vs_CI_plot(pl_params):
    '''
    Generate plots for KDE probabilities versus contamination indexes.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, xarr, yarr, zarr, \
        rad, gal_name = pl_params
    siz = np.asarray(rad) * 5

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

    zarr, zsigma, aarr, asigma, marr, msigma, rad_pc, kde_prob, cont_ind = \
        [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'marr',
                                'msigma', 'rad_pc', 'kde_prob', 'cont_ind']]

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
            kde_prob[0], zarr[0][0], rad_pc[0], 'SMC'],
        # LMC
        [gs, 2, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab[0], cont_ind[1],
            kde_prob[1], aarr[1][0], rad_pc[1], 'LMC'],
        [gs, 3, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab[1], cont_ind[1],
            kde_prob[1], zarr[1][0], rad_pc[1], 'LMC']
    ]

    for pl_params in prob_CI_pl_lst:
        prob_vs_CI_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_prob_vs_CI.png', dpi=300)


def make_dist_2_cents(in_params):
    '''
    Plot ASteCA distances to center of either MC galaxy.
    '''

    zarr, aarr, earr, marr, rad_pc, cont_ind, dist_cent = \
        [in_params[_] for _ in ['zarr', 'aarr', 'earr', 'marr',
                                'rad_pc', 'cont_ind', 'dist_cent']]

    # Define names of arrays being plotted.
    x_lab, yz_lab = '$dist_{center}\,[pc]$', \
        ['$log(age/yr)_{ASteCA}$', '$[Fe/H]_{ASteCA}$', '$M\,(M_{\odot})$',
            '$E(B-V)_{ASteCA}$']
    # xmax = 8000 if fixed dist_mod is used. Else xmax = 18000. Check
    # dist_2_cloud_center() in get_params.py to see which one is used.
    xmin, xmax = 0, 8000

    fig = plt.figure(figsize=(16, 25))
    gs = gridspec.GridSpec(4, 2)

    prob_CI_pl_lst = [
        # SMC
        [gs, 0, xmin, xmax, 6., 10.5, x_lab, yz_lab[0], yz_lab[1],
            dist_cent[0], aarr[0][0], zarr[0][0], rad_pc[0], 'SMC'],
        [gs, 1, xmin, xmax, -2.4, 0.4, x_lab, yz_lab[1], yz_lab[2],
            dist_cent[0], zarr[0][0], marr[0][0], rad_pc[0], 'SMC'],
        [gs, 2, xmin, xmax, 0., 30000, x_lab, yz_lab[2], yz_lab[3],
            dist_cent[0], marr[0][0], earr[0][0], rad_pc[0], 'SMC'],
        [gs, 3, xmin, xmax, -0.01, 0.11, x_lab, yz_lab[3], yz_lab[0],
            dist_cent[0], earr[0][0], aarr[0][0], rad_pc[0], 'SMC'],
        # LMC
        [gs, 4, xmin, xmax, 6., 10.5, x_lab, yz_lab[0], yz_lab[1],
            dist_cent[1], aarr[1][0], zarr[1][0], rad_pc[1], 'LMC'],
        [gs, 5, xmin, xmax, -2.4, 0.4, x_lab, yz_lab[1], yz_lab[2],
            dist_cent[1], zarr[1][0], marr[1][0], rad_pc[1], 'LMC'],
        [gs, 6, xmin, xmax, 0., 30000, x_lab, yz_lab[2], yz_lab[3],
            dist_cent[1], marr[1][0], earr[1][0], rad_pc[1], 'LMC'],
        [gs, 7, xmin, xmax, -0.01, 0.31, x_lab, yz_lab[3], yz_lab[0],
            dist_cent[1], earr[1][0], aarr[1][0], rad_pc[1], 'LMC']
    ]

    for pl_params in prob_CI_pl_lst:
        prob_vs_CI_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/as_dist_2_cent.png', dpi=300)


def cross_match_plot(pl_params):
    '''
    Generate plots for the cross-matched age and mass values.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, z_lab, indexes, labels, \
        mark, cols, text_box, databases = pl_params

    a, e_a, b, e_b = indexes

    xy_font_s = 16
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
            leg = plt.legend(loc='upper left', markerscale=1., scatterpoints=1,
                             fontsize=xy_font_s - 5)
            leg.get_frame().set_alpha(0.85)
    plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')  # 1:1 line
    # Text box.
    ob = offsetbox.AnchoredText(text_box, loc=4, prop=dict(size=xy_font_s - 3))
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
                'Low mass', 'Mass']

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
        # Mass cross_match (low mass)
        [gs, 4, xymin[1], xymax[1], xymin[1], xymax[1], x_lab[1], y_lab[1],
            z_lab[1], indexes[1], labels[4], mark[4], cols[4], text_box[4],
            databases[4]],
        # Mass cross_match (all)
        [gs, 5, xymin[1], xymax[2], xymin[1], xymax[2], x_lab[1], y_lab[1],
            z_lab[1], indexes[1], labels[4], mark[4], cols[4], text_box[5],
            databases[5]]
    ]

    for pl_params in cross_match_lst:
        cross_match_plot(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('figures/cross_match.png', dpi=300)


def cross_match_age_ext_plot(pl_params):
    '''
    Generate plots for the cross-matched age and mass values.
    '''
    gs, i, xmin, xmax, ymin, ymax, x_lab, y_lab, data, labels, mark, cols, \
        kde_cont = pl_params

    xy_font_s = 16

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
                     fontsize=xy_font_s - 5)
    # Set the alpha value of the legend.
    leg.get_frame().set_alpha(0.65)
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

    # P99 ASteCA minus database diffs.
    diffs_db_ages_p99 = np.array(p99[4]) - np.array(p99[2])
    # P99 liter minus database diffs.
    diffs_lit_db_ages_p99 = np.array(p99[6]) - np.array(p99[2])
    # Same for extinctions.
    diffs_db_exts_p99 = np.array(p99[13]) - np.array(p99[12])
    diffs_lit_db_exts_p99 = np.array(p99[15]) - np.array(p99[12])

    # P00 ASteCA minus database diffs.
    diffs_db_ages_p00 = np.array(p00[4]) - np.array(p00[2])

    # C06
    diffs_db_ages_c06 = np.array(c06[4]) - np.array(c06[2])
    diffs_lit_db_ages_c06 = np.array(c06[6]) - np.array(c06[2])
    diffs_db_exts_c06 = np.array(c06[13]) - np.array(c06[12])
    diffs_lit_db_exts_c06 = np.array(c06[15]) - np.array(c06[12])

    # G10
    diffs_db_ages_g10 = np.array(g10[4]) - np.array(g10[2])
    diffs_lit_db_ages_g10 = np.array(g10[6]) - np.array(g10[2])
    diffs_db_exts_g10 = np.array(g10[13]) - np.array(g10[12])
    diffs_lit_db_exts_g10 = np.array(g10[15]) - np.array(g10[12])

    # Calculate std, means and medians for the age differences.
    txt = ['SMC', 'LMC', 'P99', 'P00', 'C06', 'G10']
    dbs = [diffs_lit_ages_smc, diffs_lit_ages_lmc, diffs_db_ages_p99,
           diffs_db_ages_p00, diffs_db_ages_c06, diffs_db_ages_g10]
    for i, db in enumerate(dbs):
        print '{}, std = {:.3f}'.format(txt[i], np.std(db))
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
    ext_DB_data = [[p99[13], p99[14], p99[12], []],
                   [c06[13], c06[14], c06[12], []],
                   [g10[13], g10[14], g10[12], []]]
    ext_lit_DB_data = [[p99[15], p99[16], p99[12], []],
                       [c06[15], c06[16], c06[12], []],
                       [g10[15], g10[16], g10[12], []]]
    age_lit_DB_data = [[p99[6], p99[7], p99[2], p99[3]],
                       [p00[6], p00[7], p00[2], p00[3]],
                       [c06[6], c06[7], c06[2], c06[3]],
                       [g10[6], g10[7], g10[2], g10[3]]]

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

    fig = plt.figure(figsize=(13.5, 25))
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
    plt.savefig('figures/cross_match_age_ext.png', dpi=300)


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
        plt.savefig(fig_name, dpi=150)

        # Crop image.
        cmd.save_crop_img(fig_name)
