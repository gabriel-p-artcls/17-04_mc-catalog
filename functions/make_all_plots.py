
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

from functions.ra_dec_map import ra_dec_plots
from functions.kde_2d import kde_map


def as_vs_lit_plots(pl_params):
    '''
    Generate ASteCA vs literature values plots.
    '''
    gs, i, xmin, xmax, x_lab, y_lab, z_lab, xarr, xsigma, yarr, ysigma, \
    zarr = pl_params

    xy_font_s = 16
    cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i], aspect='equal')
    #ax.set_aspect('auto')
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    if x_lab == '$(m-M)_{0;\,asteca}$':
        # Introduce random scatter.
        # 10% of axis ranges.
        xax_ext = (xmax - xmin) * 0.05
        rs_x = np.random.uniform(0., xax_ext, len(xarr))
        rs_y = np.random.uniform(0., xax_ext, len(xarr))
        xarr = xarr + rs_x
        yarr = yarr + rs_y
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
        zorder=1)
    ax.minorticks_on()
    plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')  # 1:1 line
    # Plot all clusters in dictionary.
    SC = plt.scatter(xarr, yarr, marker='o', c=zarr, s=70, lw=0.25, cmap=cm,
        zorder=3)
    # Only plot y error bar if it has a value assigned in the literarute.
    for j, xy in enumerate(zip(*[xarr, yarr])):
        if ysigma:  # Check if list is not empty (radii list)
            if ysigma[j] > -99.:
                plt.errorbar(xy[0], xy[1], xerr=xsigma[j], yerr=ysigma[j],
                    ls='none', color='k', elinewidth=0.8, zorder=1)
            else:
                plt.errorbar(xy[0], xy[1], xerr=xsigma[j], ls='none',
                    color='k', elinewidth=0.8, zorder=1)
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="5%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    zpad = 10 if z_lab == '$E_{(B-V)}$' else 5
    cbar.set_label(z_lab, fontsize=xy_font_s - 2, labelpad=zpad)


def make_as_vs_lit_plot(i, galax, k, in_params):
    '''
    Prepare parameters and call function to generate ASteca vs literature
    SMC and LMC plots.
    '''

    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, rarr = \
    [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
    'esigma', 'darr', 'dsigma', 'rarr']]

    # Generate ASteca vs literature plots.
    fig = plt.figure(figsize=(17, 26))  # create the top-level container
    #gs = gridspec.GridSpec(2, 4, width_ratios=[1, 0.35, 1, 0.35])
    gs = gridspec.GridSpec(4, 2)

    if galax == 'SMC':
        dm_min, dm_max = 18.6, 19.4
    else:
        dm_min, dm_max = 18.1, 18.9

    as_lit_pl_lst = [
    [gs, 0, -1.8, 0.45, '$[Fe/H]_{asteca}$', '$[Fe/H]_{lit}$', '$E_{(B-V)}$',
        zarr[k][0], zsigma[k][0], zarr[k][1], zsigma[k][1], earr[k][0]],
    [gs, 1, 5.9, 10.6, '$log(age)_{asteca}$', '$log(age)_{lit}$',
        '$E_{(B-V)}$', aarr[k][0], asigma[k][0], aarr[k][1], asigma[k][1],
        earr[k][0]],
    [gs, 2, -0.05, 0.3, '$E(B-V)_{asteca}$', '$E(B-V)_{lit}$',
        '$log(age)_{asteca}$', earr[k][0], esigma[k][0], earr[k][1],
        esigma[k][1], aarr[k][0]],
    [gs, 3, dm_min, dm_max, '$(m-M)_{0;\,asteca}$', '$(m-M)_{0;\,lit}$',
        '$log(age)_{asteca}$', darr[k][0], dsigma[k][0], darr[k][1],
        dsigma[k][1], aarr[k][0]],
    [gs, 4, 0., 500., '$rad_{asteca}$', '$rad_{lit}$',
        '$log(age)_{asteca}$', rarr[k][0], [], rarr[k][1], [], aarr[k][0]]
        ]
    #
    for pl_params in as_lit_pl_lst:
        as_vs_lit_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig(galax + '_as_vs_lit_' + i + '.png', dpi=300)


def kde_plots(pl_params):
    '''
    Generate KDE plots.
    '''
    gs, i, x_lab, y_lab, xarr, xsigma, yarr, ysigma = pl_params
    # Generate maps.
    z, ext = kde_map(np.array(xarr), np.array(xsigma), np.array(yarr),
        np.array(ysigma))

    # 10% of axis ranges.
    xax_ext = (ext[1] - ext[0]) * 0.1
    yax_ext = (ext[3] - ext[2]) * 0.1
    # Extend plotting range.
    ext2 = [ext[0] - xax_ext, ext[1] + xax_ext, ext[2] - yax_ext,
        ext[3] + yax_ext]

    # Make plot.
    ax = plt.subplot(gs[i])
    xy_font_s = 21
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)

    cm = plt.cm.get_cmap('RdYlBu_r')
    #cm = plt.cm.gist_earth_r
    ax.imshow(z, cmap=cm, extent=ext2)
    ax.set_aspect('auto')
    #plt.errorbar(xarr, yarr, xerr=xsigma, yerr=ysigma, fmt=None, color='k')
    # Random scatter.
    rs_x = np.random.uniform(0., 0.1 * xax_ext, len(xarr))
    rs_y = np.random.uniform(0., 0.1 * yax_ext, len(xarr))
    plt.scatter(xarr + rs_x, yarr + rs_y, marker='*', color='#6b6868', s=40,
        lw=0.5, facecolors='none')
    ax.set_xlim([ext[0] - xax_ext, ext[1] + xax_ext])
    ax.set_ylim([ext[2] - yax_ext, ext[3] + yax_ext])


def make_kde_plots(i, galax, k, in_params):
    '''
    Prepare parameters and call function to generate SMC and LMC KDE plots.
    '''
    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, marr, msigma = \
    [in_params[_] for _ in ['zarr', 'zsigma', 'aarr', 'asigma', 'earr',
    'esigma', 'darr', 'dsigma', 'marr', 'msigma']]

    fig = plt.figure(figsize=(14, 25))  # create the top-level container
    gs = gridspec.GridSpec(4, 2)       # create a GridSpec object

    kde_pl_lst = [
    [gs, 0, '$log(age)$', '$[Fe/H]$', aarr[k][0], asigma[k][0], zarr[k][0],
        zsigma[k][0]],
    [gs, 1, '$log(age)$', '$M\,(M_{\odot})$', aarr[k][0], asigma[k][0],
        marr[k][0], msigma[k][0]],
    [gs, 2, '$(m-M)_0$', '$E_{(B-V)}$', darr[k][0], dsigma[k][0],
        earr[k][0], esigma[k][0]],
    [gs, 3, '$M\,(M_{\odot})$', '$[Fe/H]$', marr[k][0], msigma[k][0],
        zarr[k][0], zsigma[k][0]]
    #[gs, 4, '$log(age)$', '$M\,(M_{\odot})$', aarr[k][0], asigma[k][0],
        #marr[k][0], msigma[k][0]],
    #[gs, 5, '$log(age)$', '$M\,(M_{\odot})$', aarr[k][0], asigma[k][0],
        #marr[k][0], msigma[k][0]]
    ]
    #
    for pl_params in kde_pl_lst:
        kde_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig(galax + '_kde_maps_' + i + '.png', dpi=300)


def make_ra_dec_plots(i, in_params):
    '''
    Prepare parameters and call function to generate RA vs DEC positional
    plots for the SMC and LMC.
    '''

    ra, dec, zarr, aarr, earr, darr, marr, rad_pc = [in_params[_] for _ in
    ['ra', 'dec', 'zarr', 'aarr', 'earr', 'darr', 'marr', 'rad_pc']]

    # Put both SMC and LMC clusters into a single list.
    ra = ra[0] + ra[1]
    dec = dec[0] + dec[1]
    zarr = zarr[0][0] + zarr[1][0]
    aarr = aarr[0][0] + aarr[1][0]
    earr = earr[0][0] + earr[1][0]
    darr = darr[0][0] + darr[1][0]
    marr = marr[0][0] + marr[1][0]
    rad_pc = rad_pc[0] + rad_pc[1]

    fig = plt.figure(figsize=(20, 20))
    fig.clf()

    ra_dec_pl_lst = [
    [fig, 321, ra, dec, zarr, '$[Fe/H]$'],
    [fig, 322, ra, dec, aarr, '$log(age)$'],
    [fig, 323, ra, dec, earr, '$E_{(B-V)}$'],
    [fig, 324, ra, dec, darr, '$(m-M)_0$'],
    [fig, 325, ra, dec, marr, '$M\,(M_{\odot})$'],
    [fig, 326, ra, dec, rad_pc, '$r_{clust}\,[pc]$']
    ]

    for pl_params in ra_dec_pl_lst:
        ra_dec_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('RA_DEC_maps_ ' + i + '.png', dpi=300)


def make_lit_ext_plot(in_params):
    '''
    Prepare parameters and call function to generate RA vs DEC positional
    plots for the SMC and LMC.
    '''

    earr, esigma, ext_sf, ext_mcev = [in_params[_] for _ in ['earr', 'esigma',
    'ext_sf', 'ext_mcev']]

    # Define values to pass.
    xmin, xmax = -0.02, 0.4
    x_lab, y_lab, z_lab = '$E(B-V)_{asteca}$', '$E(B-V)_{MCEV}$', \
    '$E(B-V)_{SF}$'

    fig = plt.figure(figsize=(16, 25))  # create the top-level container
    gs = gridspec.GridSpec(4, 2)       # create a GridSpec object

    ext_pl_lst = [
        # SMC
        [gs, 0, xmin, xmax, x_lab, y_lab, z_lab, earr[0][0], esigma[0][0],
            ext_mcev[0][0], ext_mcev[0][1], ext_sf[0][0]],
        # LMC
        [gs, 1, xmin, xmax, x_lab, y_lab, z_lab, earr[1][0], esigma[1][0],
            ext_mcev[1][0], ext_mcev[1][1], ext_sf[1][0]]
    ]

    for pl_params in ext_pl_lst:
        as_vs_lit_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('ext_lit_plot.png', dpi=300)


def int_cols_plots(pl_params):
    '''
    Generate ASteCA integrated colors plots.
    '''
    gs, i, xmin, xmax, x_lab, y_lab, z_lab, xarr, xsigma, yarr, zarr, rad = \
    pl_params
    siz = np.asarray(rad) * 5

    xy_font_s = 16
    cm = plt.cm.get_cmap('RdYlBu_r')

    ax = plt.subplot(gs[i])
    #ax.set_aspect('auto')
    plt.xlim(xmin, xmax)
    #plt.ylim(xmin, xmax)
    plt.xlabel(x_lab, fontsize=xy_font_s)
    plt.ylabel(y_lab, fontsize=xy_font_s)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
        zorder=1)
    ax.minorticks_on()
    # Plot all clusters in dictionary.
    SC = plt.scatter(xarr, yarr, marker='o', c=zarr, s=siz, lw=0.25, cmap=cm,
        zorder=3)
    # Only plot y error bar if it has a value assigned in the literarute.
    plt.errorbar(xarr, yarr, xerr=xsigma, ls='none', color='k',
        elinewidth=0.8, zorder=1)
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    zpad = 10 if z_lab == '$E_{(B-V)}$' else 5
    cbar.set_label(z_lab, fontsize=xy_font_s - 2, labelpad=zpad)


def make_int_cols_plot(in_params):
    '''
    Prepare parameters and call function to generate integrated color vs Age
    (colored by mass) plots for the SMC and LMC.
    '''

    aarr, asigma, marr, int_colors, rad_pc = [in_params[_] for _ in ['aarr',
    'asigma', 'marr', 'int_colors', 'rad_pc']]

    # Define values to pass.
    xmin, xmax = 6.5, 10.4
    x_lab, y_lab, z_lab = '$log(age/yr)_{asteca}$', \
    '$(C-T_{1})_{0;\,asteca}$', '$M\,(M_{\odot})$'

    fig = plt.figure(figsize=(16, 25))  # create the top-level container
    gs = gridspec.GridSpec(4, 1)       # create a GridSpec object

    ext_pl_lst = [
        # SMC
        [gs, 0, xmin, xmax, x_lab, y_lab, z_lab, aarr[0][0], asigma[0][0],
            int_colors[0], marr[0][0], rad_pc[0]],
        ## LMC
        #[gs, 1, xmin, xmax, x_lab, y_lab, z_lab, aarr[1][0], asigma[1][0],
            #int_colors[1], marr[1][0], rad_pc[1]]
    ]

    for pl_params in ext_pl_lst:
        int_cols_plots(pl_params)

    # Output png file.
    fig.tight_layout()
    plt.savefig('int_colors_plot.png', dpi=300)