
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage.filters import gaussian_filter
from ..inp import input_params as g


def pl_bf_synth_cl(gs, x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, x_ax, y_ax,
                   synth_clst, syn_b_edges, cp_r, cp_e, shift_isoch):
    '''
    Best fit synthetic cluster obtained.
    '''
    ax = plt.subplot(gs[4:6, 8:10])
    # Set plot limits
    plt.xlim(x_min_cmd, x_max_cmd)
    plt.ylim(y_min_cmd, y_max_cmd)
    # Set axis labels
    plt.xlabel('$' + x_ax + '$', fontsize=18)
    plt.ylabel('$' + y_ax + '$', fontsize=18)
    # Set minor ticks
    ax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    # # Plot grid.
    # if g.bf_params[2] == 'dolphin':
    #     for x_ed in syn_b_edges[0]:
    #         # vertical lines
    #         ax.axvline(x_ed, linestyle=':', color='k', zorder=1)
    #     for y_ed in syn_b_edges[1]:
    #         # horizontal lines
    #         ax.axhline(y_ed, linestyle=':', color='k', zorder=1)
    #     # Add text box
    #     text = '$({};\,{})$'.format(g.bf_params[2], g.bf_params[3])
    #     ob = offsetbox.AnchoredText(text, pad=.2, loc=1, prop=dict(size=12))
    #     ob.patch.set(alpha=0.85)
    #     ax.add_artist(ob)
    # else:
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=1,
            zorder=1)
    if synth_clst.any():
        # Plot synthetic cluster.
        plt.scatter(synth_clst[0], synth_clst[2], marker='o', s=40,
                    c='#4682b4', lw=0.5, zorder=4)
        text1 = '$N_{{synth}} = {}$'.format(len(synth_clst[0]))
        text1 = text1 + '\n' + r'$M\,=\,{}\,M_{{\odot}}$'.format(int(cp_r[4]))
    else:
        text1 = '$N_{{synth}} = {}$'.format(0)
    # Add text box
    ob = offsetbox.AnchoredText(text1, pad=.2, loc=2, prop=dict(size=14))
    ob.patch.set(alpha=0.85)
    ax.add_artist(ob)
    # Plot isochrone.
    plt.plot(shift_isoch[0], shift_isoch[1], 'r', lw=1., zorder=5)

    # Add text box to the right of the synthetic cluster.
    ax_t = plt.subplot(gs[4:6, 10:12])
    ax_t.axis('off')  # Remove axis from frame.
    # Map isochrones set selection to proper name.
    iso_print = g.tracks_dict.get(g.ps_params[2])
    t1 = r'$Synthetic\;cluster\;parameters$' + '\n' + \
        r'$[Tracks:\;{}]$'.format(iso_print.replace(' ', '\;')) + '\n\n'
    t2 = r'$z\qquad\; =\, {} \pm {}$'.format(cp_r[0], cp_e[0]) + '\n'
    t3 = r'$log(age)\, =\, {} \pm {}$'.format(cp_r[1], cp_e[1]) + '\n'
    t4 = r'$E_{{(B-V)}}\;\,=\, {} \pm {}$'.format(cp_r[2], cp_e[2]) + '\n'
    t5 = r'$(m-M)_o = {} \pm {}$'.format(cp_r[3], cp_e[3]) + '\n'
    t6 = r'$M\,(M_{{\odot}})\quad\;=\,{} \pm {}$'.format(
        cp_r[4], cp_e[4]) + '\n'
    t7 = r'$b_{{frac}}\quad\,\, =\, {} \pm {}$'.format(cp_r[5], cp_e[5]) + '\n'
    text = t1 + t2 + t3 + t4 + t5 + t6 + t7
    ob = offsetbox.AnchoredText(text, pad=1, loc=6, prop=dict(size=13))
    ob.patch.set(alpha=0.85)
    ax_t.add_artist(ob)


def pl_ga_lkl(gs, l_min_max, lkl_old, model_done, new_bs_indx):
    '''
    Likelihood evolution for the GA.
    '''
    # Genetic algorithm parameters.
    N_b = g.bf_params[-1]
    n_pop, n_gen, fdif, p_cross, cr_sel, p_mut, n_el, n_ei, n_es = g.ga_params

    ax = plt.subplot(gs[6:8, 0:4])
    plt.xlim(-0.5, len(lkl_old[0]) + int(0.01 * len(lkl_old[0])))
    plt.ylim(l_min_max[0], l_min_max[1])
    # Set minor ticks
    ax.minorticks_on()
    ax.tick_params(axis='y', which='major', labelsize=9)
    ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.6)
    plt.xlabel('Generation', fontsize=12)
    plt.ylabel('Likelihood', fontsize=12)
    # Add text box.
    text1 = '$N_{total} = %.2e\,;\,N_{btst} = %d$' '\n' % \
        (len(model_done[0]), N_b)
    text2 = '$n_{gen}=%d\,;\,n_{pop}=%d$' '\n' % (n_gen, n_pop)
    text3 = '$f_{dif}=%0.2f\,;\,cr_{sel}=%s$' '\n' % (fdif, cr_sel)
    text4 = '$p_{cross}=%0.2f\,;\,p_{mut}=%0.2f$' '\n' % (p_cross, p_mut)
    text5 = '$n_{el}=%d\,;\,n_{ei}=%d\,;\,n_{es}=%d$' % (n_el, n_ei, n_es)
    text = text1 + text2 + text3 + text4 + text5
    ob = offsetbox.AnchoredText(text, loc=1, prop=dict(size=12))
    ob.patch.set(boxstyle='square,pad=0.', alpha=0.85)
    ax.add_artist(ob)

    # Plot likelihood minimum and mean lines.
    ax.plot(range(len(lkl_old[0])), lkl_old[0], lw=1., c='black',
            label='$L_{{min}}={:.1f}$'.format(min(lkl_old[0])))
    ax.plot(range(len(lkl_old[0])), lkl_old[1], lw=1., c='blue',
            label='$L_{mean}$')
    # Plot line marking a new best solution found.
    for lin in new_bs_indx:
        lw_lin = 2. if lin > 0.05 * len(lkl_old[0]) else 0.5
        plt.axvline(x=lin, linestyle='--', lw=lw_lin, color='green')
    # Legend.
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels, loc='upper left', numpoints=1,
                    fontsize=12)
    leg.get_frame().set_alpha(0.7)


def pl_2_param_dens(gs, _2_params, min_max_p, cp_r, cp_e, model_done):
    '''
    Parameter vs parameters solutions density map.
    '''
    # Define parameters for upper and lower plots.
    if _2_params == 'age-metal':
        ax, cp, d_map, mx, my = plt.subplot(gs[6:8, 4:6]), 'r', 'Blues', 0, 1
        x_label, y_label = '$z$', '$log(age)$'
    elif _2_params == 'dist-ext':
        ax, cp, d_map, mx, my = plt.subplot(gs[6:8, 6:8]), 'b', 'Reds', 2, 3
        x_label, y_label = '$E_{(B-V)}$', '$(m-M)_o$'
    elif _2_params == 'metal-dist':
        ax, cp, d_map, mx, my = plt.subplot(gs[6:8, 8:10]), 'r', 'Blues', 0, 3
        x_label, y_label = '$z$', '$(m-M)_o$'
    elif _2_params == 'mass-binar':
        ax, cp, d_map, mx, my = plt.subplot(gs[6:8, 10:12]), 'b', 'Reds', 4, 5
        x_label, y_label = '$M\,(M_{{\odot}})$', '$b_{frac}$'

    # Parameter values and errors.
    xp, e_xp = cp_r[0][mx], cp_e[mx]
    yp, e_yp = cp_r[0][my], cp_e[my]
    # Axis limits.
    xp_min, xp_max = min_max_p[mx]
    yp_min, yp_max = min_max_p[my]
    # Special axis ticks for metallicity.
    if _2_params in {'age-metal', 'metal-dist'}:
        z_xmin, z_step = min_max_p[-1]
        ax.xaxis.set_ticks(np.arange(z_xmin, xp_max, z_step))

    plt.xlim(xp_min, xp_max)
    plt.ylim(yp_min, yp_max)
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.minorticks_on()
    # Plot best fit point.
    plt.scatter(xp, yp, marker='o', c=cp, s=30)
    # Check if errors in both dimensions are defined.
    if all([i > 0. for i in [e_xp, e_yp]]):
        # Plot ellipse error.
        plt.gca()
        ellipse = Ellipse(xy=(xp, yp), width=2 * e_xp, height=2 * e_yp,
                          edgecolor=cp, fc='None', lw=1.)
        ax.add_patch(ellipse)
    # Else plot an error bar in the corresponding dimension.
    elif e_xp < 0. and e_yp > 0.:
        plt.errorbar(xp, yp, yerr=e_yp, color=cp)
    elif e_yp < 0. and e_xp > 0.:
        plt.errorbar(xp, yp, xerr=e_xp, color=cp)
    # Plot density map.
    hist, xedges, yedges = np.histogram2d(zip(*model_done[0])[mx],
                                          zip(*model_done[0])[my], bins=100)
    # H_g is the 2D histogram with a Gaussian filter applied
    h_g = gaussian_filter(hist, 2, mode='constant')
    plt.imshow(h_g.transpose(), origin='lower',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
               cmap=plt.get_cmap(d_map), aspect='auto')


def pl_lkl_scatt(gs, ld_p, min_max_p, cp_r, cp_e, model_done):
    '''
    Parameter likelihood scatter plot.
    '''
    # Define parameters for upper and lower plots.
    if ld_p == '$z$':
        ax, cp, ci, zlab = plt.subplot(gs[8:10, 0:2]), 0, 1, '$log(age)$'
        plt.ylabel('Likelihood', fontsize=12)
    elif ld_p == '$log(age)$':
        ax, cp, ci, zlab = plt.subplot(gs[8:10, 2:4]), 1, 2, '$E_{{(B-V)}}$'
    elif ld_p == '$E_{{(B-V)}}$':
        ax, cp, ci, zlab = plt.subplot(gs[8:10, 4:6]), 2, 3, '$(m-M)_o$'
    elif ld_p == '$(m-M)_o$':
        ax, cp, ci, zlab = plt.subplot(gs[8:10, 6:8]), 3, 1, '$log(age)$'
    elif ld_p == '$M\,(M_{{\odot}})$':
        ax, cp, ci, zlab = plt.subplot(gs[8:10, 8:10]), 4, 5, '$b_{{frac}}$'
    elif ld_p == '$b_{{frac}}$':
        ax, cp, ci, zlab = plt.subplot(gs[8:10, 10:12]), 5, 4,\
            '$M\,(M_{{\odot}})$'

    # Parameter values and errors.
    xp, e_xp = cp_r[0][cp], cp_e[cp]
    # Set x axis limit.
    xp_min, xp_max = min_max_p[cp]
    # Special axis ticks for metallicity.
    if ld_p == '$z$':
        z_xmin, z_step = min_max_p[-1]
        ax.xaxis.set_ticks(np.arange(z_xmin, xp_max, z_step))
    plt.xlim(xp_min, xp_max)
    # Set minor ticks
    ax.minorticks_on()
    ax.tick_params(axis='y', which='major', labelsize=9)
    plt.xlabel(ld_p, fontsize=16)
    # Add textbox.
    text = (ld_p + '$ = {} \pm {}$').format(xp, e_xp)
    ob = offsetbox.AnchoredText(text, pad=0.1, loc=2, prop=dict(size=12))
    ob.patch.set(alpha=0.8)
    ax.add_artist(ob)
    plt.axvline(x=xp, linestyle='--', color='red', zorder=2)
    # Plot scatter points over likelihood density map.
    cm = plt.cm.get_cmap('viridis')
    col_arr = [float(_) for _ in zip(*model_done[0])[ci]]
    SC = plt.scatter(zip(*model_done[0])[cp], model_done[1], marker='o',
                     c=col_arr, s=25, edgecolors='k',
                     lw=0.2, edgecolor='w', cmap=cm, zorder=3)
    if e_xp > 0.:
        # Plot error bars only if errors where assigned.
        plt.axvline(x=xp + e_xp, linestyle='--', color='blue')
        plt.axvline(x=xp - e_xp, linestyle='--', color='blue')
    # Set y axis limit.
    min_lik = min(model_done[1])
    if min_lik > 0:
        min_y, max_y = min_lik - min_lik*0.1, 2.5*min_lik
    else:
        min_y, max_y = min_lik + min_lik*0.1, -2.5*min_lik
    plt.ylim(min_y, max_y)
    # Position colorbar.
    the_divider = make_axes_locatable(ax)
    color_axis = the_divider.append_axes("right", size="2%", pad=0.1)
    # Colorbar.
    cbar = plt.colorbar(SC, cax=color_axis)
    cbar.set_label(zlab, fontsize=12, labelpad=4)
    cbar.ax.tick_params(labelsize=8)


def plot(N, *args):
    '''
    Handle each plot separately.
    '''

    plt_map = {
        0: [pl_bf_synth_cl, 'synthetic cluster'],
        1: [pl_ga_lkl, 'GA likelihood evolution'],
        2: [pl_2_param_dens, 'age vs metallicity density map'],
        3: [pl_2_param_dens, 'distance vs extinction density map'],
        4: [pl_2_param_dens, 'z vs distance density map'],
        5: [pl_2_param_dens, 'mass vs binarity density map'],
        6: [pl_lkl_scatt, 'z likelihood scatter'],
        7: [pl_lkl_scatt, 'age likelihood scatter'],
        8: [pl_lkl_scatt, 'extinction likelihood scatter'],
        9: [pl_lkl_scatt, 'distance likelihood scatter'],
        10: [pl_lkl_scatt, 'mass likelihood scatter'],
        11: [pl_lkl_scatt, 'binarity likelihood scatter']
    }

    fxn = plt_map.get(N, None)[0]
    if fxn is None:
        raise ValueError("  ERROR: there is no plot {}.".format(N))

    try:
        fxn(*args)
    except:
        import traceback
        print traceback.format_exc()
        print("  WARNING: error when plotting {}.".format(plt_map.get(N)[1]))
