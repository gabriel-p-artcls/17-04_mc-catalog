
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.angle_helper as angle_helper
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D
from mpl_toolkits.axisartist import SubplotHost
from mpl_toolkits.axisartist import GridHelperCurveLinear


def curvelinear_test2(fig, rect=111):
    """
    Polar projection, but in a rectangular box.
    """

    # see demo_curvelinear_grid.py for details
    tr = Affine2D().translate(0, 90) + Affine2D().scale(np.pi / 180., 1.) + \
        PolarAxes.PolarTransform()

    extreme_finder = angle_helper.ExtremeFinderCycle(10, 60,
                                                     lon_cycle=360,
                                                     lat_cycle=None,
                                                     lon_minmax=None,
                                                     lat_minmax=(-90, np.inf),
                                                     )
    # Changes theta gridline count
    grid_locator1 = angle_helper.LocatorHMS(12)
    grid_locator2 = angle_helper.LocatorDMS(6)
    tick_formatter1 = angle_helper.FormatterHMS()
    tick_formatter2 = angle_helper.FormatterDMS()

    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        grid_locator2=grid_locator2,
                                        tick_formatter1=tick_formatter1,
                                        tick_formatter2=tick_formatter2
                                        )

    ax1 = SubplotHost(fig, rect, grid_helper=grid_helper)

    # make ticklabels of right and top axis visible.
    ax1.axis["right"].major_ticklabels.set_visible(True)
    ax1.axis["top"].major_ticklabels.set_visible(True)
    ax1.axis["bottom"].major_ticklabels.set_visible(True)
    # let right and bottom axis show ticklabels for 1st coordinate (angle)
    ax1.axis["right"].get_helper().nth_coord_ticks = 0
    ax1.axis["bottom"].get_helper().nth_coord_ticks = 0

    #
    fig.add_subplot(ax1)

    grid_helper = ax1.get_grid_helper()

    # You may or may not need these - they set the view window explicitly
    # rather than using the default as determined by matplotlib with extreme
    # finder.
    ax1.set_aspect(1.)
    ax1.set_xlim(-4, 25)  # moves the origin left-right in ax1
    ax1.set_ylim(-2.5, 30)  # moves the origin up-down

    ax1.set_ylabel('$DEC\,(^{\circ})$')
    ax1.set_xlabel('$RA\,(h)$')
    ax1.grid(True)
    # ax1.grid(linestyle='--', which='x') # either keyword applies to both
    # ax1.grid(linestyle=':', which='y')  # sets of gridlines

    return ax1, tr


def ra_dec_plots(pl_params):
    '''
    Generate RA vs DEC plots.
    '''

    fig, gs, ra, dec, data_arr, rad_pc, z_lab = pl_params

    # tr.transform_point((x, 0)) is always (0,0)
    ax1, tr = curvelinear_test2(fig, gs)
    # Get transformed data.
    ra_dec_tr = tr.transform(zip(ra, dec))

    # Define colormap.
    cm = plt.cm.get_cmap('RdYlBu_r')

    # Size relative to the clusters actual size in pc.
    siz = np.asarray(rad_pc) * 4.
    SC = ax1.scatter(ra_dec_tr[:, 0], ra_dec_tr[:, 1], marker='o', s=siz,
                     c=data_arr, cmap=cm, lw=0.1, zorder=9)

    # Plot clouds center. Central coords stored in degrees.
    c_SMC = [13.1875, -72.82861111]
    c_LMC = [80.2375, -69.47805556]
    # Get transformed data.
    clouds_cent = tr.transform([c_SMC, c_LMC])
    plt.scatter(clouds_cent[:, 0], clouds_cent[:, 1], marker='v', s=50,
                c='g', edgecolor='w', lw=0.5, zorder=10)
    # Colorbar
    cbar = plt.colorbar(SC, shrink=1., pad=0.05)
    cbar.ax.tick_params(labelsize=8)
    # cbar.set_clim(0., 0.4)
    cbar.set_label(z_lab, fontsize=12)
