

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.offsetbox as offsetbox
from matplotlib.ticker import MultipleLocator
import os
import re


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def read_met_file(met_f, age_values, tr):
    '''
    Read a given metallicity file and return the isochrones for the ages
    within the age range.
    '''

    # Read line start format and columns indexes for the selected set of
    # Girardi isochrones.
    line_start = "#\tIsochrone\tZ = " if tr == 0 else "#\tIsochrone  Z = "
    # Define reg expression to isolate the age of an isochrone.
    age_format = r"Age = \t(.+?) yr"
    idxs = [4, 3] if tr == 0 else [5, 4]

    # Initialize list that will hold all the isochrones for this
    # metallicity value.
    metal_isoch = []

    # Open the metallicity file.
    with open(met_f, mode="r") as f_iso:

        # Define empty lists.
        isoch_teff, isoch_llo = [], []

        # Initial value for age to avoid 'not defined' error.
        age = -99.

        # Iterate through each line in the file.
        for line in f_iso:

            # Identify beginning of a defined isochrone.
            if line.startswith(line_start):

                # Save stored values if these exist.
                # Skip first age for which the lists will be empty.
                if isoch_teff:
                    # Store color, magnitudes and masses for this
                    # isochrone.
                    metal_isoch.append([isoch_teff, isoch_llo])
                    # Reset lists.
                    isoch_teff, isoch_llo = [], []

                # Read age value for this isochrone.
                age0 = re.findall(age_format, line)  # Find age in line.
                age = np.around(np.log10(float(age0[0])), 2)

            # If age value falls inside the given range, store the
            # isochrone's data.
            if age in age_values:

                # Save mag, color and mass values for each isochrone.
                if not line.startswith("#"):
                    reader = line.split()
                    # log(Teff).
                    isoch_teff.append(float(reader[idxs[0]]))
                    # log(L/Lo)
                    isoch_llo.append(float(reader[idxs[1]]))

        # Save the last isochrone when EOF is reached.
        else:
            # If list is not empty.
            if isoch_teff:
                metal_isoch.append([isoch_teff, isoch_llo])

    return metal_isoch


def get_data():
    """
    """
    # Path to data files.
    mar, par = 'mar2008_ubvrijhk/', 'parsec11_ubvrijhk/'
    age_values = [7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]

    # Read data file
    mar_data, par_data = [], []
    for tr, tracks in enumerate([mar, par]):
        for m_f in os.listdir(tracks):
            met_f = tracks + m_f
            metal_isoch = read_met_file(met_f, age_values, tr)
            if tr == 0:
                mar_data.append(metal_isoch)
            else:
                par_data.append(metal_isoch)

    return mar_data, par_data


def make_plot(mar_data, par_data):
    """
    """

    xmin, xmax, ymin, ymax = 3.3, 4.8, -1., 4.8

    # Generate plot.
    fig = plt.figure(figsize=(20, 25))
    gs = gridspec.GridSpec(2, 2)

    met = ['0.001', '0.004', '0.015', '0.03']
    xy_font_s = 21
    for i in range(4):
        ax = plt.subplot(gs[i], aspect='auto')
        plt.xticks(fontsize=xy_font_s - 6)
        plt.yticks(fontsize=xy_font_s - 6)
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.xlabel(r'$\log(T_{eff})$', fontsize=xy_font_s)
        plt.ylabel(r'$\log(L/L_{\odot})$', fontsize=xy_font_s)
        ax.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
                zorder=1)
        ax.minorticks_on()
        # Text box.
        ob = offsetbox.AnchoredText(r'z = {}'.format(met[i]), loc=9,
                                    prop=dict(size=xy_font_s - 4))
        ob.patch.set(alpha=0.85)
        ax.add_artist(ob)
        col = ['r', 'b', 'g', 'k', 'c', 'm', 'y']
        l = ['7.0', '7.5', '8.0', '8.5', '9.0', '9.5', '10.0']
        for j, track in enumerate(mar_data[i]):
            plt.plot(track[0], track[1], c=col[j], label=l[j])
        # Legend.
        leg = plt.legend(loc='lower left', markerscale=2., scatterpoints=1,
                         fontsize=xy_font_s - 4)
        # Set the alpha value of the legend.
        leg.get_frame().set_alpha(0.85)
        ax.set_aspect('auto')
    # Position colorbar.
        for j, track in enumerate(par_data[i]):
            plt.plot(track[0], track[1], c=col[j], ls='--')
        plt.gca().invert_xaxis()

    # Output png file.
    fig.tight_layout()
    plt.savefig('mar_vs_par_isochs.png', dpi=300, bbox_inches='tight')


def main():
    '''
    Call each function.
    '''
    mar_data, par_data = get_data()

    make_plot(mar_data, par_data)

    print '\nEnd.'


if __name__ == "__main__":
    main()
