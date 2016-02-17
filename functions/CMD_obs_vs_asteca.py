
import Image
import ImageOps
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from scipy import stats
from read_photom_files import get_data as gd
import glob


def kde_limits(phot_x, phot_y):
    '''
    Return photometric diagram limits taken from a 2D KDE.
    '''
    xmin, xmax = min(phot_x), max(phot_x)
    ymin, ymax = min(phot_y), max(phot_y)
    # Stack photometric data.
    values = np.vstack([phot_x, phot_y])
    # Obtain Gaussian KDE.
    kernel = stats.gaussian_kde(values)
    # Grid density (number of points).
    gd = 25
    gd_c = complex(0, gd)
    # Define x,y grid.
    x, y = np.mgrid[xmin:xmax:gd_c, ymin:ymax:gd_c]
    positions = np.vstack([x.ravel(), y.ravel()])
    # Evaluate kernel in grid positions.
    k_pos = kernel(positions)

    # Generate 30 contour lines.
    cs = plt.contour(x, y, np.reshape(k_pos, x.shape), 30)
    # Extract (x,y) points delimitating each line.
    x_v, y_v = np.asarray([]), np.asarray([])
    # Only use the outer curve.
    col = cs.collections[0]
    # If more than one region is defined by this curve (ie: the main sequence
    # region plus a RC region or some other detached region), obtain x,y from
    # all of them.
    for lin in col.get_paths():
        x_v = np.append(x_v, lin.vertices[:, 0])
        y_v = np.append(y_v, lin.vertices[:, 1])

    return x_v, y_v


def diag_limits(r_path, cl):
    '''
    Define plot limits for *all* photometric diagrams.
    '''
    y_axis = 0
    path_no_ext = r_path + 'asteca-project/asteca/input/dont_read/'\
        + 'MC_all/' + cl + '.*'
    try:
        data_file = glob.glob(path_no_ext)[0]
    except IndexError:
        print ("The folder: {}asteca-project/asteca/ does not"
               " exist".format(r_path))
        raise SystemExit(0)
    phot_data = gd(data_file)
    phot_x, phot_y = phot_data[5], phot_data[3]

    x_v, y_v = kde_limits(phot_x, phot_y)

    # Define diagram limits.
    x_min_cmd, x_max_cmd = min(x_v) - 1.25, max(x_v) + 1.25
    y_min_cmd = max(y_v) + 1.25
    # If photometric axis y is a magnitude, make sure the brightest star
    # is always plotted.
    if y_axis == 0:
        y_max_cmd = min(phot_y) - 1.
    else:
        y_max_cmd = min(y_v) - 1.

    return x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol, except those that start with
    '#>>>' which indicates the run where the membership data is stored.
    '''
    for line in f:
        if line.strip().startswith('#>>>'):
            yield line
        elif not line.strip().startswith('#'):
            yield line


def get_DB_age_ext(r_path, cl, db):
    '''
    Read age and extinction values (and Galaxy) for the 'cl' cluster matched
    in the 'db' database.
    '''
    f_path = r_path + 'mc-catalog/ages_mass_lit/matched_clusters.dat'
    # Read data file
    with open(f_path) as f:
        for line in skip_comments(f):
            l = line.split()
            if l[0] == db and l[2] == cl:
                a, e, gal = line.split()[3], line.split()[15], line.split()[1]
                return a, e, gal


def get_cl_run(cl):
    '''
    Read the "run" that should be used for this cluster, out of the several
    ones made.
    '''
    # Path to data file.
    out_file = 'asteca_output_final.dat'

    run = ''
    # Read data file
    with open(out_file) as f:

        for line in skip_comments(f):
            l = line.split(' ')
            if l[0] == cl:
                return run
            elif l[0] == '#>>>':
                run = l[1]

    return run


def get_input_folder(r_path, cl, run):
    '''
    Find which 'input_XX' folder for this cluster in this "run" contains
    its membership file.
    '''
    name = cl + '.png'
    path = r_path + 'mc-catalog/runs/' + run + '_run/'
    for root, dirs, files in os.walk(path):
        if name in files:
            full_path = os.path.join(root, name)

    # Extract name of input folder.
    inpt = full_path.split('/')[-2]

    return inpt


def get_memb_data(r_path, run, inpt, cl):
    '''
    Read the cluster membership file. Divide into stars used in the best fit
    process and stars that were not.
    '''
    # Path where the members file is stored.
    cl_path = r_path + 'mc-catalog/runs/' + run + \
        '_run/output/' + inpt + '/' + cl + '_memb.dat'
    # Read data file
    with open(cl_path) as f:
        cl_reg_fit = [[], [], []]
        cl_reg_no_fit = [[], [], []]

        for line in skip_comments(f):
            if line.split()[8] == '1':
                # Store stars used in fit.
                cl_reg_fit[0].append(float(line.split()[5]))
                cl_reg_fit[1].append(float(line.split()[3]))
                cl_reg_fit[2].append(float(line.split()[7]))
            else:
                cl_reg_no_fit[0].append(float(line.split()[5]))
                cl_reg_no_fit[1].append(float(line.split()[3]))
                cl_reg_no_fit[2].append(float(line.split()[7]))

    # Create new list with inverted values so higher prob stars are on top.
    cl_reg_fit = [i[::-1] for i in cl_reg_fit]

    return cl_reg_fit, cl_reg_no_fit


def move_isoch(isochrone, e, d):
    '''
    Receives an isochrone of a given age and metallicity and modifies
    its color and magnitude values according to given values for the extinction
    E(B-V) (e) and distance modulus (d).
    '''
    iso_moved = [[], []]

    # For Washington system.
    #
    # E(C-T1) = 1.97*E(B-V) = (C-T1) - (C-T)o
    # M_T1 = T1 + 0.58*E(B-V) - (m-M)o - 3.2*E(B-V)
    #
    # (C-T1) = (C-T1)o + 1.97*E(B-V)
    # T1 = M_T1 - 0.58*E(B-V) + (m-M)o + 3.2*E(B-V)
    #
    e, d = float(e), float(d)
    V_Mv = d + 3.2 * e
    iso_moved = [np.array(isochrone[0]) + 1.97 * e,
                 np.array(isochrone[1]) - 0.58 * e + V_Mv]

    return iso_moved


def get_isoch(r_path, DB_asteca, z, a, e, d):
    '''
    Read a given metallicity file and return the isochrones for the age passed,
    moved according to the extinction and distance modulus values.
    '''
    if DB_asteca == 'DB':
        # Use Marigo isochrones.
        met_f = r_path + 'mc-catalog/functions/' + str(z) + '.dat'
        line_start, imass_idx = "#\tIsochrone\tZ = ", 1
        # T1, C
        mag1_idx, mag2_idx = 9, 7
    else:
        # Use PARSEC isochrones.
        line_start, imass_idx = "#\tIsochrone  Z = ", 2
        mag1_idx, mag2_idx = 10, 8
        met_f = r_path + 'asteca-project/asteca/isochrones/' + \
            'parsec11_washington/' + str(z) + '.dat'
    age_format = r"Age = \t(.+?) yr"
    cmd_select = 4

    # Initialize list that will hold all the isochrones for this
    # metallicity value.
    metal_isoch = []

    # Open the metallicity file.
    with open(met_f, mode="r") as f_iso:

        # Define empty lists.
        isoch_col, isoch_mag, isoch_mas = [], [], []

        # Initial value for age to avoid 'not defined' error.
        age = -99.

        # Iterate through each line in the file.
        for line in f_iso:

            # Identify beginning of a defined isochrone.
            if line.startswith(line_start):

                # Save stored values if these exist.
                # Skip first age for which the lists will be empty.
                if isoch_col:
                    # Store color, magnitudes and masses for this
                    # isochrone.
                    metal_isoch.append([isoch_col, isoch_mag,
                                        isoch_mas])
                    # Reset lists.
                    isoch_col, isoch_mag, isoch_mas = [], [], []

                # Read age value for this isochrone.
                age0 = re.findall(age_format, line)  # Find age in line.
                age = np.around(np.log10(float(age0[0])), 2)

            # If age value falls inside the given range, store the
            # isochrone's data.
            if float(age) == float(a):

                # Save mag, color and mass values for each isochrone.
                if not line.startswith("#"):
                    reader = line.split()
                    # Color.
                    # Generate colors correctly <-- HARDCODED, FIX
                    if cmd_select in {2, 5, 9}:
                        isoch_col.append(float(reader[mag1_idx]) -
                                         float(reader[mag2_idx]))
                    else:
                        isoch_col.append(float(reader[mag2_idx]) -
                                         float(reader[mag1_idx]))
                    # Magnitude.
                    isoch_mag.append(float(reader[mag1_idx]))
                    # Mass
                    isoch_mas.append(float(reader[imass_idx]))

        # Save the last isochrone when EOF is reached.
        else:
            # If list is not empty.
            if isoch_col:
                # Store colors, magnitudes and masses for this
                # isochrone.
                metal_isoch.append([isoch_col, isoch_mag, isoch_mas])

    isoch = move_isoch(metal_isoch[0], e, d)

    return isoch


def get_asteca_params(cl):
    '''
    Return the metallicity, age, extinction and distance modulus for the
    'cl' cluster, obtained by ASteca.
    '''
    # Path to data file.
    out_file = 'asteca_output_final.dat'

    # Valid metallicity files names.
    met_vals = ['0.000100', '0.000121', '0.000152', '0.000191', '0.000241',
                '0.000303', '0.000382', '0.000481', '0.000605', '0.000762',
                '0.000959', '0.001207', '0.001520', '0.001914', '0.002409',
                '0.003033', '0.003818', '0.004807', '0.006051', '0.007618',
                '0.009591', '0.012074', '0.015200']
    # Read data file
    with open(out_file) as f:

        for line in skip_comments(f):
            l = line.split()
            if l[0] == cl:
                as_z, as_a, as_e, as_d = l[20], l[22], l[24], l[26]
                # Replace 0. values with minimum value.
                as_z = '0.0001' if float(as_z) < 0.0001 else as_z
                # Find closest metallicity values from list.
                z_idx = min(range(len(met_vals)),
                            key=lambda i: abs(float(met_vals[i]) -
                                              float(as_z)))
                # Pass valid string.
                as_z_str = met_vals[z_idx]
                return as_z, as_z_str, as_a, as_e, as_d


def save_crop_img(fig_name):
    '''
    Crop image.
    '''
    image = Image.open(fig_name)
    image_inv = ImageOps.invert(image.convert('RGB'))
    image_inv.load()
    # imageBox = image.getbbox()
    imageBox = image_inv.getbbox()
    cropped = image.crop(imageBox)
    fig_name_2 = fig_name[:-4] + '_crop.png'
    cropped.save(fig_name_2)
    os.remove(fig_name)
