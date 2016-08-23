

import numpy as np
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


def mar_par_data():
    """
    """
    # Path to data files.
    mar, par = 'mar2008_ubvrijhk/', 'parsec11_ubvrijhk/'
    age_values = [7.5, 8.0, 8.5, 9.0, 9.5]

    # Read data file
    mar_data, par_data = [], []
    for tr, tracks in enumerate([mar, par]):
        path = os.path.join(os.path.dirname(__file__), tracks)
        for m_f in os.listdir(path):
            met_f = path + m_f
            metal_isoch = read_met_file(met_f, age_values, tr)
            if tr == 0:
                mar_data.append(metal_isoch)
            else:
                par_data.append(metal_isoch)

    return mar_data, par_data


if __name__ == "__main__":
    mar_par_data()
