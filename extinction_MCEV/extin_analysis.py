
import numpy as np


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def get_data():
    '''
    Read RA, DEC, extinction values from table.
    '''

    tab_name = 'ra_dec_exts_mult_matches.dat'

    with open(tab_name, 'r') as f:
        ext_pars = []

        for line in skip_comments(f):
            # Skip coords witn no values.
            if line.split()[0] != '--':
                ra_c, dec_c, ra, dec = float(line.split()[0]), \
                    float(line.split()[1]), float(line.split()[4]), \
                    float(line.split()[5])
                # Convert extinction from E(V-I) to E(B-V)
                E_BV, e_EBV = float(line.split()[2]) / 1.38, \
                    float(line.split()[3]) / 1.38
                # Append all values.
                ext_pars.append([ra_c, dec_c, E_BV, e_EBV, ra, dec])

    # Return zipped list.
    ext_zip = zip(*ext_pars)

    return ext_zip


def match_coords(ext_pars):
    '''
    Store all the indexes pointing to each cluster in separate lists.
    '''

    coords_match = [[] for _ in range(210)]

    ra_old, dec_old, st_indx = 0., 0., -1
    for i, ra in enumerate(ext_pars[4]):
        dec = ext_pars[5][i]

        if ra_old == ra and dec_old == dec:
            coords_match[st_indx].append(i)
        else:
            ra_old, dec_old, st_indx = ra, dec, st_indx + 1
            coords_match[st_indx].append(i)

    return coords_match


def get_ext_values(ext_pars, coords_match):
    '''
    Obtain the closest extinction value and its distance (in degrees),
    the average extinction and its standard deviation, and the maximum
    extinction value.
    Convert from E(V-I) to E(B-V) according to:
    E(V-I) = 1.38 *  E(B-V)
    '''

    clusts_exts = []
    # Iterate for every cluster.
    for clust in coords_match:

        # Get original ra, dec for this cluster.
        ra, dec = ext_pars[4][clust[0]], ext_pars[5][clust[0]]

        # Get closest extinction value and its distance.
        closest_idx, dist_min = 0, 1000.
        for ext_idx in clust:
            ra_c, dec_c = ext_pars[0][ext_idx], ext_pars[1][ext_idx]
            # Get angular distance.
            cos_d = np.sin(np.deg2rad(dec)) * np.sin(np.deg2rad(dec_c)) + \
                np.cos(np.deg2rad(dec)) * np.cos(np.deg2rad(dec_c)) * \
                np.cos(np.deg2rad(ra - ra_c))
            # Arccos, radians to decimal degrees.
            dist = np.rad2deg(np.arccos(cos_d))
            if dist < dist_min:
                dist_min = dist
                closest_idx = ext_idx

        # Get average extinction value and standard deviation.
        avr_ext_lst = []
        for ext_idx in clust:
            avr_ext_lst.append(ext_pars[2][ext_idx])
        avr_ext, std_dev = np.mean(avr_ext_lst), np.std(avr_ext_lst)

        # Get maximum extinction value.
        max_ext, ext_old = 0., -1000.
        for ext_idx in clust:
            if ext_pars[2][ext_idx] > ext_old:
                ext_old = ext_pars[2][ext_idx]
                max_ext = ext_pars[2][ext_idx]

        # Store all values.
        clusts_exts.append([ra, dec, ext_pars[2][closest_idx], dist_min,
            avr_ext, std_dev, max_ext])

    return clusts_exts


def print_file(clusts_exts):
    '''
    Print data to file.
    '''

    with open('cls_exts_match.dat', 'w') as f:
        f.write("#RA_(deg)       DEC_(deg)    E_BV_close  dist(deg)  E_BV_avrg"
        "  E_BV_std_dev  E_BV_max\n")
        for line in clusts_exts:
            f.write("{:<15} {:<15} {:>8.3f} {:>8.3f} {:>8.3f} {:>8.3f} "
            "{:>8.3f}\n".format(*line))


def main():
    '''
    Take the table downloaded from the MCEV service
    (http://dc.zah.uni-heidelberg.de/mcextinct/q/cone/form) and match the
    extinction coordinates with the clusters coordinates.
    '''

    # Get data from table.
    ext_pars = get_data()

    # Match values for a single coordinate.
    coords_match = match_coords(ext_pars)

    # Get several ext values for each star.
    clusts_exts = get_ext_values(ext_pars, coords_match)

    # Print to file.
    print_file(clusts_exts)

    print '\nEnd.'


if __name__ == "__main__":
    main()
