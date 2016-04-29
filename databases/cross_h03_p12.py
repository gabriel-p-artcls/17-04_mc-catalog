
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy import units as u


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def find_dup_cls_in_database(db_name, db):
    '''
    Check for duplicate clusters in the database.

    db = [gal, names, log_age, e_age, mass, e_mass, quality]
    (where 'quality' is present only for the H03 database)
    '''

    dup_names, dup_vals = [[] for _ in db], [[] for _ in db]
    indx_dup = []
    for i, el1 in enumerate(db):
        # If this cluster has already been stored as a duplicate, do not
        # process again.
        if i not in indx_dup:

            names1 = el1[1]
            for j, el2 in enumerate(db[(i + 1):]):
                names2 = el2[1]
                if any(n1 == n2 for n1 in names1 for n2 in names2):

                    # Store index of duplicated entry so that it won't be
                    # processed again by the outer loop.
                    indx_dup.append(j + i + 1)
                    # Store extra names for the cluster, not already stored.
                    for n1 in names1:
                        if n1 not in dup_names[i]:
                            dup_names[i].append(n1)
                    for n2 in names2:
                        if n2 not in dup_names[i]:
                            dup_names[i].append(n2)

                    # Store age, e_age, mass, and e_mass values.
                    dup_vals[i].append(db[j + i + 1][2:6:])
                    print db_name, i, db[i], j + i + 1, dup_names[i],\
                        dup_vals[i]

    print '\nEntries to be removed from {} ({}): {}'.format(
        db_name, len(indx_dup), np.asarray(sorted(indx_dup)[::-1]))

    print '\nAveraged values for duplicated clusters:'
    # Average values for all duplicated clusters in database.
    for i, vals in enumerate(dup_vals):
        # skip non-duplicated clusters.
        if vals:
            # Store all values for the cluster in a flat list.
            if len(vals) > 1:
                dup_ages = list(zip(*vals)[0])
                dup_e_ages = list(zip(*vals)[1])
                dup_masses = list(zip(*vals)[2])
                dup_e_masses = list(zip(*vals)[3])
            else:
                dup_ages = [vals[0][0]]
                dup_e_ages = [vals[0][1]]
                dup_masses = [vals[0][2]]
                dup_e_masses = [vals[0][3]]

            # Obtain the average for age and mass.
            age_avrg = np.mean([[db[i][2]] + dup_ages])
            mass_avrg = np.mean([[db[i][4]] + dup_masses])
            # print 'ages:', [db[i][2]] + dup_ages
            # print 'e_ages:', [db[i][3]] + dup_e_ages
            # print 'masses:', [db[i][4]] + dup_masses
            # print 'e_masses:', [db[i][5]] + dup_e_masses

            # Obtain the new errors as the maximum value between all the errors
            # and half the distance between the min and max value for each
            # parameter.
            age_half_dist = abs(max([db[i][2]] + dup_ages) -
                                min([db[i][2]] + dup_ages)) / 2.
            mass_half_dist = abs(max([db[i][4]] + dup_masses) -
                                 min([db[i][4]] + dup_masses)) / 2.
            # print dup_e_ages, age_half_dist
            # print dup_e_masses, mass_half_dist
            age_err = max(max([db[i][3]] + dup_e_ages), age_half_dist)
            mass_err = max(max([db[i][5]] + dup_e_masses), mass_half_dist)

            # Store age, e_age, mass, and e_mass final values for this cluster.
            # Replace in first appearance of list.
            db[i][2:6:] = [age_avrg, age_err, mass_avrg, mass_err]
            print db[i]

    # Remove all duplicated entries from database.
    for i in sorted(indx_dup)[::-1]:
        db.pop(i)

    return db


def mag_2_mass(M_V_10):
    '''
    Convert absolute magnitude at 10 Myr to mass using Eq (1) in Hunter et al.
    2003.
    '''
    return 10 ** (6. + 0.4 * (-14.55 - M_V_10))


def h03_age_errors(gal, age):
    '''
    Assign errors to Hunter et al. (2003) age values.
    '''
    # Errors in logarithmic scale as defined in the article. Last value (0.15)
    # is for all cluster with ages outside of the defined ranges.
    errs_LMC_SMC = [
        [0.13, 0.13, 0.52, 0.54, 1.06, 0.17, 0.09, 0.15],
        [0.15, 0.31, 0.58, 0.24, 0.24, 0.16, 0.12, 0.15]
    ]

    # Identify galaxy.
    g = 0 if gal == 'LMC' else 1

    # Identify age range (See Table 1 in H03).
    if age == 5:
        i = 0
    elif age == 10:
        i = 1
    elif age == 20:
        i = 2
    elif 30 <= age <= 60:
        i = 3
    elif 70 <= age <= 80:
        i = 4
    elif 90 <= age <= 100:
        i = 5
    elif 200 <= age <= 10000:
        i = 6
    else:
        i = 7

    # Assign error.
    age_err = errs_LMC_SMC[g][i]

    return age_err


def read_hunter():
    '''
    Read Hunter et al. (2003) database.

    Return
    ------

    h03 = [elem1, elem2, ...] <-- One element per cluster.
    elemX = [GAL, [names], age, mass, quality]
    GAL <-- 'LMC' or 'SMC'.
    [names] <-- List of names as strings.
    log(age) <-- Logarithmic age of the cluster.
    e_age <-- Age error.
    mass <-- Mass value for the cluster.
    -1. <-- Since there's no mass error assigned.
    quality <-- quality flag: 0=good, 1=probable, 2=questionable
    '''

    # Path to data file.
    h03_file = 'hunter_03.dat'

    # Read data file
    with open(h03_file) as f:

        h03 = []
        for line in skip_comments(f):
            lin = line.split()
            gal = lin[0]
            # Store all names as separate uppercase strings.
            names = [_.upper() for _ in lin[2].split(',')]

            # Remove empty entries ('') and ESO* clusters that introduce
            # spurious duplicates.
            idx_rm = []
            for i, name in enumerate(names):
                if name == '':
                    idx_rm.append(i)
                if name[:3] == 'ESO':
                    idx_rm.append(i)
            if idx_rm:
                for i in idx_rm[::-1]:
                    names.pop(i)

            age = float(lin[21])
            # Error in age is already in logarithmic scale.
            e_age = h03_age_errors(gal, age)
            # Convert to log(age)
            log_age = np.log10(age * 10 ** 6)
            M_V_10 = float(lin[22])
            # Convert absolute magnitude to mass.
            mass = mag_2_mass(M_V_10)
            c = SkyCoord(ra=Angle(str(lin[3]) + ' hours'),
                         dec=Angle(str(lin[4]) + ' degrees'), frame='icrs')
            h03.append([gal, names, log_age, e_age, mass, -1.,
                       [c.ra.deg, c.dec.deg]])

    h03 = find_dup_cls_in_database('H03', h03)

    return h03


def read_popescu():
    '''
    Read Popescu et al. (2012) LMC database (correlated with H03 and G10).

    Return
    ------

    p12 = [gal, names, log_age, e_age, mass, e_mass]
    '''

    # Path to data file.
    p12_h03_file = 'popescu_12_LMC.dat'
    p12_g10_file = 'popescu_12_LMC_glatt.dat'

    p12 = []

    # Read data file
    with open(p12_h03_file) as f:

        for line in skip_comments(f):
            lin = line.split()
            gal = 'LMC'
            # Store all names as separate uppercase strings.
            names = [_.upper() for _ in lin[0].split('.')]
            log_age = float(lin[10])
            # Obtain average error from upper and lower estimates.
            e_age = (float(lin[11]) + float(lin[12])) / 2.
            # Mass.
            mass = float(lin[13])
            e_mass = (float(lin[14]) + float(lin[15])) / 2.
            ra, dec = line[143:].split()
            c = SkyCoord(ra=float(ra)*u.degree,
                         dec=float(dec)*u.degree, frame='icrs')
            p12.append([gal, names, log_age, e_age, mass, e_mass,
                        [c.ra.deg, c.dec.deg]])

    # Read data file
    with open(p12_g10_file) as f:

        for line in skip_comments(f):
            lin = line.split()
            gal = 'LMC'
            # Store all names as separate uppercase strings.
            names = [_.upper() for _ in lin[0].split('.')]
            log_age = float(lin[13])
            # Obtain average error from upper and lower estimates.
            e_age = (float(lin[14]) + float(lin[15])) / 2.
            # Mass.
            mass = float(lin[16])
            e_mass = (float(lin[17]) + float(lin[18])) / 2.
            ra, dec = line[163:].split()
            c = SkyCoord(ra=float(ra)*u.degree,
                         dec=float(dec)*u.degree, frame='icrs')
            p12.append([gal, names, log_age, e_age, mass, e_mass,
                       [c.ra.deg, c.dec.deg]])

    p12 = find_dup_cls_in_database('P12', p12)

    # # Check for duplicates.
    # dup_found = False
    # for i, cl1 in enumerate(p12):
    #     names1 = cl1[1]
    #     if any(n1 == n2 for n1 in names1 for n2 in names):
    #         dup_found = True
    #         dup_indx = i
    #         # Store values for duplicated cluster.
    #         dup_cls.append()

    # if dup_found:
    #     # Average values.
    #     print 'P12 dup:', cl1, clust_data
    #     a = [cl1[2:], clust_data[2:]]
    #     clust_data[2:] = map(mean_lsts, zip(*a))
    #     print 'P12 avrg:', clust_data
    #     p12[dup_indx] = clust_data
    # else:
    #     p12.append(clust_data)

    return p12


def remove_duplicates(values):
    output = []
    seen = set()
    for value in values:
        # If value has not been encountered yet,
        # ... add it to both list and set.
        if value not in seen:
            output.append(value)
            seen.add(value)
    return output


def write_out_data(h03, p12, idx, d2d):
    '''
    Write matched clusters data to output file.
    '''

    match_cl = []
    for i, j in enumerate(idx):
        # Reject matches with a distance of more than X arcsec.
        if d2d[i].deg < 0.00556:  # 20 arcsec ~ 0.00556
            names = ', '.join(remove_duplicates(p12[i][1] + h03[j][1]))
            match_cl.append(['LMC', p12[i][-1][0], p12[i][-1][1],
                            h03[j][2], h03[j][3], p12[i][2], p12[i][3],
                            h03[j][4], -1., p12[i][4], p12[i][5], names])

    # Read data file
    with open('matched_H03_P12.dat', 'w') as f_out:
        f_out.write("#\n# Age1: log(age)_H03\n# Age2: log(age)_P12\n")
        f_out.write("#\n# Mass1: Mass_H03\n# Mass2: Mass_P12\n#\n")
        f_out.write("#GAL   RA        DEC         Age1   e_age  Age2   "
                    "e_age   Mass1   e_m1  Mass2    e_m2     NAMES\n")
        for clust in match_cl:
            f_out.write(
                "{:<4} {:>10.5f} {:>10.5f} {:>6.2f} {:>6.2f} "
                "{:>6.2f} {:>6.2f} {:>8.0f} {:>3.0f} {:>8.0f} "
                "{:>8.0f}     {}\n".format(*clust))


def main():
    """
    Cross match Hunter et al. (2003) and Popescu et al. (2012) databases.
    """

    # Read Hunter et al. (2003) data.
    h03 = read_hunter()
    # Read Popescu et al. (2012) data.
    p12 = read_popescu()

    # Store ra,dec coordinates.
    ra1, dec1 = zip(*zip(*p12)[-1])
    ra2, dec2 = zip(*zip(*h03)[-1])
    # Create catalogues in decimal degrees.
    cp12 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
    ch03 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
    # Cross-match all clusters.
    # http://docs.astropy.org/en/stable/coordinates/matchsep.html
    # idx are indices into ch03 that are the closest objects to each of the
    # coordinates in cp12.
    # d2d are the on-sky distances between them.
    idx, d2d, d3d = cp12.match_to_catalog_sky(ch03)

    # Write cross-matched data to file.
    # write_out_data(h03, p12, idx, d2d)

    # Some stats on this data.
    print '\n\n'
    n = 0
    m_le, m_gt = [[], []], [[], []]
    for i, j in enumerate(idx):
        if d2d[i].deg < 0.00556:  # 20 arcsec ~ 0.00556
            n += 1
            if .5*(h03[j][4]+p12[i][4]) <= 500000.:
                m_le[0].append(h03[j][4])
                m_le[1].append(p12[i][4])
            else:
                m_gt[0].append(h03[j][4])
                m_gt[1].append(p12[i][4])
            if .5*(h03[j][4]+p12[i][4]) > 100000.:
                print p12[i][1], p12[i][-1][0], p12[i][-1][1],\
                    h03[j][2], p12[i][2],  h03[j][4], p12[i][4]
        else:
            # OCs in P12 not matched in H03.
            # print d2d[i].deg, p12[i]
            pass

    print '\nOCs in H03, P12:', len(h03), len(p12)
    print 'OCs matched:', n

    f, (ax1, ax2) = plt.subplots(1, 2)

    ax1.set_xlabel('0.5*(P12+H03)>5000')
    ax1.set_ylabel('H03-P12')
    ax1.scatter(np.log((np.array(m_le[0])+np.array(m_le[1]))/2.),
                np.array(m_le[0]) - np.array(m_le[1]), c='r')

    ax2.set_xlabel('0.5*(P12+H03)<5000')
    ax2.set_ylabel('H03-P12')
    ax2.scatter(np.log((np.array(m_le[0])+np.array(m_le[1]))/2.),
                (np.array(m_le[0]) - np.array(m_le[1])) /
                (np.array(m_le[0])+np.array(m_le[1])), c='b')

    plt.show()

if __name__ == "__main__":
    main()
