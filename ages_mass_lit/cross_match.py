
from pyexcel_ods import ODSBook
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def get_asteca_data():
    '''
    Read the ASteCA output data file 'asteca_output_final.dat' and store each
    data column for each cluster.
    '''

    # Path to data file.
    out_file = '../asteca_output_final.dat'

    # Read data file
    with open(out_file) as f:
        as_names, as_pars = [], []

        for line in skip_comments(f):
            as_names.append(line.split()[0])
            # Read clusters parameters obtained by ASteCA.
            as_pars.append(line.split()[1:])

    return as_names, as_pars


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


def slices(s, args):
    '''
    Take a long string 's', a list with column widths 'args' and return the
    string sliced into as many smaller strings as column widths are passed.
    '''
    position = 0
    for length in args:
        yield s[position:position + length]
        position += length


def get_liter_data():
    '''
    Read the data file with the literature values for each cluster as a
    dictionary.
    '''

    # Read .ods file with literature data.
    cl_file = ODSBook('../lista_unica_cumulos.ods')
    # Store as dictionary and then as list.
    cl_dict = cl_file.sheets()["S-LMC"]

    # Indexes of coord, age and extinction columns in .ods literature file.
    ra_i, dec_i, age_i, e_age_i, ext_i, e_ext_i = cl_dict[0].index(u'ra_deg'),\
        cl_dict[0].index(u'dec_deg'), cl_dict[0].index(u'log(age)'), \
        cl_dict[0].index(u'e_log(age)'), cl_dict[0].index(u'E(B-V) (lit)'), \
        cl_dict[0].index(u'e_E(B-V)')
    # Index of the cluster's name in the .ods file.
    name_idx = cl_dict[0].index(u'Name')

    names_ra_dec, ra, dec, ages, e_age, exti, e_exti = [], [], [], [], [], [],\
        []
    for cl in cl_dict:
        names_ra_dec.append(str(cl[name_idx]))
        ra.append(cl[ra_i])
        dec.append(cl[dec_i])
        ages.append(cl[age_i])
        e_age.append(cl[e_age_i])
        exti.append(cl[ext_i])
        e_exti.append(cl[e_ext_i])
    # remove first line (column names) and 4 last lines (empty string)
    del names_ra_dec[-4:], ra[-4:], dec[-4:], ages[-4:], e_age[-4:], \
        exti[-4:], e_exti[-4:]
    del names_ra_dec[0], ra[0], dec[0], ages[0], e_age[0], exti[0], e_exti[0]

    # Create the RA, DEC catalog.
    cat_ra_dec = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')

    return names_ra_dec, cat_ra_dec, ages, e_age, exti, e_exti


def match_ra_dec_asteca(names_ra_dec, cat_ra_dec, ra, dec):
    '''
    Receive cluster center (ra, dec) coordinates in decimal degrees and use
    them to match with the closest cluster in the ASteCA database, within
    some predefined tolerance.
    '''

    # Store (ra, dec) as valid coordinate object.
    cl_coord = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
    # Find closest match in ASteCA catalog to (ra, dec) coordinates.
    i, sep2d, dist3d = match_coordinates_sky(cl_coord, cat_ra_dec,
                                             nthneighbor=1)

    # Distance to closest match in degrees.
    dist_deg = float(sep2d[0].to_string(decimal=True))

    # Match within a given tolerance.
    # 1 arcsec ~ 0.000278 deg
    # if dist_deg < 0.004167:  # 15 arcsec ~ 0.004167 (25 clusters)
    # if dist_deg < 0.00833:  # 30 arcsec ~ 0.00833 (26 clusters)
    # if dist_deg < 0.0167:  # 1 arcmin ~ 0.0167 (27 clusters)
    # if dist_deg < 0.002778:  # 10 arcsec ~ 0.002778 (24 clusters)
    if dist_deg < 0.00556:  # 20 arcsec ~ 0.00556
        name = str(names_ra_dec[i])
    else:
        name = ''

    return name, dist_deg


def read_pietr99(names_ra_dec, cat_ra_dec):
    '''
    Read Pietrzynski et al. (1999) database.

    Return
    ------

    p99 = []
    '''

    # Path to data file.
    p99_file = 'pietrz_99_SMC.dat'

    # Read data file
    with open(p99_file) as f:
        p99 = []

        for line in skip_comments(f):
            # Width of columns in file.
            col_widths = [7, 13, 14, 5, 7, 5, 4]
            lin = list(slices(line, col_widths))
            # convert coords to decimal degrees.
            c = SkyCoord(lin[1] + lin[2], unit=(u.hourangle, u.deg))
            # print c.ra.deg, c.dec.deg
            # Find match in ASteCA database.
            name, dist_deg = match_ra_dec_asteca(names_ra_dec, cat_ra_dec,
                                                 c.ra.deg, c.dec.deg)
            # Only append if a match was found.
            if name:
                print 'P99 match: ', name, c.ra.deg, c.dec.deg, dist_deg
                gal = 'SMC'
                E_BV = float(lin[4])
                log_age = float(lin[5])
                e_age = float(lin[6])
                p99.append([gal, [name], log_age, e_age, E_BV])

    return p99


def read_pietr(names_ra_dec, cat_ra_dec):
    '''
    Read Pietrzynski et al. (2000) database.

    Return
    ------

    p00 = []
    '''

    # Path to data file.
    p00_file = 'pietrz_00_LMC.dat'

    # Read data file
    with open(p00_file) as f:
        p00 = []

        for line in skip_comments(f):
            # Width of columns in file.
            col_widths = [7, 13, 14, 9, 6, 6, 11]
            lin = list(slices(line, col_widths))
            # convert coords to decimal degrees.
            c = SkyCoord(lin[1] + lin[2], unit=(u.hourangle, u.deg))
            # print c.ra.deg, c.dec.deg
            # Find match in ASteCA database.
            name, dist_deg = match_ra_dec_asteca(names_ra_dec, cat_ra_dec,
                                                 c.ra.deg, c.dec.deg)
            # Only append if a match was found.
            if name:
                print 'P00 match: ', name, c.ra.deg, c.dec.deg, dist_deg
                gal = 'LMC'
                log_age = float(lin[4])
                e_age = float(lin[5])
                p00.append([gal, [name], log_age, e_age])

    return p00


def read_rafel(names_ra_dec, cat_ra_dec):
    '''
    Read Rafelski et al. (2005) database (use GALEV model with z=0.004).

    Return
    ------

    r05 = []
    '''

    # Path to data file.
    r05_file = 'rafelski_05_SMC.dat'

    # Read data file
    with open(r05_file) as f:
        r05 = []

        for line in skip_comments(f):
            # Width of columns in file.
            col_widths = [5, 14, 16, 58, 10, 7, 28, 7, 7, 7, 10, 7, 28, 7,
                          7, 7, 10, 7, 28, 7, 7, 7, 10, 7, 28, 7, 7, 7, 10, 7,
                          28, 7, 7, 7]
            lin = list(slices(line, col_widths))
            # convert coords to decimal degrees.
            c = SkyCoord(lin[1] + lin[2], unit=(u.hourangle, u.deg))
            # print c.ra.deg, c.dec.deg
            # Find match in ASteCA database.
            name, dist_deg = match_ra_dec_asteca(names_ra_dec, cat_ra_dec,
                                                 c.ra.deg, c.dec.deg)

            # Only append if a match was found and the cluster has an age
            # assigned.
            # Use GALEv z=0.004 model age for testing if an age was assigned.
            age = float(lin[7])
            if name and age < 99999:
                print 'R05 match: ', name, c.ra.deg, c.dec.deg, dist_deg
                # For each defined age for all the models.
                age_values, age_errors = [], [[], []]
                for i_age in [7, 13, 19, 25, 31]:
                    age = age = float(lin[i_age])
                    er_age, ER_age = float(lin[i_age + 1]),\
                        float(lin[i_age + 2].rstrip('\n'))
                    # Store values.
                    age_values.append(age)
                    age_errors[0].append(er_age)
                    age_errors[1].append(ER_age)

                # Average ages and obtain error as the midpoint between the
                # minimum lower error and the maximum upper error.
                age_mean = np.mean(age_values)
                # Obtain error midpoint.
                lo_bound = np.array(age_values) - np.array(age_errors[0])
                up_bound = np.array(age_values) + np.array(age_errors[1])
                e_m_age = (max(up_bound) - min(lo_bound)) / 2.
                # Obtain log age error.
                e_age = (e_m_age / age_mean) * (1. / np.log(10))
                log_age = np.log10(age_mean * 10 ** 6)

                print name
                print age_values
                print age_errors[0]
                print age_errors[1]
                print min(lo_bound), max(up_bound)
                print age_mean, e_m_age
                print log_age, e_age, '\n'

                gal = 'SMC'
                r05.append([gal, [name], log_age, e_age])

    return r05


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
            h03.append([gal, names, log_age, e_age, mass, -1.])

    h03 = find_dup_cls_in_database('H03', h03)

    return h03


def c06_age_errors(c):
    '''
    Assign errors to Chiosi et al. (2006) age values.
    '''
    # Errors in logarithmic scale as defined in the article
    errs = [0.3, 0.4, 0.5]

    # Identify \delta log(age) range.
    if c == 1:
        i = 0
    elif c == 2:
        i = 1
    elif c == 3:
        i = 2

    # Assign error.
    age_err = errs[i]

    return age_err


def read_chiosi():
    '''
    Read Chiosi et al. (2006) LMC database.

    Return
    ------

    c06 = [gal, names, log_age, e_age, E_BV, t]
    '''

    # Path to data file.
    c06_file = 'chiosi_06.dat'

    c06 = []

    # Read data file
    with open(c06_file) as f:

        for line in skip_comments(f):
            lin = line.split()
            gal = 'SMC'
            # Store all names as separate uppercase strings.
            names = [_.upper() for _ in lin[8].split(',')]
            log_age = float(lin[4])
            c = int(lin[7])
            e_age = c06_age_errors(c)
            # E(V-I) = 1.244*E(B-V)
            E_BV = float(lin[5]) / 1.244
            c06.append([gal, names, log_age, e_age, E_BV])

    return c06


def g10_age_errors(q):
    '''
    Assign errors to Glatt et al. (2010) age values.
    '''
    # Errors in logarithmic scale as defined in the article
    errs = [0.3, 0.4, 0.5]

    # Identify \delta log(age) range.
    if q == 1:
        i = 0
    elif q == 2:
        i = 1
    elif q == 3:
        i = 2
    elif q == 9:
        i = 2

    # Assign error.
    age_err = errs[i]

    return age_err


def read_glatt():
    '''
    Read Glatt et al. (2010) database.

    Return
    ------

    g10 = []
    '''

    # Path to data file.
    g10_file = 'glatt_10.dat'

    # Read data file
    with open(g10_file) as f:
        g10 = []

        for line in skip_comments(f):
            # Width of columns in file.
            col_widths = [12, 6, 8, 1, 7, 3, 62, 150]
            lin = list(slices(line, col_widths))
            gal = lin[0][:3]
            names = [_.upper().strip() for _ in lin[7].split(',')]
            E_BV = float(lin[1])
            log_age = float(lin[4])
            q = int(lin[5])
            e_age = g10_age_errors(q)
            g10.append([gal, names, log_age, e_age, E_BV])

    return g10


# def mean_lsts(a):
#     return sum(a) / len(a)


def read_popescu_h03():
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
            p12.append([gal, names, log_age, e_age, mass, e_mass])

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
            p12.append([gal, names, log_age, e_age, mass, e_mass])

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


def match_clusts(as_names, as_pars, names_lit, lit_ages, lit_e_age, lit_ext,
                 lit_e_ext, p99, p00, h03, r05, c06, g10, p12):
    '''
    Cross match clusters processed by ASteCA to those published in several
    articles. The final list is ordered in the same way the 'as_params' list
    is.

    match_cl = [[[p00], [p00], [h03], [g10], [p12]], ..., N_clusts]

    DB = ['XXX', Gal, name, log_DB_age, e_DB_age, log_age_asteca, e_age_asteca,
    log_age_lit, e_log_age_lit, mass_DB, e_mass_DB, mass_asteca, e_mass_asteca,
    E_BV_DB, E_BV_asteca, e_E_BV_asteca, E_BV_lit, e_E_BV_lit]

    '''

    # Store all databases in a sub-list.
    match_cl = [[[], [], [], [], [], [], []] for _ in range(len(as_names))]

    # Cross-match all clusters processed by ASteCA.
    total = [0, 0, 0, 0, 0, 0, 0]
    db_names = ['P99', 'P00', 'H03', 'R05', 'C06', 'G10', 'P12']
    # For each cluster processed by ASteCA.
    for i, cl_n in enumerate(as_names):

        # For each database.
        for k, db in enumerate([p99, p00, h03, r05, c06, g10, p12]):
            # For each cluster cross-matched and stored in database.
            for cl_db in db:
                # For each cluster name associated to each cluster.
                for cl_h_n in cl_db[1]:

                    # Only H03 and P12 have defined masses.
                    if k in [2, 6]:
                        m_DB, e_m_DB = cl_db[4], cl_db[5]
                    else:
                        m_DB, e_m_DB = -1., -1.
                    # Only P99, C06 and G10 have defined extinctions.
                    if k in [0, 4, 5]:
                        ext_DB = cl_db[4]
                    else:
                        ext_DB = -1.

                    # If names match.
                    if cl_n == cl_h_n:
                        # For each cluster with literature values (Piatti et
                        # al.)
                        for j, cl_l in enumerate(names_lit):
                            # Check match in literature.
                            if cl_l == cl_n:
                                l_ext, l_age = lit_ext[j], lit_ages[j]
                                l_e_age, l_e_ext = lit_e_age[j], lit_e_ext[j]

                                # Convert '--' strings into -1. values.
                                l_ext = -1. if isinstance(l_ext, basestring) \
                                    else l_ext
                                l_e_ext = -1. if \
                                    isinstance(l_e_ext, basestring) else \
                                    l_e_ext

                                # Store cluster data.
                                match_cl[i][k] =\
                                    [db_names[k], cl_db[0], cl_n, cl_db[2],
                                        cl_db[3], as_pars[i][21],
                                        as_pars[i][22], l_age, l_e_age, m_DB,
                                        e_m_DB, as_pars[i][27], as_pars[i][28],
                                        ext_DB, as_pars[i][23], as_pars[i][24],
                                        l_ext, l_e_ext]
                                # Increase counter.
                                total[k] = total[k] + 1

    print '\nTotal clusters matched in each database:', total, \
        sum(_ for _ in total)

    return match_cl


def write_out_data(match_cl):
    '''
    Write matched clusters data to output file.
    '''

    # Read data file
    with open('matched_clusters.dat', 'w') as f_out:
        f_out.write("#\n# Age1: log(age)_DB\n# Age2: log(age)_asteca\n"
                    "# Age3: log(age)_literature\n")
        f_out.write("#\n# E_BV1: E_BV_DB\n# E_BV2: E_BV_asteca\n"
                    "# E_BV3: E_BV_lit\n")
        f_out.write("#\n# Mass1: Mass_DB\n# Mass2: Mass_asteca\n#\n")
        f_out.write("#DB   GAL      NAME   Age1  e_age  Age2  \
e_age   Age3  e_age      Mass1   e_mass    Mass2   e_mass    E_BV1    \
E_BV2   e_E_BV    E_BV3   e_E_BV\n")
        for data_base in zip(*match_cl):
            for clust in data_base:
                if clust:  # Check that list is not empty.
                    f_out.write('''{:<4} {:>4} {:>9} {:>6.2f} {:>6.2f} {:>5} \
{:>6} {:>6.2f} {:>6.2f} {:>10.2f} {:>8.0f} {:>8} {:>8} {:>8.2f} {:>8} {:>8}\
 {:>8.2f} {:>8.2f}\n'''.format(*clust))


def main():

    # Read ASteCA data.
    as_names, as_pars = get_asteca_data()

    # Read RA & DEC literature data.
    names_ra_dec, cat_ra_dec, lit_ages, lit_e_age, lit_ext, lit_e_ext = \
        get_liter_data()

    # Read Pietrzynski et al. (1999) data.
    p99 = read_pietr99(names_ra_dec, cat_ra_dec)

    # Read Pietrzynski et al. (2000) data.
    p00 = read_pietr(names_ra_dec, cat_ra_dec)

    # Read Hunter et al. (2003) data.
    h03 = read_hunter()

    # Read Rafelski et al. (2005) data.
    r05 = read_rafel(names_ra_dec, cat_ra_dec)

    # Read Chiosi et al. (2006) data.
    c06 = read_chiosi()

    # Read Glatt  et al. (2010) data.
    g10 = read_glatt()

    # Read Popescu et al. (2012) data.
    p12 = read_popescu_h03()

    # Cross-match all clusters.
    match_cl = match_clusts(as_names, as_pars, names_ra_dec, lit_ages,
                            lit_e_age, lit_ext, lit_e_ext, p99, p00, h03, r05,
                            c06, g10, p12)

    # Write to file.
    write_out_data(match_cl)


if __name__ == "__main__":
    main()
