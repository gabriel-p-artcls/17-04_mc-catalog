
import numpy as np


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


def mag_2_mass(M_V_10):
    '''
    Convert absolute magnitude at 10Myr to mass using Eq (1) in Hunter et al.
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
            age = float(lin[21])
            # Error in age is already in logarithmic scale.
            e_age = h03_age_errors(gal, age)
            # Convert to log(age)
            log_age = "{:.2f}".format(np.log10(age * 10 ** 6))
            M_V_10 = float(lin[22])
            # Convert absolute magnitude to mass.
            mass = "{:.2f}".format(mag_2_mass(M_V_10))
            quality = int(lin[23])
            h03.append([gal, names, log_age, e_age, mass, quality])

    # names_lst = []
    # for cl in h03:
    #     for name in cl[1]:
    #         names_lst.append(name)

    # print len(names_lst), names_lst, '\n'
    # j = 0
    # for i, n1 in enumerate(names_lst):
    #     for n2 in names_lst[(i + 1):]:
    #         if n1 == n2 and n1 != '':
    #             print j, i, n1
    #             j += 1
    # raw_input()

    return h03


def slices(s, args):
    '''
    Take a long string 's', a list with column widths 'args' and return the
    string sliced into as many smaller strings as column widths are passed.
    '''
    position = 0
    for length in args:
        yield s[position:position + length]
        position += length


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
            E_BV = "{:.2f}".format(float(lin[1]))
            log_age = "{:.2f}".format(float(lin[4]))
            q = int(lin[5])
            e_age = g10_age_errors(q)
            g10.append([gal, names, log_age, e_age, E_BV])

    return g10


def read_popescu_h03():
    '''
    Read Popescu et al. (2012) LMC database (correlated with H03 and G10).

    Return
    ------

    p12 = []
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
            log_age = "{:.2f}".format(float(lin[13]))
            # Obtain average error from upper and lower estimates.
            e_age = (float(lin[14]) + float(lin[15])) / 2.
            # Mass.
            mass = "{:.2f}".format(float(lin[16]))
            e_mass = (float(lin[17]) + float(lin[18])) / 2.
            p12.append([gal, names, log_age, e_age, mass, e_mass])

    return p12


def match_clusts(as_names, as_pars, h03, g10, p12):
    '''
    Cross match clusters processed by ASteCA to those published in several
    articles. The final list is ordered in the same way the 'as_params' list
    is.

    match_cl = [[[h03], [g10], [p12]], ..., N_clusts]
    h03 = [Gal, name, 'H03', log_age, e_age, log_age_asteca, e_age, mass,
    '--', mass_asteca, e_mass, '--', quality]

    g10 = [Gal, name, 'G10', log_age, e_age, log_age_asteca, e_age, '--', '--',
    '--', '--', E_BV, '--']

    p12 = [Gal, name, 'P12', log_age, e_age, log_age_asteca, e_age, mass,
    e_mass, mass_asteca, e_mass, '--', '--']

    '''

    # Store H03, G10, P12 in each sub-list respectively.
    match_cl = [[[], [], []] for _ in range(len(as_names))]

    # Cross-match all clusters processed by ASteCA.
    total = [0, 0, 0]
    for i, cl_n in enumerate(as_names):

        # Match clusters in H03.
        for cl_h in h03:
            # For each stored cluster name.
            for cl_h_n in cl_h[1]:
                # If names match.
                if cl_n == cl_h_n:
                    # Store H03 cluster data.
                    match_cl[i][0] = [cl_h[0], cl_n, 'H03', cl_h[2], cl_h[3],
                                      as_pars[i][21], as_pars[i][22], cl_h[4],
                                      '--', as_pars[i][27], as_pars[i][28],
                                      '--', cl_h[5]]
                    # Increase counter.
                    total[0] = total[0] + 1

        # Match clusters in G10.
        for cl_h in g10:
            # For each stored cluster name.
            for cl_h_n in cl_h[1]:
                # If names match.
                if cl_n == cl_h_n:
                    # Store G10 cluster data.
                    match_cl[i][1] = [cl_h[0], cl_n, 'G10', cl_h[2], cl_h[3],
                                      as_pars[i][21], as_pars[i][22], '--',
                                      '--', '--', '--', cl_h[4], '--']
                    total[1] = total[1] + 1

        # Match clusters in P12.
        for cl_h in p12:
            # For each stored cluster name.
            for cl_h_n in cl_h[1]:
                # If names match.
                if cl_n == cl_h_n:
                    # Store P12 cluster data.
                    match_cl[i][2] = [cl_h[0], cl_n, 'P12', cl_h[2], cl_h[3],
                                      as_pars[i][21], as_pars[i][22], cl_h[4],
                                      cl_h[5], as_pars[i][27], as_pars[i][28],
                                      '--', '--']
                    # Increase counter.
                    total[2] = total[2] + 1

    print '\nTotal clusters matched in each database:', total, \
        sum(_ for _ in total)

    return match_cl


def write_out_data(match_cl):
    '''
    Write matched clusters data to output file.
    '''

    # Read data file
    with open('matched_clusters.dat', 'w') as f_out:
        f_out.write("#\n# Age1: log(age)_DB\n# Age2: log(age)_asteca\n")
        f_out.write("#\n# Mass1: Mass_DB\n# Mass: Mass_asteca\n#\n")
        f_out.write("#GAL      NAME   DB   Age1  e_age  Age2  \
e_age      Mass1   e_mass    Mass2   e_mass     E_BV  quality\n#\n")
        for data_base in zip(*match_cl):
            for clust in data_base:
                if clust:  # Check that list is not empty.
                    f_out.write('''{:<5} {:>8} {:>4} {:>6} {:>6} {:>5} {:>6} \
{:>10} {:>8} {:>8} {:>8} {:>8} {:>8}\n'''.format(*[str(_) for _ in clust]))


def main():

    # Read ASteCA data.
    as_names, as_pars = get_asteca_data()

    # Read Hunter et al. (2003) data.
    h03 = read_hunter()

    # Read Glatt  et. al (2010) data.
    g10 = read_glatt()

    # Read Popescu  et. al (2012) data correlated with H03.
    p12 = read_popescu_h03()

    # Cross-match all clusters.
    match_cl = match_clusts(as_names, as_pars, h03, g10, p12)
    # print np.array(as_names[:10])
    # print np.array(match_cl[:10])

    # Write to file.
    write_out_data(match_cl)


if __name__ == "__main__":
    main()
