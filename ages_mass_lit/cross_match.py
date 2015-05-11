
import numpy as np


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def get_asteca_data(i):
    '''
    Read the ASteCA output data file 'asteca_output.dat' and store each
    data column for each cluster.
    '''

    # Path to data file.
    out_file = '../asteca_output_' + i + '.dat'

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


def age_errors(gal, age):
    '''
    Assign errors to Hunter et al. (2003) age values.
    '''
    # Errors in logarithmic scale as defined in the article. Last value (0.15)
    # is for all cluster with ages outside of the defined ranges.
    errs_LMC_SMC = [[0.13, 0.13, 0.52, 0.54, 1.06, 0.17, 0.09, 0.15],
        [0.15, 0.31, 0.58, 0.24, 0.24, 0.16, 0.12, 0.15]]

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
            e_age = age_errors(gal, age)
            # Convert to log(age)
            log_age = np.log10(age * 10 ** 6)
            M_V_10 = float(lin[22])
            # Convert absolute magnitude to mass.
            mass = mag_2_mass(M_V_10)
            quality = float(lin[23])
            h03.append([gal, names, log_age, e_age, mass, quality])
            # raw_input()

    return h03


def slices(s, args):
    position = 0
    print args
    for length in args:
        yield s[position:position + length]
        position += length


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
            # lin = line.split()
            col_widths = [7, 5, 6, 8, 1, 7, 3, 6]
            print list(slices(line, col_widths))
            raw_input()
            # gal = lin[0]
            # # Store all names as separate uppercase strings.
            # names = [_.upper() for _ in lin[2].split(',')]
            # age = float(lin[21])
            # # Error in age is already in logarithmic scale.
            # e_age = age_errors(gal, age)
            # # Convert to log(age)
            # log_age = np.log10(age * 10 ** 6)
            # M_V_10 = float(lin[22])
            # # Convert absolute magnitude to mass.
            # mass = mag_2_mass(M_V_10)
            # quality = float(lin[23])
            # g10.append([gal, names, log_age, e_age, mass, quality])

    return g10


def match_clusts(as_names, as_pars, h03, g10):
    '''
    '''

    # Store H03, G10, P12 in each sub-list respectively.
    match_cl = [[[], [], []] for _ in len(as_names)]

    # Cross-match all clusters processed by ASteCA.
    for i, cl_n in enumerate(as_names):

        # Match clusters in H03.
        for cl_h in h03:
            # For each stored cluster name.
            for cl_h_n in cl_h[1]:
                # If names match.
                if cl_n == cl_h_n:
                    # Store H03 cluster data.
                    match_cl[i][0].append(cl_h[2], cl_h[3], cl_h[4], cl_h[5])
                    print cl_n, as_pars[i][21], cl_h[2], cl_h[3],\
                        as_pars[i][27], cl_h[4]

    return match_cl


def main():

    # Read ASteCA data.
    as_names, as_pars = get_asteca_data('0')

    # Read Hunter et al. (2003) data.
    # h03 = read_hunter()

    # Read Glatt  et. al (2010) data.
    g10 = read_glatt()
    print g10[:10]

    # Cross-match all clusters.
    match_cl = match_clusts(as_names, as_pars, h03, g10)


if __name__ == "__main__":
    main()
