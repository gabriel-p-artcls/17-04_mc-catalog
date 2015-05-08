
import numpy as np


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


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
    elemX = [GAL, [names], age, mass]
    GAL <-- 'LMC' or 'SMC'.
    [names] <-- List of names as strings.
    mass <-- Mass value for the cluster.
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


def main():
    h03 = read_hunter()
    as_names, as_pars = get_asteca_data('0')

    for i, cl in enumerate(as_names):
        for cl_h in h03:
            for cl_h_n in cl_h[1]:
                if cl == cl_h_n:
                    print cl, as_pars[i][21], cl_h[2], cl_h[3], as_pars[i][27], cl_h[4]


if __name__ == "__main__":
    main()
