
from pyexcel_ods import ODSBook


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def get_asteca_data():
    '''
    Read the ASteCA output data file 'asteca_output.dat' and store each
    data column for each cluster.
    '''

    # Path to data file.
    out_file = 'asteca_output_final.dat'

    # Read data file
    with open(out_file) as f:
        as_names, as_pars = [], []

        for line in skip_comments(f):
            as_names.append(line.split()[0])
            # Read clusters parameters obtained by ASteCA.
            as_pars.append(line.split()[1:])

    return as_names, as_pars


def get_bica_database():
    '''
    Read Bica et el. (2008) database.
    '''
    # Path to data file.
    bb_file = 'databases/bb_cat.dat'

    # Read data file
    with open(bb_file) as f:
        bica_database = []
        for line in skip_comments(f):
            # Read coordinates.
            bica_database.append([float(_) for _ in line.split()])

    return bica_database


def get_amr_lit():
    '''
    Read AMR literature data.
    '''
    # Path to data Pagel & T (1998) file.
    in_file = 'AMRs/PT98_s-lmc_amr.dat'
    # Read data file
    with open(in_file) as f:
        # i=0 --> SMC ; i=1 --> LMC
        amr_PT98_t = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            if l[2] != '--':
                amr_PT98_t[0].append([float(_) for _ in [l[2], l[3]]])
            amr_PT98_t[1].append([float(_) for _ in [l[0], l[1]]])

    # Put SMC values first and LMC second.
    amr_smc_PT98, amr_lmc_PT98 = zip(*amr_PT98_t[0]), zip(*amr_PT98_t[1])

    # Path to data G98 file.
    in_file = 'AMRs/G98_lmc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_lmc_G98 = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_lmc_G98[0].append(float(l[0]))
            amr_lmc_G98[1].append(float(l[1]))

    # Path to data C08 file.
    in_file = 'AMRs/C08_lmc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_lmc_C08 = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_lmc_C08[0].append(float(l[0]))
            amr_lmc_C08[1].append(float(l[1]))

    # Path to data C08 file.
    in_file = 'AMRs/C08_smc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_smc_C08 = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_smc_C08[0].append(float(l[0]))
            amr_smc_C08[1].append(float(l[1]))

    # Path to data HZ09 file.
    in_file = 'AMRs/HZ09_lmc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_lmc_HZ09 = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_lmc_HZ09[0].append(float(l[0]))
            amr_lmc_HZ09[1].append(float(l[1]))

    # Path to data R12 file.
    in_file = 'AMRs/R12_lmc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_lmc_R12 = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_lmc_R12[0].append(float(l[0]))
            amr_lmc_R12[1].append(float(l[1]))

    # Path to data HZ04 file.
    in_file = 'AMRs/HZ04_smc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_smc_HZ04 = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_smc_HZ04[0].append(float(l[0]))
            amr_smc_HZ04[1].append(float(l[1]))

    # Path to data N09 file.
    in_file = 'AMRs/N09_smc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_smc_N09 = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_smc_N09[0].append(float(l[0]))
            amr_smc_N09[1].append(float(l[1]))

    # Path to data TB09 file.
    in_file = 'AMRs/TB09_smc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_smc_TB09_1, amr_smc_TB09_2, amr_smc_TB09_3 = [[], []], [[], []],\
            [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_smc_TB09_1[0].append(float(l[0]))
            amr_smc_TB09_1[1].append(float(l[1]))
            amr_smc_TB09_2[0].append(float(l[0]))
            amr_smc_TB09_2[1].append(float(l[2]))
            amr_smc_TB09_3[0].append(float(l[0]))
            amr_smc_TB09_3[1].append(float(l[3]))

    # Path to data C13 file.
    in_file = 'AMRs/C13_smc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_smc_C13_B, amr_smc_C13_C = [[], []], [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            if l[0] != '--':
                amr_smc_C13_B[0].append(float(l[0]))
                amr_smc_C13_B[1].append(float(l[1]))
            amr_smc_C13_C[0].append(float(l[2]))
            amr_smc_C13_C[1].append(float(l[3]))

    # Path to data Piatti & Geisler (2003) file.
    in_file = 'AMRs/PG13_s-lmc_amr.dat'
    # Read data file
    with open(in_file) as f:
        # i=0 --> SMC ; i=1 --> LMC
        amr_PG13_t = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            amr_PG13_t[0].append([float(_) for _ in [l[4], l[6]]])
            amr_PG13_t[1].append([float(_) for _ in [l[0], l[2]]])

    # Path to data M14 file.
    in_file = 'AMRs/M14_lmc_amr.dat'
    # Read data file
    with open(in_file) as f:
        amr_lmc_M14_0, amr_lmc_M14_1, amr_lmc_M14_2 = [[], []], [[], []],\
            [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            if l[0] != '--':
                amr_lmc_M14_0[0].append(float(l[0]))
                amr_lmc_M14_0[1].append(float(l[1]))
                amr_lmc_M14_1[0].append(float(l[2]))
                amr_lmc_M14_1[1].append(float(l[3]))
            amr_lmc_M14_2[0].append(float(l[4]))
            amr_lmc_M14_2[1].append(float(l[5]))

    # Put SMC values first and LMC second.
    amr_smc_PG13, amr_lmc_PG13 = zip(*amr_PG13_t[0]), zip(*amr_PG13_t[1])

    amr_lit_smc = [amr_smc_PT98, amr_smc_HZ04, amr_smc_C08, amr_smc_N09,
                   amr_smc_TB09_1, amr_smc_TB09_2, amr_smc_TB09_3,
                   amr_smc_C13_B, amr_smc_C13_C, amr_smc_PG13]
    amr_lit_lmc = [amr_lmc_PT98, amr_lmc_G98, amr_lmc_C08, amr_lmc_HZ09,
                   amr_lmc_R12, amr_lmc_PG13, amr_lmc_M14_0, amr_lmc_M14_1,
                   amr_lmc_M14_2]
    amr_lit = [amr_lit_smc, amr_lit_lmc]

    return amr_lit


def get_liter_data():
    '''
    Read the data file with the literature values for each cluster as a
    dictionary.
    '''

    # Read .ods file with literature data.
    cl_file = ODSBook('lit_OCs_data.ods')
    # Store as dictionary and then as list.
    cl_dict = cl_file.sheets()["S-LMC"]

    return cl_dict


def get_cross_match_asteca(r_path):
    '''
    Read the cross-matched clusters between ASteCA output and several
    databases.
    '''
    # Path to data file.
    in_file = 'mc-catalog/databases/matched_clusters.dat'

    # Read data file
    with open(r_path + in_file, 'r') as f:
        cross_match = [[], [], [], [], [], [], []]

        for line in skip_comments(f):
            lin = line.split()

            if lin[0] == 'P99':
                j = 0
            elif lin[0] == 'P00':
                j = 1
            elif lin[0] == 'H03':
                j = 2
            elif lin[0] == 'R05':
                j = 3
            elif lin[0] == 'C06':
                j = 4
            elif lin[0] == 'G10':
                j = 5
            elif lin[0] == 'P12':
                j = 6

            cross_match[j].append([lin[1]] + [lin[2]] +
                                  [float(_) for _ in lin[3:]])

    for i, lst in enumerate(cross_match):
        cross_match[i] = zip(*lst)

    return cross_match


def get_cross_match_h03_p12(r_path):
    """
    Read cross matched OCs between H03 and P12 databases.
    """
    # Path to data file.
    in_file = 'mc-catalog/databases/matched_H03_P12.dat'

    # Read data file
    with open(r_path + in_file, 'r') as f:
        cross_match = []

        for line in skip_comments(f):
            lin = line.split()
            cross_match.append(map(float, [lin[3], lin[5], lin[7], lin[9]]))

    # 0: log(age)_h03, 1: log(age)_p12, 2: Mass_h03, 3: Mass_p12
    c_m_h03_p12 = zip(*cross_match)

    return c_m_h03_p12
