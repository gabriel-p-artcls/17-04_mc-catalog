
from pyexcel_ods import ODSBook
import os.path


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
    bb_file = 'bb_cat.dat'

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
        amr_lit = [[], []]
        for line in skip_comments(f):
            l = line.split()
            # Read coordinates.
            if l[2] != '--':
                amr_lit[0].append([float(_) for _ in [l[2], l[3]]])
            amr_lit[1].append([float(_) for _ in [l[0], l[1]]])

    return amr_lit


def get_liter_data():
    '''
    Read the data file with the literature values for each cluster as a
    dictionary.
    '''

    # Read .ods file with literature data.
    cl_file = ODSBook('lista_unica_cumulos.ods')
    # Store as dictionary and then as list.
    cl_dict = cl_file.sheets()["S-LMC"]

    return cl_dict


def get_cross_match_data():
    '''
    Read the cross-matched clusters between ASteCA output and several
    databases.
    '''
    # Path to data file.
    in_file = '/../ages_mass_lit/matched_clusters.dat'

    # Read data file
    with open(os.path.dirname(__file__) + in_file, 'r') as f:
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
