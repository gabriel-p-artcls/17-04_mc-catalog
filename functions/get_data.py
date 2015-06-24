
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


def get_liter_data():
    '''
    Read the ASteCA output data file 'asteca_output.dat' and store each
    data column for each cluster.
    '''

    # Read .ods file with literature data.
    cl_file = ODSBook('lista_unica_cumulos.ods')
    # Store as dictionary and then as list.
    cl_dict = cl_file.sheets()["S-LMC"]

    return cl_dict
