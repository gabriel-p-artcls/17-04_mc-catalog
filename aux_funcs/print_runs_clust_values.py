

def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def get_asteca_data(i, run, clust):
    '''
    Read the ASteCA output data file 'asteca_output.dat' and store each
    data column for each cluster.
    '''

    # Path to data file.
    out_file = '../runs/' + run + '_run/' + 'asteca_output_' + run + '.dat'

    # Read data file
    with open(out_file, 'r') as f:

        for line in skip_comments(f):
            name = line.split()[0]
            if name == clust:
                print i, name, line.split()[24], line.split()[26], \
                    line.split()[22], line.split()[23]

    return


def main():
    '''
    Read the cross-matched clusters between ASteCA output and several
    databases.
    '''

    runs = ['1st', '2nd', '3rd', '4th', '5th', '6th', '7th', '8th', '9th',
            '10th', '11th', '18th', '19th']

    # 'BSDL654' 'SL218' 'H88-131' 'KMHK975' 'BSDL631' 'L35' 'SL579'
    clust = 'SL579'

    for i, run in enumerate(runs):
        get_asteca_data((i + 1), run, clust)


if __name__ == "__main__":
    main()
