

def read_data():
    '''
    Read data from 'asteca_output.dat'.
    '''

    in_data = []
    with open('../asteca_output.dat', "r") as f:
        for line in f:
            if line.strip().startswith('#'):
                in_data.append(line)
            else:
                in_data.append(line.split())

    return in_data


def read_p_values():
    '''
    Read p_value probabilities from 'asteca_output.dat' file with the p-values.
    '''

    in_p_values = []
    with open('asteca_output_pvals.dat', "r") as f:
        for line in f:
            if line.strip().startswith('#'):
                pass
            else:
                # Store cluster's name and prob value.
                in_p_values.append([line.split()[0], line.split()[18]])

    return in_p_values


def add_p_values(in_data, in_p_values):
    '''
    Match each cluster and add its p_value.
    '''

    in_data_out = []
    for line in in_data:
        # If line is cmmted out, store as whole.
        if isinstance(line, basestring):
            in_data_out.append(line)
        else:
            cl_name = line[0]
            for clust in in_p_values:
                if cl_name == clust[0]:
                    line_p_val = line
                    line_p_val[18] = clust[1]
                    in_data_out.append(line_p_val)

    return in_data_out


def out_data(in_data_out):
    '''
    Write data from 'asteca_output.dat' with p_values assigned to final merged
    file.
    '''

    # Write values to file.
    with open('asteca_output_final.dat', "a") as f_out:
        for line in in_data_out:
            if isinstance(line, basestring):
                f_out.write(line)
            else:
                f_out.write('''{:<16} {:>8} {:>8} {:>8} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>7} {:>10} \
{:>10} {:>10} {:>9} {:>7} {:>8} {:>8} \
{:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}\
{:>8} {:>2} {:>3} {:>2} {:>2} {:>2} {:>2} {:>2} {:>2} \
{:>2} {:>2} {:>3} {:>3}'''.format(*line))
                f_out.write('\n')

    return


def main():
    '''
    Add KDE p-value probabilities for each cluster since both the 1st and 2nd
    run were made with this function turned off.
    '''

    # Read data from ASteCA output file.
    in_data = read_data()

    # Read p_values from 2nd ASteCA output file.
    in_p_values = read_p_values()

    # Add p_values matching each cluster.
    in_data_out = add_p_values(in_data, in_p_values)

    # Write to final output file.
    out_data(in_data_out)

    print '\nEnd.'


if __name__ == "__main__":
    main()
