

def d_search(dat_lst, cl_name, name_idx):
    '''
    Search the list of lists obtained when reading the .ods file, for the
    index of the list that contains a given cluster.
    '''
    for i, line in enumerate(dat_lst):
        if cl_name == line[name_idx]:
            return i
    return None


def match_clusters(as_names, cl_dict):
    '''
    Return the index pointing to each cluster in the .ods list, starting
    from the cluster's name taken from the ASteCA list.
    '''

    # Column number for the cluster's name in the .ods file.
    name_idx = cl_dict[0].index(u'Name')

    names_idx = []
    for cl_name in as_names:
        cl_i = d_search(cl_dict, cl_name, name_idx)
        if cl_i is None:
            print 'WARNING: {} not found in ods file.'.format(cl_name)
        else:
            names_idx.append(cl_i)

    return names_idx
