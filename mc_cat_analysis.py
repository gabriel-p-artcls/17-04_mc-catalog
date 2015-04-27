
#from astropy import units as u
#from astropy.coordinates import SkyCoord

from functions.get_data import get_asteca_data, get_liter_data
from functions.get_params import params
from functions.make_all_plots import make_as_vs_lit_plot, make_kde_plots, \
make_ra_dec_plots, make_lit_ext_plot


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


def check_diffs(in_params):
    '''
    check differences between ASteCA values and literature values for given
    parameters.
    '''
    gal_names, zarr, aarr, earr, darr, rarr = [in_params[_] for _ in
    ['gal_names', 'zarr', 'aarr', 'earr', 'darr', 'rarr']]

    gal = ['SMC', 'LMC']
    p_n = ['metal', 'age', 'ext', 'dist', 'rad']
    # Max diff acceptable for each parameter.
    pars_diff = [0.1, 0.2, 0.1, 0.1, 30]

    # For SMC and LMC.
    for j in [0, 1]:

        # For each cluster.
        cl_count = 0
        for i, name in enumerate(gal_names[j]):
            flag_cl = False

            # For each parameter.
            for k, par in enumerate([zarr, aarr, earr, darr, rarr]):
                diff = abs(par[j][0][i] - par[j][1][i])

                if par[j][1][i] > -99.:

                    # Metallicity.
                    if k == 0 and diff > pars_diff[0]:
                        flag_cl = True
                        print '{} {} {}, {:.4f} vs {:.4f}'.format(gal[j],
                            name, p_n[k], par[j][0][i], par[j][1][i])

                    # Age.
                    if k == 1 and diff > pars_diff[1]:
                        flag_cl = True
                        print '{} {} {}, {:.2f} vs {:.2f}'.format(gal[j],
                            name, p_n[k], par[j][0][i], par[j][1][i])

                    # Extinction.
                    if k == 2 and diff > pars_diff[2]:
                        flag_cl = True
                        print '{} {} {}, {:.2f} vs {:.2f}'.format(gal[j],
                            name, p_n[k], par[j][0][i], par[j][1][i])

                    # Distance.
                    if k == 3 and diff > pars_diff[3]:
                        flag_cl = True
                        print '{} {} {}, {:.2f} vs {:.2f}'.format(gal[j],
                            name, p_n[k], par[j][0][i], par[j][1][i])

                    # Radius.
                    if k == 4 and diff > pars_diff[4]:
                        flag_cl = True
                        print '{} {} {}, {} vs {}'.format(gal[j], name,
                            p_n[k], par[j][0][i], par[j][1][i])

            if flag_cl:
                cl_count += 1

        print '\n{}, total num of clusts with diffs in params: {}'.format(
            gal[j], cl_count)


def make_plots(i, in_params):
    '''
    Make each plot sequentially.
    '''

    print 'Making SMC ASteCA vs literature plots.'
    make_as_vs_lit_plot(i, 'SMC', 0, in_params)

    print 'Making LMC ASteCA vs literature plots.'
    make_as_vs_lit_plot(i, 'LMC', 1, in_params)

    print 'Making SMC KDE maps.'
    make_kde_plots(i, 'SMC', 0, in_params)

    print 'Making LMC KDE maps.'
    make_kde_plots(i, 'LMC', 1, in_params)

    print 'Making RA vs DEC plots.'
    make_ra_dec_plots(i, in_params)

    print 'Make ASteca vs MCEV vs SandF extinction plot.'
    make_lit_ext_plot(in_params)


def main():
    '''
    Call each function.
    '''

    # Read data from ASteca output file.
    # 0 (25/04/15): SMC clusters, restricted range in extinction and dist mod,
    # semi center with 100px search area, auto rads, DA + bayesian blocks,
    # knuth binning.
    i = '0'

    as_names, as_pars = get_asteca_data(i)
    print 'ASteCA data read from output file.'

    # Read literature data.
    cl_dict = get_liter_data()
    print 'Literature data read from .ods file.'

    # Match clusters.
    names_idx = match_clusters(as_names, cl_dict)
    print 'Cluster parameters matched.'

    # Get data parameters arrays.
    in_params = params(as_names, as_pars, cl_dict, names_idx)
    print 'Dictionary of parameters obtained.\n'

    # Check for differences in ASteCA vs Lit values.
    check_diffs(in_params)

    # Make final plots.
    #make_plots(i, in_params)

    print '\nEnd.'


if __name__ == "__main__":
    main()