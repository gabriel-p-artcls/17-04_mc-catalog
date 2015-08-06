
from functions.get_data import get_asteca_data, get_liter_data, \
    get_cross_match_data
from functions.get_params import params
from functions.make_all_plots import make_as_vs_lit_plot, make_kde_plots, \
    make_ra_dec_plots, make_lit_ext_plot, make_int_cols_plot, \
    make_concent_plot, make_radius_plot, make_probs_CI_plot, \
    make_dist_2_cents, make_cross_match, make_cross_match_age_ext


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
    gal_names, zarr, aarr, earr, darr, rarr, marr, dist_cent, ra, dec = \
        [in_params[_] for _ in ['gal_names', 'zarr', 'aarr', 'earr', 'darr',
                                'rarr', 'marr', 'dist_cent', 'ra', 'dec']]

    gal = ['SMC', 'LMC']

    # For SMC and LMC.
    for j in [0, 1]:

        # For each cluster.
        cl_count = 0
        for i, name in enumerate(gal_names[j]):
            flag_cl = False

            # For each parameter.
            for k, par in enumerate([aarr]):
                diff = par[j][0][i] - par[j][1][i]
                if par[j][1][i] > -99.:

                    # # Age diffs.
                    # print '{} {}, {:.2f} vs {:.2f} , {:.2f}'.format(
                    #     gal[j], name, par[j][0][i], par[j][1][i], diff)

                    # Age.
                    if abs(diff) > 0.4:
                        flag_cl = True
                        # AsteCA vs Literature Log(age) difference.
                        rel_diff = par[j][0][i] - par[j][1][i]
                        print '{} {}, {:.2f} vs {:.2f} , {:.2f}'.format(
                            gal[j], name, par[j][0][i], par[j][1][i], rel_diff)

            #         # Extinction.
            #         if k == 2 and diff > pars_diff[2]:
            #             flag_cl = True
            #             print '{} {} {}, {:.2f} vs {:.2f}'.format(gal[j],
            #                 name, p_n[k], par[j][0][i], par[j][1][i])

            #         # Distance.
            #         if k == 3 and diff > pars_diff[3]:
            #             flag_cl = True
            #             print '{} {} {}, {:.2f} vs {:.2f}'.format(gal[j],
            #                 name, p_n[k], par[j][0][i], par[j][1][i])

            #         # Radius.
            #         if k == 4 and diff > pars_diff[4]:
            #             flag_cl = True
            #             print '{} {} {}, {} vs {}'.format(gal[j], name,
            #                 p_n[k], par[j][0][i], par[j][1][i])

            # # Mass.
            # if marr[j][0][i] > 0.:
            #     flag_cl = True
            #     print '{} {} {}'.format(gal[j], name, marr[j][0][i])

            # # Distance to center.
            # print '{} {} {:.2f} {} {} {}'.format(
            #     gal[j], name, dist_cent[j][i] / 1000., darr[j][0][i],
            #      ra[j][i], dec[j][i])

            if flag_cl:
                cl_count += 1

        print '\n* {}, Clusters with \delta log(age)>0.4: {}\n'.format(
            gal[j], cl_count)


def make_plots(in_params, cross_match):
    '''
    Make each plot sequentially.
    '''

    # for j, gal in enumerate(['SMC', 'LMC']):
    #     make_as_vs_lit_plot(gal, j, in_params)
    #     print '{} ASteCA vs literature plots done.'.format(gal)

    #     make_kde_plots(gal, j, in_params)
    #     print '{} KDE maps done.'.format(gal)

    # make_ra_dec_plots(in_params)
    # print 'RA vs DEC plots done.'

    # make_lit_ext_plot(in_params)
    # print 'ASteca vs MCEV vs SandF extinction plot done.'

    # make_int_cols_plot(in_params)
    # print 'Integrated colors plot done.'

    # make_concent_plot(in_params)
    # print 'Concentration parameter plot done.'

    # make_radius_plot(in_params)
    # print 'ASteCA radius (pc) vs parameters plot done.'

    # make_probs_CI_plot(in_params)
    # print 'ASteCA probabilities versus CI done.'

    # make_dist_2_cents(in_params)
    # print 'Distances to center of MC done.'

    # make_cross_match(cross_match)
    # print 'Cross-matched clusters done.'

    make_cross_match_age_ext(cross_match, in_params)
    print 'Age and extinction diffs for cross-matched clusters done.'


def main():
    '''
    Call each function.
    '''

    # import sys
    # print(sys.version)
    # import pip
    # pip.get_installed_distributions(local_only=False)
    # raw_input

    # Read data from ASteca output file.
    as_names, as_pars = get_asteca_data()
    print 'ASteCA data read from output file.'

    # Read literature data.
    cl_dict = get_liter_data()
    print 'Literature data read from .ods file.'

    # Read cross-matched clusters.
    cross_match = get_cross_match_data()
    print 'Cross-matched data read.'

    # Match clusters.
    names_idx = match_clusters(as_names, cl_dict)
    print 'Cluster parameters matched.'

    # Get data parameters arrays.
    in_params = params(as_names, as_pars, cl_dict, names_idx)
    print 'Dictionary of parameters obtained.\n'

    # Check for differences in ASteCA vs Lit values.
    check_diffs(in_params)

    # Make final plots.
    make_plots(in_params, cross_match)

    print '\nEnd.'


if __name__ == "__main__":
    main()
