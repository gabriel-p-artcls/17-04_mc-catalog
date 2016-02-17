
import os
from functions.get_data import get_asteca_data, get_liter_data, \
    get_bica_database, get_cross_match_data, get_amr_lit
from functions.get_params import params
from functions.match_clusters import match_clusters
from functions.galax_struct_dist import gsd
from functions.check_diffs import check_diffs
from functions.DBs_CMD import get_DBs_ASteCA_CMD_data
from functions.make_all_plots import make_as_vs_lit_plot, make_kde_plots, \
    make_ra_dec_plots, make_lit_ext_plot, make_int_cols_plot, \
    make_concent_plot, make_radius_plot, make_probs_CI_plot, \
    make_dist_2_cents, make_cross_match, make_cross_match_age_ext, \
    make_DB_ASteCA_CMDs, make_errors_plots, make_amr_plot, make_angles_plot


def make_plots(in_params, bica_coords, cross_match, amr_lit, gal_str_pars):
    '''
    Make each plot sequentially.
    '''

    # for j, gal in enumerate(['SMC', 'LMC']):
    #     make_kde_plots(gal, j, in_params)
    #     print '{} KDE maps done.'.format(gal)

    # make_as_vs_lit_plot(in_params)
    # print 'ASteCA vs literature plots done.'

    # make_errors_plots(in_params)
    # print 'Errors plots done.'

    # make_ra_dec_plots(in_params, bica_coords)
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

    # make_cross_match_age_ext(cross_match, in_params)
    # print 'Age and extinction diffs for cross-matched clusters done.'

    # make_amr_plot(in_params, amr_lit)
    # print 'AMR maps done.'

    make_angles_plot(gal_str_pars)
    print 'Inclination vs position angles plot done.'


def CMD_DBs_vs_asteca(r_path):
    """
    CMDs of clusters matched between these two databases and
    the values given  by ASteCA.
    """
    print 'Generating CMDs of DBs for matched clusters with ASteCA.'
    for db in ['G10', 'C06']:
        db_cls = get_DBs_ASteCA_CMD_data(r_path, db)
        make_DB_ASteCA_CMDs(db, db_cls)
        print ("CMDs for matched clusters between {} DB and"
               " ASteCA clusters done.".format(db))


def main():
    '''
    Call each function.
    '''
    # Root path.
    r_path = os.path.realpath(__file__)[:-29]

    # Generate CMDs of DBs vs ASteCA.
    # CMD_DBs_vs_asteca(r_path)

    # Read data from ASteca output file.
    as_names, as_pars = get_asteca_data()
    print 'ASteCA data read from .dat output file.'

    # Read literature data.
    cl_dict = get_liter_data()
    print 'Literature data read from .ods file.'

    # Match clusters.
    names_idx = match_clusters(as_names, cl_dict)
    print 'Cluster parameters matched.'

    # Get data parameters arrays.
    in_params = params(r_path, as_names, as_pars, cl_dict, names_idx)
    print 'Dictionary of parameters obtained.'

    # Check for differences in ASteCA vs Lit values.
    check_diffs(in_params)

    # Obtain galactic structure (inclination + position angles) for MCs
    gal_str_pars = gsd(in_params)
    print 'Inclination and position angles for MCs obtained.'

    # Read cross-matched clusters.
    cross_match = get_cross_match_data()
    print 'Cross-matched data read.'

    # Read Bica et al. (2008) database.
    bica_coords = get_bica_database()

    # Read AMR data from other articles.
    amr_lit = get_amr_lit()

    # Make final plots.
    print 'Plotting...\n'
    make_plots(in_params, bica_coords, cross_match, amr_lit, gal_str_pars)

    print '\nEnd.'


if __name__ == "__main__":
    main()
