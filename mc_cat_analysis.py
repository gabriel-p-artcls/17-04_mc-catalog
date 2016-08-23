
import os
from functions.get_data import get_asteca_data, get_liter_data, \
    get_bica_database, get_cross_match_asteca, get_cross_match_h03_p12,\
    get_amr_lit, get_massclean_data
from functions.amr_kde import get_amr_asteca
from functions.get_params import params
from functions.match_clusters import match_clusters
from functions.check_diffs import check_diffs
from functions.DBs_CMD import get_DBs_ASteCA_CMD_data
from functions.marigo_parsec_isochs import mar_par_data
from functions.make_all_plots import make_as_vs_lit_plot,\
    make_as_vs_lit_mass_plot, make_kde_plots, \
    make_ra_dec_plots, make_lit_ext_plot, make_int_cols_plot, \
    make_concent_plot, make_radius_plot, make_probs_CI_plot, \
    make_cross_match_ip_age, make_cross_match_ip_mass, make_cross_match_if, \
    make_DB_ASteCA_CMDs, make_errors_plots, make_amr_plot,\
    make_cross_match_h03_p12, make_age_mass_corr, make_massclean_z_plot,\
    make_massclean_mass_plot, mar_par_plot


def rpath_fig_folder():
    """
    Obtain root path and create folders where images will be stored.
    """
    # Root path.
    r_path = os.path.realpath(__file__)[:-29]

    # Create folder where the final images will be stored.
    path = 'figures/'
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

    # Create sub-folder for DBs vs ASteCA CMDs.
    path = 'figures/DB_fit/'
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

    return r_path


def get_in_params(r_path):
    """
    Obtain data from ASteCA's output and from the literature.
    Match clusters repeated in both datasets.
    Arrange all parameters in a dictionary.
    """
    # Read data from ASteca output file.
    as_names, as_pars = get_asteca_data()
    print 'ASteCA data read from .dat output file.'

    # Read literature data.
    cl_list = get_liter_data()
    print 'Literature data read from .ods file.'

    # Match clusters.
    names_idx = match_clusters(as_names, cl_list)
    print 'Cluster parameters matched.'

    # Get data parameters arrays.
    in_params = params(r_path, as_names, as_pars, cl_list, names_idx)
    print 'Dictionary of parameters obtained.'

    return in_params


def CMD_DBs_vs_asteca(r_path):
    """
    CMDs of clusters matched between these two databases and
    the values given  by ASteCA.
    """
    print 'Generating CMDs of DBs for matched clusters with ASteCA.'
    for db in ['P99', 'P00', 'C06', 'G10']:
        db_cls = get_DBs_ASteCA_CMD_data(r_path, db, [])
        if db_cls:
            make_DB_ASteCA_CMDs(db, db_cls)


def CMD_outliers(r_path, in_params):
    """
    CMDs of outlier clusters, ie: those with large age differences between
    literature values and the values given  by ASteCA.
    """
    print 'Generating CMDs of outliers.'
    db_cls = get_DBs_ASteCA_CMD_data(r_path, 'outliers', in_params)
    if db_cls:
        make_DB_ASteCA_CMDs('outliers', db_cls)


def CMD_large_mass(r_path, in_params):
    """
    CMDs of clusters with large masses in DBs that are not found by ASteCA.
    """
    print 'Generating CMDs of large-mass OCs.'
    db_cls = get_DBs_ASteCA_CMD_data(r_path, 'largemass', in_params)
    if db_cls:
        make_DB_ASteCA_CMDs('largemass', db_cls)


def CMD_LMC_large_met(r_path, in_params):
    """
    CMDs of clusters with large metallicities and ages in the LMC.
    """
    print 'Generating CMDs of large [Fe/H] and age for the LMC.'
    db_cls = get_DBs_ASteCA_CMD_data(r_path, 'largemet', in_params)
    if db_cls:
        make_DB_ASteCA_CMDs('largemet', db_cls)


def make_plots(r_path, plots, in_params, bica_coords, cross_match,
               cross_match_h03_p12, amr_lit, amr_asteca, massclean_data_pars,
               mar_data, par_data):
    '''
    Make each plot sequentially.
    '''

    if '0' in plots:
        print 'RA vs DEC plots.'
        make_ra_dec_plots(in_params, bica_coords)

    if '1' in plots:
        print '\nErrors plots.'
        make_errors_plots(in_params)

    if '2' in plots:
        print '\nASteCA vs literature plots.'
        make_as_vs_lit_plot(in_params)

    if '3' in plots:
        print '\nASteCA vs literature mass plot.'
        make_as_vs_lit_mass_plot(in_params)

    if '4' in plots:
        print "\nCMDs for outlier clusters."
        CMD_outliers(r_path, in_params)

    if '5' in plots:
        print "\nCMDs for matched clusters between DBs and ASteCA clusters."
        CMD_DBs_vs_asteca(r_path)

    if '6' in plots:
        print '\nCross-matched isochrone fitting clusters.'
        make_cross_match_if(cross_match, in_params)

    if '7' in plots:
        print '\nCross-matched integrated photometry clusters, ages.'
        make_cross_match_ip_age(cross_match)
        print '\nCross-matched integrated photometry clusters, masses.'
        make_cross_match_ip_mass(cross_match)

    if '8' in plots:
        print "\nCross match BA plot for H03 vs P12."
        make_cross_match_h03_p12(cross_match_h03_p12)

    if '9' in plots:
        print "\nAge vs mass delta plots for ASteCA, P12, H03."
        make_age_mass_corr(cross_match, cross_match_h03_p12)

    if '10' in plots:
        print "\nCMDs for large mass clusters."
        CMD_large_mass(r_path, in_params)

    if '11' in plots:
        print '\nKDE maps.'
        make_kde_plots(in_params)

    if '12' in plots:
        print '\nAMR maps.'
        make_amr_plot(in_params, amr_lit, amr_asteca)

    if '13' in plots:
        print '\nMASSCLEAN z plot.'
        make_massclean_z_plot(massclean_data_pars)
        print '\nMASSCLEAN mass plot.'
        make_massclean_mass_plot(massclean_data_pars)

    if '14' in plots:
        print '\nMarigo vs PARSEC plot.'
        mar_par_plot(mar_data, par_data)

    if '15' in plots:
        print '\nASteCA radius (pc) vs parameters plot.'
        make_radius_plot(in_params)

    if '16' in plots:
        print '\nASteCA vs MCEV vs SandF extinction plot.'
        make_lit_ext_plot(in_params)

    if '17' in plots:
        print '\nIntegrated colors plot.'
        make_int_cols_plot(in_params)

    if '18' in plots:
        print '\nConcentration parameter plot.'
        make_concent_plot(in_params)

    if '19' in plots:
        print '\nASteCA probabilities versus CI.'
        make_probs_CI_plot(in_params)

    if '20' in plots:
        print "\nCMDs for large [Fe/H] LMC clusters."
        CMD_LMC_large_met(r_path, in_params)


def main():
    '''
    Call each function.
    '''
    # Root path.
    r_path = rpath_fig_folder()

    # Obtain data from ASteCA's output and from the literature.
    in_params = get_in_params(r_path)

    # Check for differences in ASteCA vs Lit values.
    check_diffs(in_params)

    # Define which plots to produce.
    plots = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
             '12', '13', '14', '15', '16', '17', '18', '19', '20']

    bica_coords, cross_match, cross_match_h03_p12, amr_lit, amr_asteca, \
        massclean_data_pars, mar_data, par_data = [], [], [], [], [], [], [],\
        []
    # Only obtain data if the plot is being generated.
    if '0' in plots:
        # Read Bica et al. (2008) database.
        bica_coords = get_bica_database()
        print 'Bica et al. (2008) data read.'
    if any(i in ['6', '7', '9'] for i in plots):
        # Read cross-matched ASteCA clusters.
        cross_match = get_cross_match_asteca(r_path)
        print 'Cross-matched ASteCA data read.'
    if any(i in ['8', '9'] for i in plots):
        # Read cross-matched H03,P12 clusters.
        cross_match_h03_p12 = get_cross_match_h03_p12(r_path)
        print 'Cross-matched H03,P12 data read.'
    if '12' in plots:
        # Read AMR data from other articles.
        amr_lit = get_amr_lit()
        print 'AMR data from literature read.'
        amr_asteca = get_amr_asteca(in_params)
        print 'ASteCA AMR for both MCs obtained.'
    if '13' in plots:
        massclean_data_pars = get_massclean_data()
        print 'MASSCLEAN data read.'
    if '14' in plots:
        mar_data, par_data = mar_par_data()

    # Make final plots.
    print 'Plotting...\n'
    make_plots(r_path, plots, in_params, bica_coords, cross_match,
               cross_match_h03_p12, amr_lit, amr_asteca, massclean_data_pars,
               mar_data, par_data)

    print '\nEnd.'


if __name__ == "__main__":
    main()
