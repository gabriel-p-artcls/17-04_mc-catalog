
import functions.CMD_obs_vs_asteca as cmd


def get_cross_match_OCs(db):
    """
    Read OCs that could be cross-matched from file.
    """
    if db == 'outliers':
        mc_cls = [['KMHK975', 'SL579', 'BSDL631', 'KMHK979', 'H88-316',
                   'L39', 'L35', 'H86-188', 'B134', 'K47']]
        # Isochrones used in the analysis of the above clusters.
        isochs = [['M08', 'M08', 'M08', 'M08', 'M08', 'M08', 'G02', 'M08',
                   'M08', 'M08']]
    elif db == 'largemass':
        mc_cls = [['NGC419', 'NGC1917', 'NGC1751', 'L27', 'SL244']]
        # Use this list to pass galaxies info.
        isochs = [['SMC', 'LMC', 'LMC', 'SMC', 'LMC']]
    else:
        isochs = []  # dummy
        # Path to data file.
        out_file = 'databases/matched_clusters.dat'
        # Read data file
        with open(out_file) as f:
            cl_names = []
            for line in f:
                if not line.startswith('#'):
                    l = line.split()
                    if l[0] == db:
                        cl_names.append(l[2])
        # Split list into chunks of N clusters.
        N = 10
        mc_cls = [cl_names[x:x+N] for x in xrange(0, len(cl_names), N)]

    return mc_cls, isochs


def get_lit_params(r_path, cl, db, in_params, isochs, i, j):
    """
    Obtain fundamental parameters for 'cl' cluster given by the literature.
    Also return the isochrone used in the literature to analyse it.
    """

    # Obtain age and extinction. In the case of DBs, get from
    # the 'matched_clusters.dat' file.
    db_a, db_e, gal = cmd.get_DB_age_ext(r_path, cl, db, in_params)

    # Set metallicity, distance modulus values, and isochrone set used.
    if db == 'P99':
        [db_z, db_d] = [0.004, 18.65]
        # Uses Padova isochrones (Bertelli et al. 1994)
        isoc = 'G02'
    elif db == 'P00':
        [db_z, db_d] = [0.008, 18.24]
        # Uses Padova isochrones (Bertelli et al. 1994)
        isoc = 'G02'
    elif db == 'C06':
        [db_z, db_d] = [0.008, 18.9]
        isoc = 'G02'
    elif db == 'G10':
        [db_z, db_d] = [0.004, 18.9] if gal == 'SMC' else [0.008, 18.5]
        # Uses Padova isochrones (Girardi et al. 1995)
        isoc = 'G02'
    elif db == 'outliers':
        [db_z, db_d] = [0.004, 18.9] if gal == 'SMC' else [0.008, 18.5]
        isoc = isochs[i][j]

    # Obtain DB isochrone.
    lit_isoch = cmd.get_isoch(r_path, db, isoc, db_z, db_a, db_e, db_d)

    return gal, db_z, db_a, db_e, db_d, lit_isoch


def get_CMD_data(r_path, db, in_params, mc_cls, isochs):
    '''
    Obtain cross-matched OCs data, used to generate figures containing CMDs
    of databases vs ASteCA for the clusters matched in the selected database.
    '''

    db_cls = [[] for _ in mc_cls]
    for i, cl_lst in enumerate(mc_cls):
        for j, cl in enumerate(cl_lst):

            # Find photometric file for cluster.
            phot_data = cmd.find_phot_file(r_path, cl)

            # Obtain ASteCA parameters.
            as_z, as_z_str, as_a, as_e, as_d, as_m = cmd.get_asteca_params(cl)
            # Obtain ASteCA isochrone.
            asteca_isoch = cmd.get_isoch(r_path, 'AS', '', as_z_str, as_a,
                                         as_e, as_d)

            # Literature values.
            if db != 'largemass':
                gal, db_z, db_a, db_e, db_d, lit_isoch = get_lit_params(
                    r_path, cl, db, in_params, isochs, i, j)
            else:
                # For 'largemass' OCs, this list passes galaxies info instead
                # of isochrones.
                gal = isochs[i][j]
                # Take the rest of the info from ASteCA values.
                db_z, db_a, db_e, db_d, lit_isoch = as_z, as_a, as_e, as_d,\
                    asteca_isoch

            # Fetch which run holds this cluster's membership data.
            run = cmd.get_cl_run(cl)
            # Fetch what 'input_XX' folder in the above run contains the
            # membership file.
            inpt = cmd.get_input_folder(r_path, cl, run)

            # Membership data for cluster.
            cl_reg_fit, cl_reg_no_fit, synth_stars = cmd.get_memb_data(
                r_path, run, inpt, cl)

            # Obtain CMD limits for cluster.
            x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = cmd.diag_limits(
                phot_data)

            db_cls[i].append([x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, cl,
                              db, gal, cl_reg_fit, cl_reg_no_fit, synth_stars,
                              lit_isoch, asteca_isoch, db_z, db_a, db_e, db_d,
                              as_z, as_a, as_e, as_d, as_m])
            print '{} {} data obtained'.format(db, cl)

    return db_cls


def get_DBs_ASteCA_CMD_data(r_path, db, in_params):
    """
    Gather information to plot CMDs of cross-matched OCs in databases, and
    age outliers.
    """
    # Read OCs names and set of isochrones used.
    mc_cls, isochs = get_cross_match_OCs(db)
    # Obtain data for each OC.
    db_cls = get_CMD_data(r_path, db, in_params, mc_cls, isochs)

    return db_cls
