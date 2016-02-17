
import os
import functions.CMD_obs_vs_asteca as cmd


def get_DBs_ASteCA_CMD_data(r_path, db):
    '''
    Generate figures containing CMDs of databases vs ASteCA for the clusters
    matched in the selected database.
    '''
    if db == 'G10':
        mc_cls = [['B112', 'SL218', 'KMHK979', 'KMHK229', 'HW22', 'HW42',
                   'HW63', 'KMHK378', 'SL446A', 'L91'],
                  ['SL674', 'SL290', 'HW40', 'HW31', 'BSDL268', 'HW41',
                   'SL162', 'SL230', 'SL555', 'SL132'],
                  ['HS264', 'BRHT4B', 'SL96', 'L63', 'HS38', 'NGC1839',
                   'NGC294', 'B34', 'NGC1793', 'L72'],
                  ['NGC2093', 'KMHK112', 'BS265', 'SL678', 'SL35', 'B39',
                   'L50', 'L30', 'SL397', 'NGC1863'],
                  ['BRHT45A', 'HW55', 'NGC1838', 'KMHK1055', 'SL444', 'L62',
                  'SL505', 'L34', 'H88-320', 'HS412'],
                  ['LW54', 'L58', 'L49', 'SL510', 'SL551', 'BSDL631', 'L45',
                  'H88-316', 'BS35', 'L35'],
                  ['SL579']]
    elif db == 'C06':
        mc_cls = [['B47', 'H86-70', 'L63', 'L62', 'B39', 'BS121', 'BS88',
                   'NGC294', 'L19', 'L34'],
                  ['L30', 'B34', 'L72', 'NGC419', 'BS35', 'L35']]

    # Create folder where the final images will be stored.
    path = 'figures/DB_fit/'
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

    db_cls = [[] for _ in mc_cls]
    for i, cl_lst in enumerate(mc_cls):
        for cl in cl_lst:

            # Obtain CMD limits for cluster.
            x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = cmd.diag_limits(
                r_path, cl)

            # Obtain age and extinction from 'matched_clusters.dat' file.
            db_a, db_e, gal = cmd.get_DB_age_ext(r_path, cl, db)
            if gal == 'SMC':
                if db == 'G10':
                    db_z, db_d = 0.004, 18.9
                elif db == 'C06':
                    db_z, db_d = 0.008, 18.9
            elif gal == 'LMC' and db == 'G10':
                db_z, db_d = 0.008, 18.5

            # Fetch which run holds this cluster's membership data.
            run = cmd.get_cl_run(cl)
            # Fetch what 'input_XX' folder in the above run contains the
            # membership file.
            inpt = cmd.get_input_folder(r_path, cl, run)
            # Membership data for cluster.
            cl_reg_fit, cl_reg_no_fit = cmd.get_memb_data(r_path, run, inpt,
                                                          cl)

            # Obtain DB isochrone.
            lit_isoch = cmd.get_isoch(r_path, 'DB', db_z, db_a, db_e, db_d)

            # Obtain ASteCA parameters.
            as_z, as_z_str, as_a, as_e, as_d = cmd.get_asteca_params(cl)
            # Obtain ASteCA isochrone.
            asteca_isoch = cmd.get_isoch(r_path, 'AS', as_z_str, as_a, as_e,
                                         as_d)

            db_cls[i].append([x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, cl,
                              db, gal, cl_reg_fit, cl_reg_no_fit, lit_isoch,
                              asteca_isoch, db_z, db_a, db_e, db_d, as_z, as_a,
                              as_e, as_d])
            print '{} {} data obtained'.format(db, cl)

    return db_cls
