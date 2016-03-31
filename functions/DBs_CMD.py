
import os
import functions.CMD_obs_vs_asteca as cmd


def get_DBs_ASteCA_CMD_data(r_path, db, in_params):
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

    elif db == 'outliers':
        mc_cls = [['KMHK975', 'SL579', 'BSDL631', 'KMHK979', 'H88-316',
                   'L39', 'L35', 'H86-188', 'B134', 'K47']]
        # Isochrones used in the analysis of the above clusters.
        isochs = [['M08', 'M08', 'M08', 'M08', 'M08', 'M08', 'G02', 'M08',
                   'M08', 'M08']]

    db_cls = [[] for _ in mc_cls]
    for i, cl_lst in enumerate(mc_cls):
        for j, cl in enumerate(cl_lst):

            # Find photometric file for cluster.
            phot_data = cmd.find_phot_file(r_path, cl)

            # Obtain age and extinction from 'matched_clusters.dat' file.
            db_a, db_e, gal = cmd.get_DB_age_ext(r_path, cl, db, in_params)
            # Set metallicity and distance modulus values.
            if gal == 'SMC':
                if db in ['G10', 'outliers']:
                    db_z, db_d = 0.004, 18.9
                elif db == 'C06':
                    db_z, db_d = 0.008, 18.9
            elif gal == 'LMC' and db in ['G10', 'outliers']:
                db_z, db_d = 0.008, 18.5

            if db == 'G10':
                [db_z, db_d] = [0.004, 18.9] if gal == 'SMC' else [0.008, 18.5]
                # Actually uses Padova isochrones (Girardi et al. 1995)
                isoc = 'G02'
            elif db == 'C06':
                [db_z, db_d] = [0.008, 18.9]
                isoc = 'G02'
            elif db == 'outliers':
                [db_z, db_d] = [0.004, 18.9] if gal == 'SMC' else [0.008, 18.5]
                isoc = isochs[i][j]

            # Obtain DB isochrone.
            lit_isoch = cmd.get_isoch(r_path, db, isoc, db_z, db_a, db_e, db_d)

            # Obtain ASteCA parameters.
            as_z, as_z_str, as_a, as_e, as_d = cmd.get_asteca_params(cl)
            # Obtain ASteCA isochrone.
            asteca_isoch = cmd.get_isoch(r_path, 'AS', '', as_z_str, as_a,
                                         as_e, as_d)

            # Fetch which run holds this cluster's membership data.
            run = cmd.get_cl_run(cl)
            # Fetch what 'input_XX' folder in the above run contains the
            # membership file.
            inpt = cmd.get_input_folder(r_path, cl, run)
            # Membership data for cluster.
            cl_reg_fit, cl_reg_no_fit = cmd.get_memb_data(r_path, run, inpt,
                                                          cl)

            # Obtain CMD limits for cluster.
            x_max_cmd, x_min_cmd, y_min_cmd, y_max_cmd = cmd.diag_limits(
                phot_data)

            db_cls[i].append([x_min_cmd, x_max_cmd, y_min_cmd, y_max_cmd, cl,
                              db, gal, cl_reg_fit, cl_reg_no_fit, lit_isoch,
                              asteca_isoch, db_z, db_a, db_e, db_d, as_z, as_a,
                              as_e, as_d])
            print '{} {} data obtained'.format(db, cl)

    return db_cls
