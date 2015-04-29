import numpy as np


def float_str(val):
    '''
    Return float depending on whether the value is a float or a string.
    '''
    try:
        val_f = float(val)
    except:
        val_f = -99.9

    return val_f


def rad_in_pc(float_lst):
    '''
    Convert radius from pixel to arcsec to parsecs.
    '''
    # Unpack
    r_px, scale, dis_mod, e_bv = float_lst
    # Convert to parsecs.
    rad_arcsec = r_px * scale
    Av = 3.1 * e_bv
    d_pc = 10 ** (0.2 * (dis_mod + 5 - Av))
    r_pc = d_pc * np.tan(np.deg2rad(rad_arcsec / 3600.))

    return r_pc


def correct_int_col_extin(int_col, extinc):
    '''
    Corret integrated color for the extinction.
    '''
    E_CT1_E_BV = 1.97
    int_col_cor = int_col - E_CT1_E_BV * float_str(extinc)

    return int_col_cor


def params(as_names, as_pars, cl_dict, names_idx):
    '''
    Return ASteCA output and literature parameters values.
    '''

    # Indexes of columns in ASteCA output file.
    a_zi, a_zei, a_ai, a_aei, a_ei, a_eei, a_di, a_dei, a_mi, a_mei, a_rad, \
    a_int_c, a_nmemb = 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 4, 18, 13

    # Indexes of columns in .ods literature file.
    ra_i, dec_i, gal_i, l_zi, l_zei, l_ai, l_aei, l_ei, l_eei, l_di, l_dei, \
    l_rad, l_scale, l_e_sandf, l_e_e_sandf, l_e_mcev, l_e_e_mcev = \
    cl_dict[0].index(u'ra_deg'), cl_dict[0].index(u'dec_deg'), \
    cl_dict[0].index(u'Galaxia'), cl_dict[0].index(u'[Fe/H] (dex)'), \
    cl_dict[0].index(u'e_Fe/H'), cl_dict[0].index(u'log(age)'), \
    cl_dict[0].index(u'e_log(age)'), cl_dict[0].index(u'E(B-V) (lit)'), \
    cl_dict[0].index(u'e_E(B-V)'), cl_dict[0].index(u'(m-M)o (mag)'), \
    cl_dict[0].index(u'e_(m-M)o'), cl_dict[0].index(u'rad (eye)'), \
    cl_dict[0].index(u'arcsec/pixel'), cl_dict[0].index(u'E_B_V_SandF'), \
    cl_dict[0].index(u'stdev_E_B_V_SandF'), cl_dict[0].index(u'E_BV_max'), \
    cl_dict[0].index(u'E_BV_std_dev')

    # Initialize empty lists. The first sub-list in the paramaters list
    # corresponds to clusters in the SMC and the second to those in the LMC.
    gal_names = [[], []]
    int_colors = [[], []]
    n_memb = [[], []]
    ra, dec = [[], []], [[], []]
    ext_sf, ext_mcev = [[[], []], [[], []]], [[[], []], [[], []]]
    # Cluster radius in parsecs.
    rad_pc = [[], []]
    zarr, zsigma, aarr, asigma, earr, esigma, darr, dsigma, marr, marr, \
    msigma, rarr = ([[[], []], [[], []]] for i in range(12))

    for i, as_p in enumerate(as_pars):

        # j = 0 for SMC clusters and 1 for LMC.
        j = 0 if cl_dict[names_idx[i]][gal_i] == 'SMC' else 1

        # Store names for each galaxy.
        gal_names[j].append(as_names[i])
        # Store coordinates.
        ra[j].append(cl_dict[names_idx[i]][ra_i])
        dec[j].append(cl_dict[names_idx[i]][dec_i])
        # Integrated color.
        int_col_no_corr = float_str(as_pars[i][a_int_c])
        # Correct for extinction.
        int_col_corr = correct_int_col_extin(int_col_no_corr, as_pars[i][a_ei])
        # Store extinction corrected integrated color.
        int_colors[j].append(int_col_corr)
        # Approx number of members.
        n_memb[j].append(float_str(as_pars[i][a_nmemb]))
        # Store literature E(B-V) values: Schlafly & Finkbeiner (SandF) and
        # MCEV.
        ext_sf[j][0].append(cl_dict[names_idx[i]][l_e_sandf])
        ext_sf[j][1].append(cl_dict[names_idx[i]][l_e_e_sandf])
        ext_mcev[j][0].append(cl_dict[names_idx[i]][l_e_mcev])
        ext_mcev[j][1].append(cl_dict[names_idx[i]][l_e_e_mcev])
        # Store radius value in parsecs.
        float_lst = []
        for el in [as_p[a_rad], cl_dict[names_idx[i]][l_scale], as_p[a_di],
            as_p[a_ei]]:
                # Store in list as floats.
            float_lst.append(float_str(el))
        r_pc = rad_in_pc(float_lst)
        #print as_names[i], as_p[a_rad], cl_dict[names_idx[i]][l_scale], r_pc
        rad_pc[j].append(r_pc)

        # Organize param values, ASteCA first, lit second.
        met = [as_pars[i][a_zi], cl_dict[names_idx[i]][l_zi]]
        smet = [as_pars[i][a_zei], cl_dict[names_idx[i]][l_zei]]
        age = [as_pars[i][a_ai], cl_dict[names_idx[i]][l_ai]]
        sage = [as_pars[i][a_aei], cl_dict[names_idx[i]][l_aei]]
        ext = [as_pars[i][a_ei], cl_dict[names_idx[i]][l_ei]]
        sext = [as_pars[i][a_eei], cl_dict[names_idx[i]][l_eei]]
        dis = [as_pars[i][a_di], cl_dict[names_idx[i]][l_di]]
        sdis = [as_pars[i][a_dei], cl_dict[names_idx[i]][l_dei]]
        mass = [as_pars[i][a_mi], -1.]
        smass = [as_pars[i][a_mei], -1.]
        rads = [as_pars[i][a_rad], cl_dict[names_idx[i]][l_rad]]

        # Store ASteCA values (k=0) and literature values (k=1).
        for k in [0, 1]:

            # Convert z to [Fe/H] for ASteCA values.
            if k == 0:
                z, ez = float_str(met[k]), float_str(smet[k])
                if z > 0.:
                    fe_h = np.log10(z / 0.0152)
                else:
                    # Use minimum metallicity value if z=0.
                    fe_h = np.log10(0.0001 / 0.0152)
                e_fe_h = (1. / np.log(10.)) * z / ez
                # Store.
                zarr[j][k].append(fe_h)
                zsigma[j][k].append(e_fe_h)
            else:
                zarr[j][k].append(float_str(met[k]))
                zsigma[j][k].append(float_str(smet[k]))

            aarr[j][k].append(float_str(age[k]))
            asigma[j][k].append(float_str(sage[k]))
            earr[j][k].append(float_str(ext[k]))
            esigma[j][k].append(float_str(sext[k]))
            darr[j][k].append(float_str(dis[k]))
            dsigma[j][k].append(float_str(sdis[k]))
            marr[j][k].append(float_str(mass[k]))
            msigma[j][k].append(float_str(smass[k]))
            rarr[j][k].append(float_str(rads[k]))

    # Pass as dict.
    pars_dict = {'gal_names': gal_names, 'ra': ra, 'dec': dec, 'zarr': zarr,
        'zsigma': zsigma, 'aarr': aarr, 'asigma': asigma, 'earr': earr,
        'esigma': esigma, 'darr': darr, 'dsigma': dsigma, 'marr': marr,
        'msigma': msigma, 'rarr': rarr, 'ext_sf': ext_sf, 'ext_mcev': ext_mcev,
        'rad_pc': rad_pc, 'int_colors': int_colors, 'n_memb': n_memb}

    return pars_dict