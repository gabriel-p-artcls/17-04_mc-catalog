
import numpy as np
from functions.photom_dispersion import get_disp
from functions.dist_2_cent import dist_2_cloud_center


def float_str(val):
    '''
    Return float depending on whether the value is a float or a string.
    '''
    try:
        val_f = float(val)
    except:
        val_f = -9999999999.9

    return val_f


def rad_in_pc(float_lst):
    '''
    Convert radius from pixel to arcsec to parsecs.
    '''
    # Unpack
    r_px, scale, dis_mod, e_bv = float_lst
    # Convert from px to arcsecs.
    rad_arcsec = r_px * scale
    Av = 3.1 * e_bv
    # Distance to cluster in parsecs.
    d_pc = 10 ** (0.2 * (dis_mod + 5 - Av))
    # Radius in parsecs.
    r_pc = d_pc * np.tan(np.deg2rad(rad_arcsec / 3600.))

    return r_pc


def correct_int_col_extin(int_col, extinc):
    '''
    Correct integrated color for the extinction.
    '''
    E_CT1_E_BV = 1.97
    int_col_cor = int_col - E_CT1_E_BV * float_str(extinc)

    return int_col_cor


def z_to_feh(z, ez):
    '''
    Convert z to [Fe/H] for ASteCA values.
    '''
    # Use minimum metallicity value if z=0.
    z = max(0.0001, z)
    fe_h = np.log10(z / 0.0152)
    e_fe_h = (1. / np.log(10.)) * (ez / z)
    # Trim error if it's too large.
    e_fe_h = min(e_fe_h, 2.)

    return fe_h, e_fe_h


def params(r_path, as_names, as_pars, cl_dict, names_idx):
    '''
    Return ASteCA output and literature parameters values.
    '''
    # Indexes of columns in ASteCA output file.
    a_zi, a_zei, a_ai, a_aei, a_ei, a_eei, a_di, a_dei, a_mi, a_mei, a_rad, \
        a_erad, a_int_c, a_nmemb, a_CI, a_prob, a_r_core, a_e_r_core =\
        19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 4, 5, 18, 13, 11, 17, 6, 7

    # Indexes of columns in .ods literature file.
    ra_i, dec_i, gal_i, l_zi, l_zei, l_ai, l_aei, l_ei, l_eei, l_di, l_dei, \
        l_rad, l_scale, l_e_sandf, l_e_e_sandf, l_e_mcev, l_e_mcev_max, \
        l_e_e_mcev, l_mcev_dist, l_mass, l_e_mass = \
        cl_dict[0].index(u'ra_deg'), cl_dict[0].index(u'dec_deg'), \
        cl_dict[0].index(u'Galaxia'), cl_dict[0].index(u'[Fe/H] (dex)'), \
        cl_dict[0].index(u'e_Fe/H'), cl_dict[0].index(u'log(age)'), \
        cl_dict[0].index(u'e_log(age)'), cl_dict[0].index(u'E(B-V) (lit)'), \
        cl_dict[0].index(u'e_E(B-V)'), cl_dict[0].index(u'(m-M)o (mag)'), \
        cl_dict[0].index(u'e_(m-M)o'), cl_dict[0].index(u'rad (eye)'), \
        cl_dict[0].index(u'arcsec/pixel'), cl_dict[0].index(u'E_B_V_SandF'), \
        cl_dict[0].index(u'stdev_E_B_V_SandF'), \
        cl_dict[0].index(u'E_BV_closer_MCEV'), cl_dict[0].index(u'E_BV_max'), \
        cl_dict[0].index(u'E_BV_std_dev'), cl_dict[0].index(u'Dist (deg)'), \
        cl_dict[0].index(u'Mass'), cl_dict[0].index(u'e_mass')

    # Initialize empty lists. The first sub-list in the parameters list
    # corresponds to clusters in the SMC and the second to those in the LMC.
    gal_names, int_colors, n_memb, cont_ind, kde_prob, ra, dec, rad_pc, \
        erad_pc, dist_cent, e_d_cent, r_core_pc, e_r_core, phot_disp = \
        [[[], []] for _ in range(14)]
    # First sub-list stores SMC values, the second one stores LMC values.
    # First and 2nd sub-sublist store Schlafly & Finkbeiner extinction values
    # and their errors. Third and 4th store MCEV extinction values and their
    # errors.
    ext_sf = [[[], []], [[], []]]
    ext_mcev = [[[], [], [], []], [[], [], [], []]]
    # First sub-list stores ASteCA values, the second one stores literature
    # values.
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
        int_col_no_corr = float_str(as_p[a_int_c])
        # Correct for extinction.
        int_col_corr = correct_int_col_extin(int_col_no_corr, as_p[a_ei])
        # Store extinction corrected integrated color.
        int_colors[j].append(int_col_corr)
        # Approx number of members.
        n_memb[j].append(float_str(as_p[a_nmemb]))
        # Contamination index.
        cont_ind[j].append(float_str(as_p[a_CI]))
        # Cluster KDE p-value probability.
        kde_prob[j].append(float_str(as_p[a_prob]))

        # Store literature E(B-V) values: Schlafly & Finkbeiner (SandF) and
        # MCEV.
        ext_sf[j][0].append(cl_dict[names_idx[i]][l_e_sandf])
        ext_sf[j][1].append(cl_dict[names_idx[i]][l_e_e_sandf])
        ext_mcev[j][0].append(cl_dict[names_idx[i]][l_e_mcev])
        ext_mcev[j][1].append(cl_dict[names_idx[i]][l_e_mcev_max])
        ext_mcev[j][2].append(cl_dict[names_idx[i]][l_e_e_mcev])
        ext_mcev[j][3].append(cl_dict[names_idx[i]][l_mcev_dist])

        # Store radius value in parsecs.
        float_lst = []
        for el in [as_p[a_rad], cl_dict[names_idx[i]][l_scale], as_p[a_di],
                   as_p[a_ei]]:
            # Store in list as floats.
            float_lst.append(float_str(el))
        r_pc = rad_in_pc(float_lst)
        rad_pc[j].append(r_pc)
        # Repeat process for errors in radius.
        float_lst = []
        for el in [as_p[a_erad], cl_dict[names_idx[i]][l_scale], as_p[a_di],
                   as_p[a_ei]]:
            # Store in list as floats.
            float_lst.append(float_str(el))
        e_r_pc = rad_in_pc(float_lst)
        erad_pc[j].append(e_r_pc)
        # Calculate r_core and its error in parsecs, using a simple 3 rule.
        r_core_px = float_str(as_p[a_r_core])
        e_r_core_px = float_str(as_p[a_e_r_core])
        px_2_pc_scale = r_pc / float_str(as_p[a_rad])
        if r_core_px > 0.:
            r_c_pc = r_core_px * px_2_pc_scale
            e_r_c_pc = e_r_core_px * px_2_pc_scale
        else:
            # Dummy values for clusters with no r_core values.
            r_c_pc, e_r_c_pc = -10., 0.
        # Store core radius and its error in parsecs.
        r_core_pc[j].append(r_c_pc)
        e_r_core[j].append(e_r_c_pc)
        # if r_c_pc > r_pc:
        #     print j, as_names[i], as_p[a_rad], r_core_px

        # Get 3D distance from cluster to galaxy center, and its error. Both
        # values are expressed in parsecs.
        d_c, e_d_c = dist_2_cloud_center(j, ra[j][-1], dec[j][-1],
                                         as_p[a_di], float(as_p[a_dei]))
        dist_cent[j].append(d_c.value)
        e_d_cent[j].append(e_d_c.value)

        # Get photometric dispersion parameter.
        phot_disp[j].append(get_disp(r_path, as_names[i]))
        # phot_disp[j].append(0.)

        # Organize param values, ASteCA first, lit second.
        met = [as_p[a_zi], cl_dict[names_idx[i]][l_zi]]
        smet = [as_p[a_zei], cl_dict[names_idx[i]][l_zei]]
        age = [as_p[a_ai], cl_dict[names_idx[i]][l_ai]]
        sage = [as_p[a_aei], cl_dict[names_idx[i]][l_aei]]
        ext = [as_p[a_ei], cl_dict[names_idx[i]][l_ei]]
        sext = [as_p[a_eei], cl_dict[names_idx[i]][l_eei]]
        dis = [as_p[a_di], cl_dict[names_idx[i]][l_di]]
        sdis = [as_p[a_dei], cl_dict[names_idx[i]][l_dei]]
        mass = [as_p[a_mi], cl_dict[names_idx[i]][l_mass]]
        smass = [as_p[a_mei], cl_dict[names_idx[i]][l_e_mass]]
        rads = [as_p[a_rad], cl_dict[names_idx[i]][l_rad]]

        # Store ASteCA values (k=0) and literature values (k=1).
        for k in [0, 1]:

            if k == 0:
                z, ez = float_str(met[k]), float_str(smet[k])
                # Convert z to [Fe/H] for ASteCA values.
                fe_h, e_fe_h = z_to_feh(z, ez)
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

    # Pass as dictionary.
    pars_dict = {
        'gal_names': gal_names, 'ra': ra, 'dec': dec, 'zarr': zarr,
        'zsigma': zsigma, 'aarr': aarr, 'asigma': asigma, 'earr': earr,
        'esigma': esigma, 'darr': darr, 'dsigma': dsigma, 'marr': marr,
        'msigma': msigma, 'rarr': rarr, 'ext_sf': ext_sf, 'ext_mcev': ext_mcev,
        'rad_pc': rad_pc, 'erad_pc': erad_pc, 'int_colors': int_colors,
        'n_memb': n_memb, 'cont_ind': cont_ind, 'kde_prob': kde_prob,
        'dist_cent': dist_cent, 'e_d_cent': e_d_cent, 'r_core_pc': r_core_pc,
        'e_r_core': e_r_core, 'phot_disp': phot_disp
    }

    return pars_dict
