
import numpy as np
from scipy.stats import ks_2samp, pearsonr


def ccc(l1, l2):
    '''
    Concordance correlation coefficient.
    See: https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    '''
    return 2 * np.cov(l1, l2)[0, 1] / (np.var(l1) + np.var(l2) +
                                       (np.mean(l1) - np.mean(l2)) ** 2)


def check_diffs(in_params):
    '''
    check differences between ASteCA values and literature values for given
    parameters.
    '''
    gal_names, zarr, zsigma, aarr, asigma, earr, darr, dsigma, rarr, marr,\
        dist_cent, e_d_cent, ra, dec, n_memb, rad_pc, int_colors, cont_ind,\
        phot_disp = \
        [in_params[_] for _ in ['gal_names', 'zarr', 'zsigma', 'aarr',
                                'asigma', 'earr', 'darr', 'dsigma', 'rarr',
                                'marr', 'dist_cent', 'e_d_cent', 'ra', 'dec',
                                'n_memb', 'rad_pc', 'int_colors', 'cont_ind',
                                'phot_disp']]

    gal = ['SMC', 'LMC']
    age_xmc_f_all = []
    # For SMC and LMC.
    for j in [0, 1]:
        print '\n*** {} ***\n'.format(gal[j])

        # print 'Contamination index'
        # for i, name in enumerate(gal_names[j]):
        #     print name, cont_ind[j][i]

        print '\n{} clusters in age/rad range:'.format(gal[j])
        for i, name in enumerate(gal_names[j]):
            a, r = aarr[j][0][i], rad_pc[j][i]
            if a < 8.5 and r > 12.5:
                print '{}: age: {} ; rad: {} pc'.format(name, a, r)
        print ''

        # For each cluster.
        met_count = 0
        for i, name in enumerate(gal_names[j]):
            # Metallicity.
            z_diff = 0.75
            diff = zarr[j][0][i] - zarr[j][1][i]
            if zarr[j][1][i] > -99.:
                if abs(diff) > z_diff:
                    met_count += 1
                    # AsteCA minus Literature metallicity difference.
                    rel_diff = zarr[j][0][i] - zarr[j][1][i]
                    print '{} {}, {:.2f} vs {:.2f} , {:.2f}'.format(
                        gal[j], name, zarr[j][0][i], zarr[j][1][i], rel_diff)

        print '{}, Clusters with \delta z>{}: {}\n'.format(
            gal[j], z_diff, met_count)

        # For each cluster.
        age_diff = []
        print 'Gal-Clust alpha delta log(age)_lit log(age)_ASteCA'
        for i, name in enumerate(gal_names[j]):
            # Age.
            a_diff = 0.5
            # ASteCA - literature
            diff = aarr[j][0][i] - aarr[j][1][i]
            if aarr[j][1][i] > -99.:
                if abs(diff) > a_diff:
                    age_diff.append(diff)
                    # Literature minus AsteCA Log(age) difference.
                    rel_diff = aarr[j][1][i] - aarr[j][0][i]
                    print ("{}-{} & {:.5f} & {:.5f} & {:.2f} & "
                           "{:.2f} & {:.2f}\\\\".format(
                            gal[j][0], name, ra[j][i], dec[j][i],
                            aarr[j][1][i], aarr[j][0][i], rel_diff))

        print '{}, Clusters with \delta log(age)>{}: {}\n'.format(
            gal[j], a_diff, len(age_diff))

        if j == 0:
            avrg_mass = []
            print 'Masses for SMC clusters: ASteCA - Maia et al. (2013) = diff'
            for i, (ma, ml) in enumerate(zip(*[marr[0][0], marr[0][1]])):
                if abs(ml) < 5000:
                    print '{}: {} - {} = {}, {}'.format(
                        gal_names[j][i], ma, ml, ma-ml, cont_ind[0][i])
                    m_lim = 1500.
                    if abs(ma-ml) < m_lim:
                        avrg_mass.append(ma-ml)
            print 'Mean mass diff for Delta<{}: {} +- {}'.format(
                m_lim, np.mean(avrg_mass), np.std(avrg_mass))

        print ''
        err_c, err_thresh = 0, 0.2
        for ez in zsigma[j][0]:
            if ez <= err_thresh:
                err_c += 1
        perc = float(err_c)/len(zsigma[j][0])
        print 'Perc of OC with z errors below {}: {}'.format(err_thresh, perc)

        err_c, err_min, err_max = 0, 0.0029, 0.0031
        for i, e in enumerate(zsigma[j][0]):
            # convert from [Fe/H] to z
            z = 10**zarr[j][0][i] * 0.0152
            e_z = e*z*np.log(10.)
            if err_min <= e_z <= err_max:
                err_c += 1
        perc = float(err_c)/len(zsigma[j][0])
        print 'Perc of OC with {}<= z <={}: {}'.format(err_min, err_max,
                                                       perc)

        default_fe_h, feh_sum = [-0.7, -0.4], 0.
        for i, fe_h_lit in enumerate(zarr[j][1]):
            feh = default_fe_h[j]
            if (feh - 0.01) <= fe_h_lit <= (feh + 0.01):
                feh_sum += 1
        tot_feh = len(zarr[j][1]) if j == 0 else (len(zarr[j][1]) - 36)
        perc = float(feh_sum)/tot_feh
        print 'Perc of OC with lit values [Fe/H]~{}: {}'.format(
            feh, perc)

        err_c, err_thresh = 0, 0.1
        for e in asigma[j][0]:
            if e <= err_thresh:
                err_c += 1
        perc = float(err_c)/len(asigma[j][0])
        print 'Perc of OC with age errors below {}: {}\n'.format(err_thresh,
                                                                 perc)

        min_rgc, max_rgc = 3500., 5500.
        min_a, max_a = 7.5, 8.5
        for i, a in enumerate(aarr[j][0]):
            R_gc = dist_cent[j][i]
            if min_rgc < R_gc < max_rgc and min_a < a < max_a:
                print 'OCs in {} < R_gc < {} & {} < age < {}: {}'.format(
                    min_rgc, max_rgc, min_a, max_a, gal_names[j][i])
                print '  OC R_gc error:', e_d_cent[j][i]

        print ''
        e_rgc_max = 100.
        for i, R_gc in enumerate(dist_cent[j]):
            e_rgc = e_d_cent[j][i]
            if e_rgc < e_rgc_max:
                print ("OCs w e_R_gc < {}: {} / ra, dec, dm, e_dm: "
                       "{}, {}, {}, {}".format(
                            e_rgc_max, gal_names[j][i], ra[j][i],
                            dec[j][i], darr[j][0][i], dsigma[j][0][i]))

        print '\nAverage ASteCA E(B-V) for the {}: {} +- {}'.format(
            gal[j], np.mean(earr[j][0]), np.std(earr[j][0]))
        print 'Average ASteCA mass for the {}: {}'.format(
            gal[j], np.mean(marr[j][0]))
        print 'Average ASteCA radius for the {}: {}'.format(
            gal[j], np.mean(rarr[j][0]))
        print 'Average ASteCA density for the {}: {}'.format(
            gal[j], np.mean(n_memb[j]/(np.pi*np.array(rarr[j][0])**2)))

        print '\nMean literature e_log(age) for {}: {}\n'.format(
            gal[j], np.mean(asigma[j][1]))

        # Mean only for those clusters with ASteCA age values closer than 0.5
        # to literature values.
        age_xmc_f = []
        for a in list(np.array(aarr[j][0]) - np.array(aarr[j][1])):
            if abs(a) <= 0.5:
                age_xmc_f.append(a)
        age_xmc_f_all += age_xmc_f
        print 'Age mean for Delta log(age)<0.5:', np.mean(age_xmc_f)
        print 'Met vals mean/std, AS:', np.mean(zarr[j][0]), \
            np.std(zarr[j][0])
        # Filter out clusters with no metal values in the literature
        # (.ods file)
        z_xmc_lit_f, z_xmc_ast_f = [], []
        for z_ast, z_lit in zip(*[zarr[j][0], zarr[j][1]]):
            if abs(z_lit) < 10000:
                z_xmc_lit_f.append(z_lit)
                z_xmc_ast_f.append(z_ast)
        print 'Met vals mean/std, Lit:', np.mean(z_xmc_lit_f), \
            np.std(z_xmc_lit_f)
        print 'Met vals mean diff:', np.mean(np.array(z_xmc_ast_f) -
                                             np.array(z_xmc_lit_f))
        print 'Met vals median diff:', np.median(np.array(z_xmc_ast_f) -
                                                 np.array(z_xmc_lit_f)), '\n'

        # K_S test.
        # Null hypothesis: that 2 independent samples are drawn from the same
        # continuous distribution (sample sizes can be different)
        #
        # If the K-S statistic is small or the p-value is high, then we cannot
        # reject the hypothesis that the distributions of the two samples are
        # the same.
        # For two identical distributions the KS value will be small and
        # the p-value high.
        print 'Met vals CCC:', ccc(z_xmc_ast_f, z_xmc_lit_f)
        print 'Met vals PCC:', np.corrcoef(z_xmc_ast_f, z_xmc_lit_f)[0, 1]
        ks, pval = ks_2samp(z_xmc_ast_f, z_xmc_lit_f)
        print 'Met vals K-S:', ks, pval
        print 'Age vals CCC:', ccc(aarr[j][0], aarr[j][1])
        print 'Age vals PCC:', np.corrcoef(aarr[j][0], aarr[j][1])[0, 1]
        ks, pval = ks_2samp(aarr[j][0], aarr[j][1])
        print 'Age vals K-S:', ks, pval, '\n'

        print 'Filter out Piatti (2011) clusters that only have ages assigned'
        par_all, par_delta = [[[], [], [], []] for _ in range(2)]
        for i, param in enumerate([zarr, aarr, earr, darr]):
            for k, p_lit in enumerate(param[j][1]):
                # Filter out Piatti (2011) clusters that only have ages
                # assigned.
                if abs(zarr[j][1][k]) < 30000.:
                    # \delta as: ASteCA - literature values.
                    # par_all[i].append(param[j][0][k])
                    par_delta[i].append(param[j][0][k] - param[j][1][k])
        p_name = ['[Fe/H]', 'log(age)', 'E(B-V)', 'dm']
        for i, span in enumerate(par_delta):
            print 'Delta {}/{} mean+-std: {:.3f} +- {:.3f}'.format(
                gal[j], p_name[i], np.mean(span), np.std(span))
        for i, span in enumerate(par_delta[:-1]):
            r_pears = pearsonr(span, par_delta[-1])
            print 'Correlation Delta {} vs dm: {:.3f}'.format(
                p_name[i], r_pears[0])
        print ''

    print '\nAge mean for Delta log(age)<0.5 S/LMC:{}\n'.format(
        np.mean(age_xmc_f_all))
