
import numpy as np


def check_diffs(in_params):
    '''
    check differences between ASteCA values and literature values for given
    parameters.
    '''
    gal_names, zarr, zsigma, aarr, asigma, earr, darr, rarr, marr, dist_cent,\
        ra, dec, n_memb = \
        [in_params[_] for _ in ['gal_names', 'zarr', 'zsigma', 'aarr',
                                'asigma', 'earr', 'darr', 'rarr', 'marr',
                                'dist_cent', 'ra', 'dec', 'n_memb']]

    gal = ['SMC', 'LMC']
    print ''

    # For SMC and LMC.
    for j in [0, 1]:

        # For each cluster.
        met_count = 0
        for i, name in enumerate(gal_names[j]):
            # Metallicity.
            z_diff = 0.5
            diff = zarr[j][0][i] - zarr[j][1][i]
            if zarr[j][1][i] > -99.:
                if abs(diff) > z_diff:
                    met_count += 1
                    # AsteCA minus Literature metallicity difference.
                    rel_diff = zarr[j][0][i] - zarr[j][1][i]
                    print '{} {}, {:.2f} vs {:.2f} , {:.2f}'.format(
                        gal[j], name, zarr[j][0][i], zarr[j][1][i], rel_diff)

        print '\n* {}, Clusters with \delta z>{}: {}\n'.format(
            gal[j], z_diff, met_count)

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

        # For each cluster.
        age_diff = []
        for i, name in enumerate(gal_names[j]):
            # Age.
            a_diff = 0.5
            # ASteCA - literature
            diff = aarr[j][0][i] - aarr[j][1][i]
            if aarr[j][1][i] > -99.:
                if abs(diff) > a_diff:
                    age_diff.append(diff)
                    # AsteCA minus Literature Log(age) difference.
                    rel_diff = aarr[j][0][i] - aarr[j][1][i]
                    print '{} {}, {:.2f} vs {:.2f} , {:.2f}'.format(
                        gal[j], name, aarr[j][0][i], aarr[j][1][i], rel_diff)

        print '\n* {}, Clusters with \delta log(age)>{}: {}\n'.format(
            gal[j], a_diff, len(age_diff))

        err_c, err_thresh = 0, 0.1
        for e in asigma[j][0]:
            if e <= err_thresh:
                err_c += 1
        perc = float(err_c)/len(asigma[j][0])
        print 'Perc of OC with age errors below {}: {}'.format(err_thresh,
                                                               perc)

        print 'Average mass for the {}: {}'.format(gal[j], np.mean(marr[j][0]))
        print 'Average radius for the {}: {}'.format(gal[j],
                                                     np.mean(rarr[j][0]))
        print 'Average density for the {}: {}'.format(
            gal[j], np.mean(n_memb[j]/(np.pi*np.array(rarr[j][0])**2)))
