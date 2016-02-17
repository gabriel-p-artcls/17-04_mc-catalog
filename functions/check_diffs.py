

def check_diffs(in_params):
    '''
    check differences between ASteCA values and literature values for given
    parameters.
    '''
    gal_names, zarr, aarr, earr, darr, rarr, marr, dist_cent, ra, dec = \
        [in_params[_] for _ in ['gal_names', 'zarr', 'aarr', 'earr', 'darr',
                                'rarr', 'marr', 'dist_cent', 'ra', 'dec']]

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

        # For each cluster.
        age_count = 0
        for i, name in enumerate(gal_names[j]):
            # Age.
            a_diff = 0.5
            # ASteCA - literature
            diff = aarr[j][0][i] - aarr[j][1][i]
            if aarr[j][1][i] > -99.:
                if abs(diff) > a_diff:
                    age_count += 1
                    # AsteCA minus Literature Log(age) difference.
                    rel_diff = aarr[j][0][i] - aarr[j][1][i]
                    print '{} {}, {:.2f} vs {:.2f} , {:.2f}'.format(
                        gal[j], name, aarr[j][0][i], aarr[j][1][i], rel_diff)

        print '\n* {}, Clusters with \delta log(age)>{}: {}\n'.format(
            gal[j], a_diff, age_count)
