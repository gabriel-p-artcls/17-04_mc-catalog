
import numpy as np
from kde_2d import kde_map


def age_met_rel(xarr, xsigma, yarr, ysigma):
    '''
    Generate an age-metallicity relation (AMR) based on the KDE of the
    age-metallicity 2D space and a weighted average of the metallicity for
    several age values.
    '''

    # Define extension where the KDE will be obtained.
    x_min, x_max, y_min, y_max = 0., 6., -2.5, 0.
    ext = [x_min, x_max, y_min, y_max]
    # Grid density.
    gd = 20
    # Generate metallicity values as in grid. Invert list so the weighted
    # average is obtained correctly.
    met_vals = np.linspace(y_min, y_max, gd)[::-1]
    age_vals = np.linspace(x_min, x_max, gd)

    # yarr = -0.5 + 0.0001 * np.random.randn(len(yarr))
    # xsigma, ysigma = [0.05] * len(xarr), [0.1] * len(xarr)

    # Obtain age-metallicity KDE for the entire range.
    z = kde_map(np.array(xarr), np.array(xsigma), np.array(yarr),
                np.array(ysigma), ext, gd)
    # Order KDE in age columns.
    a_m_kde = zip(*z)

    # Obtain metallicity weighted average for each age value.
    met_weighted = [[], []]
    for age_col in a_m_kde:
        age_col = np.array(age_col)
        # met_w = sum(np.array(met_vals[(i + 1):-(i + 1)]) *
        #             np.array(age_col)) / sum(age_col)
        sum_age_col = sum(age_col)
        met_w = sum(met_vals * age_col) / sum_age_col
        met_weighted[0].append(met_w)
        # Obtain standard deviation.
        nume = sum(age_col * (met_w - met_vals) ** 2)
        deno = sum_age_col - (sum(age_col ** 2) / sum_age_col)
        stdev_met_w = np.sqrt(nume / deno)
        met_weighted[1].append(stdev_met_w)

    # Store list of weighted average value for the metallicity obtained for
    # each defined age, *for this KDE area*.
    # met_vals_age.append(met_weighted)

    # Average values to obtain final AMR. It consists of one metallicity
    # point for each age value defined, along with its error.
    # amr = []
    # for area_met_vals in zip(*met_vals_age):
    #     amr.append([np.mean(area_met_vals), np.std(area_met_vals)])
    # print amr
    # print zip(*met_weighted)

    return age_vals, met_weighted


def main(in_params):
    '''
    '''

    zarr, zsigma, aarr, asigma = [in_params[_] for _ in ['zarr', 'zsigma',
                                                         'aarr', 'asigma']]

    # First index k indicates the galaxy (0 for SMC, 1 for KMC), the second
    # index 0 indicates ASteCA values.
    k = 0  # SMC
    # Age in Gyrs.
    age_gyr_smc = [10 ** (np.asarray(aarr[k][0]) - 9),
                   np.asarray(asigma[k][0]) * np.asarray(aarr[k][0]) *
                   np.log(10) / 5.]
    age_vals_smc, met_weighted_smc = age_met_rel(
        age_gyr_smc[0], age_gyr_smc[1], zarr[k][0], zsigma[k][0])

    k = 1  # LMC
    # Age in Gyrs.
    age_gyr_lmc = [10 ** (np.asarray(aarr[k][0]) - 9),
                   np.asarray(asigma[k][0]) * np.asarray(aarr[k][0]) *
                   np.log(10) / 5.]
    age_vals_lmc, met_weighted_lmc = age_met_rel(
        age_gyr_lmc[0], age_gyr_lmc[1], zarr[k][0], zsigma[k][0])

    import matplotlib.pyplot as plt
    plt.xlim(0., 6)
    plt.ylim(-2.2, 0.1)
    plt.grid(b=True, which='major', color='gray', linestyle='--', lw=0.5,
             zorder=1)
    plt.plot(age_vals_smc, met_weighted_smc[0], c='r')
    plt.scatter(age_gyr_smc[0], zarr[0][0], c='r', s=25, marker='s')
    plt.plot(age_vals_lmc, met_weighted_lmc[0], c='b')
    plt.scatter(age_gyr_lmc[0], zarr[1][0], c='b', s=25, marker='o')
    y_err_min_smc = np.array(met_weighted_smc[0]) - \
        np.array(met_weighted_smc[1])
    y_err_max_smc = np.array(met_weighted_smc[0]) + \
        np.array(met_weighted_smc[1])
    plt.fill_between(age_vals_smc, y_err_min_smc, y_err_max_smc, alpha=0.1,
                     color='red')
    y_err_min_lmc = np.array(met_weighted_lmc[0]) - \
        np.array(met_weighted_lmc[1])
    y_err_max_lmc = np.array(met_weighted_lmc[0]) + \
        np.array(met_weighted_lmc[1])
    plt.fill_between(age_vals_lmc, y_err_min_lmc, y_err_max_lmc, alpha=0.1,
                     color='blue')
    plt.show()


if __name__ == "__main__":
    main()
