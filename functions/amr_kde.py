
import numpy as np
from kde_2d import kde_map


def age_met_rel(xarr, xsigma, yarr, ysigma):
    '''
    Generate an age-metallicity relation (AMR) based on the KDE of the
    age-metallicity 2D space and a weighted average of the metallicity for
    several age values.

    See: http://math.stackexchange.com/q/1457390/37846
    '''

    # Define 2D space extension where the KDE will be obtained.
    x_min, x_max = min(np.array(xarr) - np.array(xsigma)), \
        max(np.array(xarr) + np.array(xsigma))
    y_min, y_max = min(np.array(yarr) - np.array(ysigma)), \
        max(np.array(yarr) + np.array(ysigma))
    # x_min, x_max, y_min, y_max = -4., 10., -7., 5.
    ext = [x_min, x_max, y_min, y_max]
    # Grid density.
    # Metallicity step.
    met_step = 0.02
    gd = int((y_max - y_min) / met_step)

    # Generate metallicity values as in grid. Invert list so the weighted
    # average is obtained correctly.
    met_vals = np.linspace(y_min, y_max, gd)[::-1]
    age_vals = np.linspace(x_min, x_max, gd)

    # yarr = -0.5 + 0.0001 * np.random.randn(len(yarr))
    # xsigma, ysigma = [0.05] * len(xarr), [0.1] * len(xarr)

    # Obtain age-metallicity KDE for the entire range.
    z = kde_map(np.array(xarr), np.array(xsigma), np.array(yarr),
                np.array(ysigma), ext, gd)
    # Order KDE in age columns where each column is associated with an age.
    a_m_kde = zip(*z)

    # Obtain metallicity weighted average for each age value.
    met_weighted = [[], []]
    for age_col in a_m_kde:
        # # Filter out points with very small values (assign 0. value)
        # min_w = max(age_col) / 20.
        # age_col2 = []
        # N = 0
        # for _ in age_col:
        #     if _ < min_w:
        #         age_col2.append(0.)
        #         N += 1
        #     else:
        #         age_col2.append(_)
        # age_col = np.array(age_col2)
        age_col = np.array(age_col)
        # Obtain weighted metallicity for this age value.
        sum_age_col = sum(age_col)
        met_w = sum(met_vals * age_col) / sum_age_col
        met_weighted[0].append(met_w)
        # Obtain standard deviation.
        nume = sum(age_col * (met_w - met_vals) ** 2)
        deno = sum_age_col - (sum(age_col ** 2) / sum_age_col)
        stdev_met_w = np.sqrt(nume / deno)
        met_weighted[1].append(stdev_met_w)

    return age_vals, met_weighted
