
import numpy as np
from kde_map import kde_2d


def age_met_rel(xarr, xsigma, yarr, ysigma):
    '''
    Generate an age-metallicity relation (AMR) based on the KDE of the
    age-metallicity 2D space and a weighted average of the metallicity for
    several age values.

    See: http://math.stackexchange.com/q/1457390/37846
    '''

    # Define 2D space extension where the KDE will be obtained.
    x_min, x_max = min(np.array(xarr) - np.array(xsigma)), \
        max(np.array(xarr) + 0.5*np.array(xsigma))
    # These range is important since it defines the metallicity values
    # (met_vals) that will be weighted by the KDE below.
    y_min, y_max = min(np.array(yarr) - np.array(ysigma)), \
        max(np.array(yarr) + np.array(ysigma))

    # Metallicity step. THIS NUMBER AFFECTS THE SHAPE OF THE FINAL AMR.
    # We select a value of 0.3, which gives steps ~0.3 dex in the metallicity
    # grid below. This value is very similar to the average uncertainty in
    # [Fe/H]. A value that is too large (i.e. ~1.) results in a flat AMR
    # shifted towards lower [Fe/H] values. A value that is too low (i.e. ~0.01)
    # results in a noisy AMRshifted towards larger values.
    met_step = 0.3
    # Grid density.
    gd = int((y_max - y_min) / met_step)

    # Generate metallicity values as in grid. Invert list so the weighted
    # average is obtained correctly.
    met_vals = np.linspace(y_min, y_max, gd)[::-1]
    age_vals = np.linspace(x_min, x_max, gd)

    # yarr = -0.5 + 0.0001 * np.random.randn(len(yarr))
    # xsigma, ysigma = [0.05] * len(xarr), [0.1] * len(xarr)

    # Obtain age-metallicity KDE for the entire defined range.
    ext = [x_min, x_max, y_min, y_max]
    z = kde_2d(np.array(xarr), np.array(xsigma), np.array(yarr),
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

        # Metallicity values given by the KDE for the entire metallicity range,
        # for a single age value.
        age_col = np.array(age_col)

        # Obtain weighted metallicity for this *single* age value.
        sum_age_col = sum(age_col)
        met_w = sum(met_vals * age_col) / sum_age_col
        met_weighted[0].append(met_w)

        # Obtain standard deviation.
        nume = sum(age_col * (met_w - met_vals) ** 2)
        deno = sum_age_col - (sum(age_col ** 2) / sum_age_col)
        stdev_met_w = np.sqrt(nume / deno)
        met_weighted[1].append(stdev_met_w)

    return age_vals, met_weighted
