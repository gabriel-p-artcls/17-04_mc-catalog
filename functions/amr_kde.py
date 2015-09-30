
import numpy as np
from kde_2d import kde_map


def age_met_rel(xarr, xsigma, yarr, ysigma):
    '''
    Generate an age-metallicity relation based on the KDE of the
    age-metallicity 2D space and a weighted average of the metallicity for
    several age values.
    '''

    x_min, x_max, y_min, y_max = 0., 6., -2.5, 0.2

    # Define extension where the KDE will be obtained.
    ext = [x_min, x_max, y_min, y_max]

    # Grid density.
    gd = 4

    # Obtain age-metallicity KDE.
    z = kde_map(np.array(xarr), np.array(xsigma), np.array(yarr),
                np.array(ysigma), ext, gd)
    a_m_kde = zip(*z)
    print '\na_m_kde', a_m_kde
    met_vals = np.linspace(y_min, y_max, gd)
    print 'x_grid', np.linspace(x_min, x_max, gd)
    print 'y_grid', met_vals, '\n'

    # Obtain metallicity weighted average for each age value.
    met_vals_age = []
    for age_col in a_m_kde:
        print 'met in age col', age_col
        met_w = sum(np.array(met_vals) * np.array(age_col)) / sum(age_col)
        print 'weight met', met_w
        met_vals_age.append(met_w)

    raw_input()
    # Average values to obtain final AMR. It consists in one metallicity
    # point for each age value defined, along with its error.

    return amr


def main(gal, k, in_params):
    '''
    '''

    zarr, zsigma, aarr, asigma = [in_params[_] for _ in ['zarr', 'zsigma',
                                                         'aarr', 'asigma']]

    # First index k indicates the galaxy (0 for SMC, 1 for KMC), the second
    # index 0 indicates ASteCA values.
    amr = age_met_rel(aarr[k][0], asigma[k][0], zarr[k][0], zsigma[k][0])

if __name__ == "__main__":
    main()
