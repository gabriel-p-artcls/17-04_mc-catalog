

import numpy as np


def kernel(p0, p, s):
    '''
    Gaussian kernel.
    '''
    # Replace 0 error with very small value.
    s = 0.000001 if s < 0. else s
    val = (1. / s) * np.exp(-0.5 * ((p0 - p) / s)**2)

    return val


def kde_val(point, xarr, xsigma, yarr, ysigma, dim):
    '''
    Returns the value of the 1D/2D Kernel Density Estimator (KDE) evaluated in
    the point passed.
    '''
    # Evaluate KDE in 'point'.
    kde_p = 0.
    if dim == '1':
        for x, sx in zip(xarr, xsigma):
            kde_p = kde_p + kernel(point[0], x, sx)
        # Normalize.
        norm_kde = kde_p / (np.sqrt(2 * np.pi) * len(xarr))
    elif dim == '2':
        for x, sx, y, sy in zip(xarr, xsigma, yarr, ysigma):
            kde_p = kde_p + kernel(point[0], x, sx) * kernel(point[1], y, sy)
        # Normalize.
        norm_kde = kde_p / (2 * np.pi * len(xarr))

    return norm_kde


def kde_2d(xarr, xsigma, yarr, ysigma, ext, grid_dens):
    '''
    Take an array of x,y data with their errors, create a grid of points in x,y
    and return the 2D KDE density map.
    '''

    # Grid density (number of points).
    gd_c = complex(0, grid_dens)

    # Define grid of points in x,y where the KDE will be evaluated.
    x, y = np.mgrid[ext[0]:ext[1]:gd_c, ext[2]:ext[3]:gd_c]
    positions = np.vstack([x.ravel(), y.ravel()])

    # Evaluate KDE in x,y grid.
    kde_grid = []
    for p in zip(*positions):
        kde_grid.append(kde_val(p, xarr, xsigma, yarr, ysigma, '2'))

    # Re-shape values for plotting.
    z = np.rot90(np.reshape(np.array(kde_grid).T, x.shape))

    return z


def kde_1d(xarr, xsigma, ext, grid_dens):
    '''
    Take an array of x data with their errors, create a grid of points in x
    and return the 1D KDE density map.
    '''

    # Grid density (number of points).
    gd_c = complex(0, grid_dens)

    # Define grid of points in x where the KDE will be evaluated.
    x = np.mgrid[ext[0]:ext[1]:gd_c]
    positions = np.vstack([x.ravel()])

    # Evaluate KDE in x grid.
    kde_grid = []
    for p in zip(*positions):
        kde_grid.append(kde_val(p, xarr, xsigma, [], [], '1'))

    # Re-shape values for plotting.
    z = np.array(kde_grid)

    return positions[0], z


# def solve_gaussian(val, data_array, sigma_array):
#     return (1. / sigma_array) * np.exp(
#         - (val - data_array) * (val - data_array) /
#         (2 * sigma_array * sigma_array))


# def solve_kde(xlist, data_array, sigma_array):
#     kde_array = np.array([])
#     print np.ndim(kde_array)
#     for xx in xlist:
#         single_kde = solve_gaussian(xx, data_array, sigma_array)
#         if np.ndim(kde_array) == 3:
#             kde_array = np.concatenate(
#                 (kde_array, single_kde[np.newaxis, :, :]), axis=0)
#         else:
#             kde_array = np.dstack(single_kde)
#     return kde_array

# N = 300
# data_array = np.array([np.random.uniform(0., 10., N) for _ in range(2)])
# sigma_array = np.array([np.random.uniform(0., 0.5, N) for _ in range(2)])
# xlist = np.linspace(0, 1, 101)  # Adjust as needed
# kde_array = solve_kde(xlist, data_array, sigma_array)
# kde_vector = np.sum(np.sum(kde_array, axis=2), axis=1)
# mode_guess = xlist[np.argmax(kde_vector)]
