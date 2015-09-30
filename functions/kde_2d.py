

import numpy as np


def kernel(p0, p, s):
    '''
    Gaussian kernel.
    '''
    if s > 0.:
        val = (1. / np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((p0 - p) / s) ** 2)
    else:
        # Replace 0 error with very small value.
        s = 0.000001
        val = (1. / np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((p0 - p) / s) ** 2)

    return val


def kde_val(point, xarr, xsigma, yarr, ysigma):
    '''
    Returns the value of the 2D Kernel Density Estimator (KDE) evaluated in
    the point passed.
    '''
    # Evaluate KDE in 'point'.
    kde_p = 0.
    for x, sx, y, sy in zip(xarr, xsigma, yarr, ysigma):
        kde_p = kde_p + kernel(point[0], x, sx) * kernel(point[1], y, sy)
    # Normalize.
    norm_kde = kde_p / len(xarr)

    return norm_kde


def kde_map(xarr, xsigma, yarr, ysigma, ext, grid_dens):
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
        kde_grid.append(kde_val(p, xarr, xsigma, yarr, ysigma))
        print 'p, kde', p, kde_val(p, xarr, xsigma, yarr, ysigma)

    # Re-shape values for plotting.
    for el in np.reshape(np.array(kde_grid).T, x.shape):
        print 'kde_no_rot', el
    z = np.rot90(np.reshape(np.array(kde_grid).T, x.shape))
    for el in z:
        print 'kde_rot', el

    return z


# def measure(n):
#     "Return two coupled measurements."
#     m1 = np.random.normal(size=n)
#     m2 = np.random.normal(scale=0.5, size=n)
#     return m1 + m2, m1 - m2


# import matplotlib.pyplot as plt
# # Create x,y data an its errors (sigmas).
# N = 250
# # Data.
# xarr, yarr = measure(N)
# # Errors.
# xsigma, ysigma = np.random.uniform(0., 1., N), np.random.uniform(0., 1., N)


# # Plot.
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.imshow(z, cmap=plt.cm.gist_earth_r, extent=[xmin, xmax, ymin, ymax])
# #ax.set_aspect('auto')
# #plt.errorbar(xarr, yarr, xerr=xsigma, yerr=ysigma, fmt=None, color='r')
# plt.scatter(xarr, yarr, marker='.', color='r', s=5)
# ax.set_xlim([xmin, xmax])
# ax.set_ylim([ymin, ymax])
# plt.show()
