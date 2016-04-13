
import numpy as np
from scipy import stats
from kde_map import kde_1d


def non_weight_linear_fit(x, y):
    """
    Performs a *non weighted* linear fit to data.
    """
    fit = np.polyfit(x, y, 1)
    fit_nw = np.poly1d(fit)
    return fit_nw


def weight_linear_fit(xdata, ydata, ysigma=None):
    """
    http://bulldog2.redlands.edu/facultyfolder/deweerd/tutorials/fitting.txt
    http://nbviewer.ipython.org/url/bulldog2.redlands.edu/facultyfolder/
    deweerd/tutorials/LinearRegression.txt

    Performs a linear fit to data, *weighted* by errors in y.

    Parameters
    ----------
    xdata : An array of length N.
    ydata : An array of length N.
    sigma : None or an array of length N,
        If provided, it is the standard-deviation of ydata.
        This vector, if given, will be used as weights in the fit.

    Returns
    -------
    a, b   : Optimal parameter of linear fit (y = a*x + b)
    sa, sb : Uncertainties of the parameters
    """

    if ysigma is None:
        w = np.ones(len(ydata))  # Each point is equally weighted.
    else:
        w = 1.0/(ysigma**2)

    sw = sum(w)
    wx = w*xdata  # this product gets used to calculate swxy and swx2
    swx = sum(wx)
    swy = sum(w*ydata)
    swxy = sum(wx*ydata)
    swx2 = sum(wx*xdata)

    a = (sw*swxy - swx * swy)/(sw * swx2 - swx * swx)
    b = (swy*swx2 - swx*swxy)/(sw * swx2 - swx * swx)
    sa = np.sqrt(sw / (sw * swx2 - swx * swx))
    sb = np.sqrt(swx2 / (sw * swx2 - swx * swx))

    if ysigma is None:
        chi2 = sum(((a*xdata + b)-ydata)**2)
    else:
        chi2 = sum((((a*xdata + b)-ydata)/ysigma)**2)
    dof = len(ydata) - 2
    rchi2 = chi2/dof
    # print 'results of linear_fit:'
    # print '   chi squared = ', chi2
    # print '   degrees of freedom = ', dof
    # print '   reduced chi squared = ', rchi2

    return a, b, sa, sb, rchi2, dof


def confband(xd, yd, a, b, conf=0.95):
    """
    http://astropython.blogspot.com.ar/2011/12/
    calculating-and-plotting-prediction.html
    https://gist.github.com/rsnemmen/f2c03beb391db809c90f
    https://gist.github.com/rsnemmen/0eb32832c657c19e4d39

    Calculates the confidence band of the linear regression model at the
    desired confidence level, using analytical methods. The 2sigma confidence
    interval is 95% sure to contain the best-fit regression line. This is not
    the same as saying it will contain 95% of the data points.

    Arguments:
    - xd, yd:
        data arrays
    - a, b:
        linear fit parameters as in y = ax+b
    - conf:
        desired confidence level, by default 0.95 (2 sigma)

    Returns:
    Sequence(lcb, ucb, x) with the arrays holding the lower and upper
    confidence bands corresponding to the[input] x array.

    Usage:
    >> > lcb, ucb, x = nemmen.confband(all.kp, all.lg, a, b, conf=0.95)
    calculates the confidence bands for the given input arrays
    >> > pylab.fill_between(x, lcb, ucb, alpha=0.3, facecolor='gray')
    plots a shaded area containing the confidence band

    References:
    1. http://en.wikipedia.org/wiki/Simple_linear_regression, see Section
       Confidence intervals
    2. http://www.weibull.com/DOEWeb/
       confidence_intervals_in_simple_linear_regression.htm

    Author:
        Rodrigo Nemmen
    v1 Dec. 2011
    v2 Jun. 2012:
        corrected bug in computing dy
    """

    alpha = 1. - conf   # significance
    n = xd.size   # data sample size
    x = np.linspace(xd.min(), xd.max(), 100)
    # Predicted values (best-fit model)
    y = a*x+b

    # Auxiliary definitions

    # Std. deviation of an individual measurement (Bevington, eq. 6.15)
    N = np.size(xd)
    sd = 1. / (N - 2.) * np.sum((yd - a * xd - b) ** 2)
    sd = np.sqrt(sd)

    sxd = np.sum((xd - xd.mean()) ** 2)
    sx = (x - xd.mean()) ** 2  # array

    # Quantile of Student's t distribution for p=1-alpha/2
    q = stats.t.ppf(1. - alpha / 2., n - 2)

    # Confidence band
    dy = q * sd * np.sqrt(1. / n + sx / sxd)
    ucb = y + dy    # Upper confidence band
    lcb = y - dy    # Lower confidence band

    return lcb, ucb, x


# def mean_confidence_interval(data, confidence=0.95):
#     """
#     http://stackoverflow.com/a/15034143/1391441
#     """
#     a = 1.0*np.array(data)
#     n = len(a)
#     m, se = np.mean(a), stats.sem(a)
#     h = se * stats.t._ppf((1+confidence)/2., n-1)
#     return m, [m-h, m+h]


def monte_carlo_conf_bands(xarr, xsigma, x_rang, grid_dens):
    """
    Calculate 95% confidence bands for a 1D KDE curve.

    Good explanation: http://www.graphpad.com/guides/prism/6/statistics/
    index.htm?stat_more_about_confidence_interval.htm
    """

    # Draw N_mc KDE values for each x point in the grid, perturbing the input
    # data (xarr) through a Gaussian distribution.
    N_mc, y_vals = 50, []
    for _ in xrange(N_mc):
        xarr_ran = np.random.normal(xarr, xsigma)
        x, y = kde_1d(xarr_ran, np.array(xsigma), x_rang, grid_dens)
        y_vals.append(y)

    # Obtain the 95% confidence interval for each x point in the 1D grid.
    kde_vals = zip(*y_vals)
    low_conf_curve, up_conf_curve = [], []
    for x_pos in kde_vals:
        mu, sigma, N = np.mean(x_pos), np.std(x_pos), len(x_pos)
        conf_int = stats.norm.interval(0.95, loc=mu, scale=sigma/np.sqrt(N))
        # mu, conf_int = mean_confidence_interval(x_pos)
        low_conf_curve.append(conf_int[0])
        up_conf_curve.append(conf_int[1])

    # import matplotlib.pyplot as plt
    # x, y = kde_1d(xarr, np.array(xsigma), [x_rang[0], x_rang[1]],
    #               grid_dens)
    # plt.plot(x, y, color='g')
    # plt.plot(x, low_conf_curve, color='r')
    # plt.plot(x, up_conf_curve, color='b')
    # plt.show()

    return low_conf_curve, up_conf_curve


# if __name__ == "__main__":
#     xarr = np.random.uniform(0., 10., 300)
#     xsigma = np.random.uniform(0., 0.5, 300)
#     x_rang = [0., 10.]
#     grid_dens = 20
#     monte_carlo_conf_bands(xarr, xsigma, x_rang, grid_dens)
    # main()
