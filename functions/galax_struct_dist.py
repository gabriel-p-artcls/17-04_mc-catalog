
from astropy.coordinates import Angle, SkyCoord, Distance
from astropy import units as u
import numpy as np
import deproj_dist as dd
from MCs_data import MCs_data
import scipy.interpolate
import scipy.optimize as optimize
from tabulate import tabulate


def inc_PA_grid(N, inc_rang, pa_rang):
    '''
    Generate a grid of inclination and position angles for the ranges and
    size of grid set.
    '''
    inc_lst = np.linspace(inc_rang[0], inc_rang[1], N)
    pa_lst = np.linspace(pa_rang[0], pa_rang[1], N)

    return inc_lst, pa_lst


def gal_data(ra, dec, dist_cent, e_d_cent, aarr, darr, dsigma, j):
    """
    Return data for the selected galaxy: j=0 --> SMC ; j=1 --> LMC.
    """
    # Equatorial coordinates for clusters in this galaxy.
    ra_g, dec_g = ra[j], dec[j]
    # 3D distances from clusters to center of galaxies obtained by ASteCA
    # + astropy, in Kpc.
    d_d_g, e_dd_g = np.asarray(dist_cent[j]) / 1000.,\
        np.asarray(e_d_cent[j]) / 1000.
    # ASteCA ages for clusters that belong to this galaxy.
    age_g = aarr[j][0]
    # ASteCA distance moduli and their errors.
    dm_g, e_dm_g = darr[j][0], dsigma[j][0]

    # Retrieve center coordinates and distance (in parsecs) to galaxy.
    gal_cent, gal_dist, e_dm_dist = MCs_data(j)

    return ra_g, dec_g, d_d_g, e_dd_g, age_g, dm_g, e_dm_g, gal_cent, gal_dist


def get_rho_phi(ra, dec, gal_cent):
    '''
    Obtain projected angular distance from center of galaxy to cluster, and
    position angle of cluster.
    '''
    coords = SkyCoord(zip(*[ra, dec]), unit=(u.deg, u.deg))
    # Obtain angular projected distance and position angle for the cluster.
    rho, Phi, phi = dd.rho_phi(coords, gal_cent)

    return rho, phi


def dist_filter(r_min, ra_g, dec_g, age_g, d_d_g, e_dd_g, dm_g, e_dm_g,
                gal_cent):
    """
    Filter clusters based on their projected angular distances 'rho'.

    Values are used as: (r_min, r_max]
    """

    # Obtain angular projected distance and position angle for the
    # clusters in the galaxy.
    rho_g, phi_g = get_rho_phi(ra_g, dec_g, gal_cent)

    ra_f, dec_f, age_f, d_d_f, e_dd_f, dm_f, e_dm_f, rho_f, phi_f =\
        [], [], [], [], [], [], [], [], []
    for i, d in enumerate(d_d_g):
        # if r_min < rho_g[i].degree <= r_max:
        if r_min < rho_g[i].degree:
            ra_f.append(ra_g[i])
            dec_f.append(dec_g[i])
            age_f.append(age_g[i])
            d_d_f.append(d_d_g[i])
            e_dd_f.append(e_dd_g[i])
            dm_f.append(dm_g[i])
            e_dm_f.append(e_dm_g[i])
            rho_f.append(rho_g[i].degree)
            phi_f.append(phi_g[i].degree)

    rho_f = Angle(rho_f, unit=u.deg)
    phi_f = Angle(phi_f, unit=u.deg)

    return ra_f, dec_f, age_f, d_d_f, e_dd_f, dm_f, e_dm_f, rho_f, phi_f


def get_deproj_dist(gal_dist, inc, pa, rho, phi):
    '''
    Obtain deprojected distance between cluster and the center of the MC
    in kpc.
    '''
    dist_kpc = dd.deproj_dist(pa, inc, gal_dist, rho, phi).value

    return dist_kpc


def i_PA_dist_vals(rho, phi, inc_lst, pa_lst, gal_dist):
    '''
    Calculate deprojected distance values for each (i, PA) point in the
    defined grid.
    '''

    # Create empty list with correct shape.
    dep_dist_i_PA_vals = [[[] for _ in inc_lst] for _ in pa_lst]
    for i, inc in enumerate(inc_lst):
        for j, pa in enumerate(pa_lst):
            # Assign 'degrees' units before passing.
            inc, pa = Angle(inc, unit=u.degree), Angle(pa, unit=u.degree)

            # Obtain deprojected distances for all the clusters, in kpc,
            # using the values of inclination and position angles passed.
            dep_dist_kpc = get_deproj_dist(gal_dist, inc, pa, rho, phi)

            # Store deprojected distance values.
            dep_dist_i_PA_vals[i][j] = dep_dist_kpc

    return dep_dist_i_PA_vals


def plane_equation(inc_lst, pa_lst):
    '''
    Calculate the equation representing the inclined plane for given
    inclination and position angle values:

    a*x + b*y + c*z + d = 0

    The coefficients for this x',y' plane equation is obtained via:

    z' = 0

    where x',y',z' are taken from van der Marel & Cioni (2001).
    '''

    a_lst, b_lst, c_lst = [], [], []
    for i, inc in enumerate(inc_lst):
        # Assign 'degrees' units before passing.
        inc = Angle(inc, unit=u.degree)
        for j, pa in enumerate(pa_lst):
            # Assign 'degrees' units before passing.
            pa = Angle(pa, unit=u.degree)

            # Convert PA (N-->E) to theta (W-->E).
            theta = pa + Angle('90.d')

            # Obtain coefficients for the inclined x',y' plane.
            a_lst.append(-1.*np.sin(theta.radian) * np.sin(inc.radian))
            b_lst.append(np.cos(theta.radian) * np.sin(inc.radian))
            c_lst.append(np.cos(inc.radian))
            # This equals 1, so don't pass.
            # sq = np.sqrt(a**2 + b**2 + c**2)

    plane_abc = [a_lst, b_lst, c_lst]

    return plane_abc


def ccc(l1, l2):
    '''
    Concordance correlation coefficient.
    See: https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    '''
    ccc_val = 2 * np.cov(l1, l2)[0, 1] / (np.var(l1) + np.var(l2) +
                                          (np.mean(l1) - np.mean(l2)) ** 2)

    return ccc_val


def ccc_map(dep_dist_i_PA_vals, rand_dist_kpc, N_grid):
    '''
    Obtain CCC value comparing the deprojected distances calculated for
    each plane defined by an inclination and position angle, with the values
    for the same distance obtained via the distance moduli given by ASteCA.
    '''
    ccc_lst = []
    for i in range(N_grid):
        for j in range(N_grid):
            # Retrieve deprojected distance values already calculated, using
            # the van der Marel & Cioni (2001) equations.
            dep_dist_kpc = dep_dist_i_PA_vals[i][j]

            # Calculate CCC between deprojected distances using the ASteca
            # distance moduli + astropy, and those obtained via the van
            # der Marel & Cioni (2001) equations.
            c = ccc(dep_dist_kpc, rand_dist_kpc)

            # Store CCC values.
            ccc_lst.append(c)

    # Pass as numpy array.
    ccc_lst = np.asarray(ccc_lst)

    return ccc_lst


def interp_dens_map(inc_lst, pa_lst, xi, yi, z):
    '''
    Interpolate density map of z values into a finer grid.
    '''
    # Define interpolating function.
    N = len(inc_lst)
    z = z.reshape(N, N)
    rbs = scipy.interpolate.RectBivariateSpline(inc_lst, pa_lst, z)

    # Get z values on finer grid.
    zi = rbs(xi, yi)

    return zi


def xyz_coords(rho, phi, D_0, r_dist):
    '''
    Calculate coordinates in the (x,y,z) system. The x,y plane is parallel
    the the plane of the sky, and the z axis points towards the observer
    (see van der Marel & Cioni (2001) and van der Marel et al. (2002)).
    '''
    # Convert distance moduli to Kpc.
    d_kpc = Distance((10**((r_dist + 5.)*0.2)) / 1000., unit=u.kpc)

    # Obtain coordinates.
    x = d_kpc * np.sin(rho.radian) * np.cos(phi.radian)
    y = d_kpc * np.sin(rho.radian) * np.sin(phi.radian)
    z = D_0.kpc*u.kpc - d_kpc*np.cos(rho.radian)

    return x, y, z


def fix_plane_perp_dist(plane_abc, x, y, z):
    '''
    Calculate the distance to each inclined plane for all the clusters in
    the galaxy being analysed. Plane equation in the form:

    a*x + b*y + c*z + d = 0

    Pass the averaged sum of the absolute values of each distance, for each
    inclination and position angle values.

    http://mathworld.wolfram.com/Point-PlaneDistance.html
    '''
    # Unpack lists of inclined planes coefficients.
    a_lst, b_lst, c_lst = plane_abc
    # Calculate the averaged sum of (positive) distances to each inclined
    # plane.
    pl_dists_kpc = np.sum(abs(np.outer(a_lst, x) + np.outer(b_lst, y) +
                              np.outer(c_lst, z)), axis=1)/len(x)

    return pl_dists_kpc


def perp_error(params, xyz):
    """
    Return the averaged sum of the absolute values for the perpendicular
    distance of the points in 'xyz', to the plane defined by the
    coefficients 'a,b,c,d'.
    """
    a, b, c, d = params
    x, y, z = xyz
    length = np.sqrt(a**2+b**2+c**2)
    return (np.abs(a*x + b*y + c*z + d).sum()/length)/len(x)


def minimize_perp_distance(x, y, z, N_min):
    """
    Find coefficients of best plane fit to given points. Plane equation is
    in the form:

    a*x + b*t + c*z + d = 0

    and the minimization is done for the perpendicular distances to the plane.

    Source: http://stackoverflow.com/a/35118683/1391441
    """
    def unit_length(params):
        a, b, c, d = params
        return a**2 + b**2 + c**2 - 1

    # Remove units from the x,y,z coordinates of the points passed.
    x, y, z = x.value, y.value, z.value

    # Random initial guess for the coefficients.
    initial_guess = np.random.uniform(-10., 10., 4)

    # # Similar as the block below, but using a simpler random sampling
    # # algorithm. In several tests both methods gave the same coefficients,
    # # with this one occasionally failing.
    # # This algorithm is 4-5 times faster than the Basin-Hopping one.
    # # Leave it here to remember I tried.

    # results = []
    # for _ in xrange(50):

    #     # Generate a new random initial guess every 10 runs.
    #     if _ in [10, 20, 30, 40]:
    #         initial_guess = np.random.uniform(-10., 10., 4)

    #     # Constrain the vector perpendicular to the plane be of unit length
    #     cons = ({'type': 'eq', 'fun': unit_length})
    #     sol = optimize.minimize(perp_error, initial_guess, args=[x, y, z],
    #                             constraints=cons)

    #     # Extract minimized coefficients.
    #     abcd = list(sol.x)
    #     # Update initial guess.
    #     initial_guess = abcd
    #     # Save coefficients and the summed abs distance to the plane.
    #     results.append([abcd, perp_error(list(sol.x), [x, y, z])])

    # # Select the coefficients with the minimum value for the summed absolute
    # # distance to plane.
    # abcd = min(results, key=lambda x: x[1])[0]

    # Use Basin-Hopping to obtain the best fir coefficients.
    cons = ({'type': 'eq', 'fun': unit_length})
    min_kwargs = {"constraints": cons, "args": [x, y, z]}
    sol = optimize.basinhopping(
        perp_error, initial_guess, minimizer_kwargs=min_kwargs, niter=N_min)
    abcd = list(sol.x)

    return abcd


def angle_betw_planes(plane_abcd):
    """
    The dihedral angle (inclination) is the angle between two planes

    http://mathworld.wolfram.com/DihedralAngle.html

    http://stackoverflow.com/q/2827393/1391441
    """

    # Normal vector to the a,b,c,d plane.
    a, b, c, d = plane_abcd

    # Counter clockwise rotation angle around the z axis.
    # 0. <= theta <= 180.
    if d != 0.:
        x_int, y_int = -d/a, -d/b
        if y_int > 0.:
            t = 180. - np.arctan2(y_int, x_int)*180./np.pi
        elif y_int < 0.:
            t = abs(np.arctan2(y_int, x_int)*180./np.pi)
    else:
        x_int, y_int = 1., -a/b
        if y_int > 0.:
            t = np.arctan2(y_int, x_int)*180./np.pi
        elif y_int < 0.:
            t = np.arctan2(y_int, x_int)*180./np.pi + 180.
    theta = Angle(t, unit='degree')

    # Set theta to correct range [90., 270.] (given that the position angle's
    # range is [0., 180.])
    if theta.degree < 90.:
        n = int((90. - theta.degree)/180.) + 1
        theta = theta + n*Angle('180.d')
    # Pass position angle (PA) instead of theta, to match the other methods.
    # We use the theta = PA + 90 convention.
    PA = theta - Angle('90.d')

    # Clockwise rotation angle around the x' axis.
    # Obtain the inclination with the correct sign.
    if c != 0.:
        # Normal vector to inclined plane.
        v1 = [a/c, b/c, 1.]
        # theta = 90.d
        if b == 0.:
            if a/c > 0.:
                sign = 'neg'
            elif a/c < 0.:
                sign = 'pos'
            else:
                sign = 'no'
        else:
            # 90. < theta < 270.
            if b/c > 0.:
                sign = 'neg'
            elif b/c < 0.:
                sign = 'pos'

        # Set correct sign for inclination angle.
        if sign == 'pos':
            # Vector normal to the z=0 plane, ie: the x,y plane.
            # This results in a positive i value.
            v2 = [0, 0, 1]
        elif sign == 'neg':
            # This results in a negative i value.
            v2 = [0, 0, -1]
        else:
            # This means i=0.
            v2 = v1

        i = np.arccos(np.dot(v1, v2) / np.sqrt(np.dot(v1, v1)*np.dot(v2, v2)))
        inc = Angle(i*u.radian, unit='degree')

    else:
        inc = Angle(90., unit='degree')
    # 0. <= inc <= 180.

    # Set inclination angles to correct ranges [-90, 90]
    if inc.degree > 90.:
        inc = Angle(i*u.radian, unit='degree') - Angle('180.d')

    return inc.degree, PA.degree


def best_fit_angles(method, best_angles_pars):
    '''
    Obtain inclination and position angles associated with the maximum CCC
    value in this map.
    '''
    if method == 'deproj_dists':
        xi, yi, zi = best_angles_pars
        # Obtain index of maximum density value.
        max_ccc_idx = np.unravel_index(zi.argmax(), zi.shape)
        # Calculate "best" angle values.
        inc_b, pa_b = xi[max_ccc_idx[0]], yi[max_ccc_idx[1]]
        # Max CCC value in map.
        # z_b = np.amax(zi)
    elif method == 'perp_d_fix_plane':
        xi, yi, zi = best_angles_pars
        # Obtain index of minimum sum of absolute distances to plane value.
        min_d_idx = np.unravel_index(zi.argmin(), zi.shape)
        # Calculate "best" angle values.
        inc_b, pa_b = xi[min_d_idx[0]], yi[min_d_idx[1]]
        # Min sum of abs(distances 2 plane) value in map.
        # z_b = np.amin(zi)
    elif method == 'perp_d_free_plane':
        plane_abcd = best_angles_pars
        inc_b, pa_b = angle_betw_planes(plane_abcd)

    return inc_b, pa_b


def draw_rand_dep_dist(d_g, e_d_g):
    '''
    Take each 3D distance between a cluster and the center of the galaxy,
    and normally draw a random value using its error as the standard deviation.
    '''
    rand_dist = []
    for mu, std in zip(*[d_g, e_d_g]):
        r_dist = np.random.normal(mu, std)
        # Don't store negative distances.
        rand_dist.append(max(0., r_dist))

    return rand_dist


def draw_rand_dist_mod(dm_g, e_dm_g):
    """
    Draw random distance moduli assuming a normal distribution of the errors.
    """
    r_dist = np.random.normal(np.asarray(dm_g), np.asarray(e_dm_g))

    return r_dist


def monte_carlo_errors(N_maps, method, params):
    """
    Obtain N_maps angles-CCC density maps, by randomly sampling
    the distance to each cluster before calculating the CCC.
    """

    # Unpack params.
    if method == 'deproj_dists':
        N_grid, inc_lst, pa_lst, xi, yi, d_f, e_d_f, dep_dist_i_PA_vals =\
            params
    elif method == 'perp_d_fix_plane':
        inc_lst, pa_lst, xi, yi, dm_f, e_dm_f, rho_f, phi_f, gal_dist,\
            plane_abc = params
    elif method == 'perp_d_free_plane':
        dm_f, e_dm_f, rho_f, phi_f, gal_dist, N_min = params

    inc_pa_mcarlo = []
    milestones = [0.1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for _ in range(N_maps):

        if method == 'deproj_dists':
            # Draw random deprojected distances (in Kpc), obtained via the
            # ASteCA distance moduli + astropy.
            rand_dist = draw_rand_dep_dist(d_f, e_d_f)
            # Obtain density map of CCC (z) values.
            z = ccc_map(dep_dist_i_PA_vals, rand_dist, N_grid)

            # Obtain (finer) interpolated angles-CCC density map.
            # Rows in zi correspond to inclination values.
            # Columns correspond to position angle values.
            zi = interp_dens_map(inc_lst, pa_lst, xi, yi, z)

            # Store values of inclination and position angles for the
            # interpolated map.
            best_angles_pars = [xi, yi, zi]

        elif method == 'perp_d_fix_plane':
            # Draw random distances moduli, obtained via ASteCA.
            rand_dist = draw_rand_dist_mod(dm_f, e_dm_f)
            # Positions in the (x,y,z) system.
            x, y, z = xyz_coords(rho_f, phi_f, gal_dist, rand_dist)
            # Obtain density map (z), composed of the sum of the absolute
            # values of the distances to each plane.
            z = fix_plane_perp_dist(plane_abc, x, y, z)
            zi = interp_dens_map(inc_lst, pa_lst, xi, yi, z)
            best_angles_pars = [xi, yi, zi]

        elif method == 'perp_d_free_plane':
            # Draw random distances moduli, obtained via ASteCA.
            rand_dist = draw_rand_dist_mod(dm_f, e_dm_f)
            # Obtain randomly drawn coords in the (x, y, z) sustem.
            x, y, z = xyz_coords(rho_f, phi_f, gal_dist, rand_dist)
            # Store params used to obtain the Monte Carlo errors.
            best_angles_pars = minimize_perp_distance(x, y, z, N_min)

        inc_b, pa_b = best_fit_angles(method, best_angles_pars)

        inc_pa_mcarlo.append([inc_b, pa_b])

        # Print percentage done.
        percentage_complete = (100. * (_ + 1) / N_maps)
        while len(milestones) > 0 and \
                percentage_complete >= milestones[0]:
            # print " {:>3}% done".format(milestones[0])
            # Remove that milestone from the list.
            milestones = milestones[1:]

    return inc_pa_mcarlo


def eigsorted(cov):
    '''
    Eigenvalues and eigenvectors of the covariance matrix.
    '''
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:, order]


def cov_ellipse(points, nstd=1):
    """
    Generate an `nstd` sigma ellipse based on the mean and covariance of a
    point "cloud".

    source: http://stackoverflow.com/a/12321306/1391441

    Parameters
    ----------
        points : An Nx2 array of the data points.
        nstd : The radius of the ellipse in numbers of standard deviations.
               Defaults to 1 standard deviations.
    """
    points = np.asarray(points)

    # Location of the center of the ellipse.
    mean_pos = points.mean(axis=0)

    # The 2x2 covariance matrix to base the ellipse on.
    cov = np.cov(points, rowvar=False)

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)

    return mean_pos, width, height, theta


def ccc_sum_d_for_best_fit(gal_dist, rho_f, phi_f, d_d_f, cl_x, cl_y, cl_z,
                           best_angles_pars, inc_b, pa_b, method):
    """
    For the best fit angles obtained with this method, calculate:

    1. The CCC coefficient comparing the deprojected distances on the fitted
    plane with the ASteCA+astropy distances.
    2. The sum of absolute values of the perpendicular distances to the plane.

    Return also the deprojected distances on the fitted plane, for plotting.
    """
    # Deprojected distances obtained using the best-fit angles.
    # THIS ASSUMES THAT THE INCLINED PLANE IS OBTAINED BY ROTATING TWICE
    # THE (x,y,z) SYSTEM. THE VALUES GIVEN FOR THE ANGLES OBTAINED WITH THE
    # FREE PLANE BEST FIT METHOD WILL THUS REPRESENT AN APPROXIMATION OF
    # THE ACTUAL VALUES IF THE (x',y') PLANE INTERSECTED THE (x,y) PLANE
    # THROUGH THE ORIGIN (ie: IF  d=0 IN THE PLANE'S EQUATION)
    dep_dist_kpc = get_deproj_dist(
        gal_dist, Angle(inc_b, unit=u.degree),
        Angle(pa_b, unit=u.degree), rho_f, phi_f)
    # CCC value for the best fit angles.
    ccc_b = ccc(dep_dist_kpc, d_d_f)

    xyz = [cl_x.value, cl_y.value, cl_z.value]
    if method in ['deproj_dists', 'perp_d_fix_plane']:
        # Convert best fit PA to theta.
        theta = pa_b + 90.
        # Plane coefficients according to Eq (6) in vdM&C01 for z'=0.
        a = -1.*np.sin(np.deg2rad(theta))*np.sin(np.deg2rad(inc_b))
        b = np.cos(np.deg2rad(theta))*np.sin(np.deg2rad(inc_b))
        c = np.cos(np.deg2rad(inc_b))
        d = 0.
        abcd = [a, b, c, d]
        sum_d_b = perp_error(abcd, xyz)
    elif method == 'perp_d_free_plane':
        sum_d_b = perp_error(best_angles_pars, xyz)

    return dep_dist_kpc, ccc_b, sum_d_b


def gsd(in_params):
    '''
    Calculate the best match for the inclination and position angles of the
    MCs, based on the distance assigned to each cluster by ASteCA.
    '''
    ra, dec, aarr, dist_cent, e_d_cent, darr, dsigma = [
        in_params[_] for _ in ['ra', 'dec', 'aarr', 'dist_cent',
                               'e_d_cent', 'darr', 'dsigma']]

    # Define ranges for the grid of inclination and position angles.
    inc_rang, pa_rang = [-89., 89.], [1., 179.]

    # Grid limits (for plotting).
    xmin, xmax = inc_rang[0] - 0.1, inc_rang[1] + 0.1
    ymin, ymax = pa_rang[0] - 0.1, pa_rang[1] + 0.1

    # Obtain grid of inclination and position angles of size: N x N
    N_grid = 20  # FIXME
    inc_lst, pa_lst = inc_PA_grid(N_grid, inc_rang, pa_rang)

    # Obtain *denser/finer* grid of inclination and position angles.
    N_f = 200
    xi, yi = inc_PA_grid(N_f, inc_rang, pa_rang)

    # Number of iterations for the minimization algorithm that searches for the
    # best plane fit to the clusters, with no constrains imposed.
    N_min = 5  # FIXME

    # Number of Monte Carlo runs, where the distance to each
    # cluster is randomly drawn from a normal probability distribution.
    # This is used to assign errors to the best fit angle values.
    N_maps = 5  # FIXME

    # Store parameters needed for plotting the density maps and 1:1 diagonal
    # plots of deprojected distance in Kpc, as well as the r_min plots.
    gal_str_pars, rho_plot_pars, table = [[], []], [[], []], []
    for j in [0, 1]:  # SMC, LMC = 0, 1
        print 'Galaxy:', j

        # Retrieve data for this galaxy.
        ra_g, dec_g, d_d_g, e_dd_g, age_g, dm_g, e_dm_g, gal_cent, gal_dist =\
            gal_data(ra, dec, dist_cent, e_d_cent, aarr, darr, dsigma, j)

        # Input minimum projected angular distance values to use in filter.
        # The value ia used as: (r_min...]
        rho_lst = [0., 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4.]  # FIXME
        # rho_lst = [0., 0.5, 1.]  # DELETE
        # Select index of r_min value to plot.
        r_idx_save = 2
        for r_idx, r_min in enumerate(rho_lst):

            # Filter clusters by distance to center of galaxy.
            ra_f, dec_f, age_f, d_d_f, e_dd_f, dm_f, e_dm_f, rho_f, phi_f =\
                dist_filter(r_min, ra_g, dec_g, age_g, d_d_g, e_dd_g,
                            dm_g, e_dm_g, gal_cent)

            # Calculate deprojected distances for all clusters in this galaxy,
            # for all inclination and position angles defined in the grid.
            # These values only depend on the coordinates of the clusters, the
            # rotation angles that define each inclined plane, and the distance
            # to the galaxy.
            dep_dist_i_PA_vals = i_PA_dist_vals(rho_f, phi_f, inc_lst, pa_lst,
                                                gal_dist)

            # Obtain coefficients for the plane equation defined by each
            # combination of rotation (inc, PA) angles.
            plane_abc = plane_equation(inc_lst, pa_lst)

            # Obtain coordinates of filtered clusters in the (x,y,z) system.
            # Used by two of the methods below.
            cl_x, cl_y, cl_z = xyz_coords(rho_f, phi_f, gal_dist,
                                          np.asarray(dm_f))

            # Run once for each method defined.
            inc_best, pa_best, mc_inc_std, mc_pa_std, in_mcarlo, pa_mcarlo,\
                ccc_sum_d_best = [], [], [], [], [], [], []
            for method in ['deproj_dists', 'perp_d_fix_plane',
                           'perp_d_free_plane']:
                if method == 'deproj_dists':
                    # Store params used to obtain the Monte Carlo errors.
                    params = [N_grid, inc_lst, pa_lst, xi, yi, d_d_f, e_dd_f,
                              dep_dist_i_PA_vals]

                    # Store CCC density map obtained using the distance values
                    # with no random sampling.
                    ccc_vals = ccc_map(dep_dist_i_PA_vals, d_d_f, N_grid)
                    # Interpolate density map into finer/denser grid.
                    dens_vals_interp = interp_dens_map(inc_lst, pa_lst, xi, yi,
                                                       ccc_vals)
                    # Data to pass to extract best fit angles.
                    best_angles_pars = [xi, yi, dens_vals_interp]

                elif method == 'perp_d_fix_plane':
                    # Store params used to obtain the Monte Carlo errors.
                    params = [inc_lst, pa_lst, xi, yi, dm_f, e_dm_f, rho_f,
                              phi_f, gal_dist, plane_abc]

                    # Store density map obtained using the distance values with
                    # no random sampling.
                    pl_dists_kpc = fix_plane_perp_dist(plane_abc, cl_x, cl_y,
                                                       cl_z)
                    # Interpolate density map into finer/denser grid.
                    dens_vals_interp = interp_dens_map(inc_lst, pa_lst, xi, yi,
                                                       pl_dists_kpc)
                    # Data to pass to extract best fit angles.
                    best_angles_pars = [xi, yi, dens_vals_interp]

                elif method == 'perp_d_free_plane':
                    # Store params used to obtain the Monte Carlo errors.
                    params = [dm_f, e_dm_f, rho_f, phi_f, gal_dist, N_min]

                    # Obtain a,b,c,d coefficients for the best fit plane,
                    # given the set of clusters passed in the (x,y,z) system.
                    # The best fit is obtained minimizing the perpendicular
                    # distance to the plane.
                    best_angles_pars = minimize_perp_distance(cl_x, cl_y,
                                                              cl_z, N_min)

                    # Obtain density map from the 'perp_d_fix_plane' method.
                    # Used for plotting since this method doesn't generate a
                    # density map.
                    pl_dists_kpc = fix_plane_perp_dist(plane_abc, cl_x, cl_y,
                                                       cl_z)
                    dens_vals_interp = interp_dens_map(inc_lst, pa_lst, xi, yi,
                                                       pl_dists_kpc)

                # Best fit angles for the density map with no random sampling.
                inc_b, pa_b = best_fit_angles(method, best_angles_pars)
                # Save best fit angles obtained with all methods. Their
                # average is the final value for each rotation angle.
                inc_best.append(inc_b)
                pa_best.append(pa_b)

                # Retrieve the CCC value and the sum of the abs values of the
                # deprojected distances, for the distances obtained with the
                # best fit rotation angles versus the ones given by
                # ASteCA + astropy.
                dep_dist_kpc, ccc_b, sum_d_b = ccc_sum_d_for_best_fit(
                    gal_dist, rho_f, phi_f, d_d_f, cl_x, cl_y, cl_z,
                    best_angles_pars, inc_b, pa_b, method)
                ccc_sum_d_best.append([ccc_b, sum_d_b])

                # Obtain distribution of rotation angles via Monte Carlo random
                # sampling.
                inc_pa_mcarlo = monte_carlo_errors(N_maps, method, params)
                # Standard deviations for the angles for *each* method.
                mc_inc_std.append(np.std(zip(*inc_pa_mcarlo)[0]))
                mc_pa_std.append(np.std(zip(*inc_pa_mcarlo)[1]))
                # Save inclination and position angles obtained via the Monte
                # Carlo process. Combining values from *all* methods, we
                # calculate the combined standard deviation for each angle.
                in_mcarlo = in_mcarlo + list(zip(*inc_pa_mcarlo)[0])
                pa_mcarlo = pa_mcarlo + list(zip(*inc_pa_mcarlo)[1])

                # Calculate mean inclination and position angles, along with
                # the 1 sigma error ellipsoid (for plotting).
                # The theta angle rotates the ellipse counter-clockwise
                # starting from the positive x axis.
                mean_pos, width, height, theta = cov_ellipse(inc_pa_mcarlo)
                # Calculate standard deviation for inclination and position
                # angles.
                i_pa_std = np.std(np.asarray(inc_pa_mcarlo), axis=0)

                # Store parameters for density maps and 1:1 diagonal plots.
                if r_idx == r_idx_save:

                    gal_str_pars[0].append([xmin, xmax, ymin, ymax, xi, yi,
                                           dens_vals_interp, [inc_b, pa_b],
                                           i_pa_std, width, height, theta])

                    # Append dummy values at the end.
                    gal_str_pars[1].append([0.1, 7.95, 0.01, 7.95, d_d_f,
                                           dep_dist_kpc, age_f, ccc_sum_d_best,
                                           '', '', '', ''])

            # Number of clusters used in this run. Used for plotting.
            N_clust = len(ra_f)

            # # Mean and standard deviation for the rotation angles.
            # inc_std = np.std(in_mcarlo)
            # pa_std = np.std(pa_mcarlo)
            # # Pass best fit CCC values for colorbar.
            # ccc_best = zip(*ccc_sum_d_best)[0]
            # # Store parameters for plotting.
            # rho_plot_pars[j].append([r_min, N_clust, inc_best, inc_std,
            #                          pa_best, pa_std, ccc_best])

            # Mean and standard deviation for the rotation angles.
            inc_mean, inc_std = np.mean(inc_best), np.std(in_mcarlo)
            pa_mean, pa_std = np.mean(pa_best), np.std(pa_mcarlo)
            perp_sum_mean = np.mean(zip(*ccc_sum_d_best)[1])
            # Store parameters for plotting.
            rho_plot_pars[j].append([r_min, N_clust, inc_mean, inc_std,
                                     pa_mean, pa_std, perp_sum_mean])

            print 'rho min=', r_min
            print 'CCC/Perp sum=', ccc_sum_d_best
            print 'PA Best, MC std:', pa_best, mc_pa_std
            print 'Inc best, MC std:', inc_best, mc_inc_std
            print 'PA combined MC:', pa_mean, pa_std
            print 'Inc combined MC:', inc_mean, inc_std, '\n'

            table.append([
                "{:.1f}".format(r_min),
                r"${:.1f}pm{:.1f}$".format(pa_best[0], mc_pa_std[0]),
                r"${:.1f}pm{:.1f}$".format(inc_best[0], mc_inc_std[0]),
                r"${:.1f}pm{:.1f}$".format(pa_best[1], mc_pa_std[1]),
                r"${:.1f}pm{:.1f}$".format(inc_best[1], mc_inc_std[1]),
                r"${:.1f}pm{:.1f}$".format(pa_best[2], mc_pa_std[2]),
                r"${:.1f}pm{:.1f}$".format(inc_best[2], mc_inc_std[2]),
                r"${:.1f}pm{:.1f}$".format(pa_mean, pa_std),
                r"${:.1f}pm{:.1f}$".format(inc_mean, inc_std)])

        print tabulate(table, tablefmt="latex")

    return gal_str_pars, rho_plot_pars


if __name__ == "__main__":
    ra = [[15.5958333333333, 12.1375, 10.9333333333333, 17.225, 15.1416666666667, 15.0041666666667, 14.3416666666667, 13.5625, 357.245833333333, 15.1041666666667, 11.3583333333333, 16.0916666666667, 14.4583333333333, 16.8333333333333, 22.6583333333333, 9.425, 18.875, 18.2583333333333, 13.3541666666667, 24.0041666666667, 11.6833333333333, 23.75, 23.3083333333333, 17.5541666666667, 25.4291666666667, 4.60416666666667, 19.5666666666667, 25.5916666666667, 11.8, 15.2833333333333, 21.2333333333333, 11.475, 29.1833333333333, 18.0166666666667, 5.66666666666667, 14.45, 22.8833333333333, 14.4458333333333, 13.275, 25.6166666666667, 11.2166666666667, 13.1458333333333, 10.7458333333333, 13.8875, 15.9708333333333, 14.425, 6.17916666666667, 20.7, 18.2125, 27.3666666666667, 5.3625, 10.35, 23.6083333333333, 5.76666666666667, 17.0791666666667, 15.1458333333333, 27.5791666666667, 11.9583333333333, 16.0166666666667, 12.3625, 17.6958333333333, 15.2333333333333, 15.4916666666667, 11.5041666666667, 14.325, 15.2041666666667, 17.0583333333333, 14.0583333333333, 16.8791666666667, 16.7375, 13.6291666666667, 12.5875, 14.3333333333333, 15.35, 17.2583333333333, 12.325, 15.1375, 10.8875, 12.1541666666667, 11.775, 16.2583333333333, 11.7291666666667, 10.9083333333333, 12.05, 11.8541666666667, 12.0041666666667, 15.7958333333333, 17.2541666666667, 15.0583333333333], [76.4166666666667, 74.225, 82.9916666666667, 78.8541666666667, 76.1416666666667, 88.0458333333333, 72.3083333333333, 76.5083333333333, 78.8125, 87.7, 76.9458333333333, 77.3, 76.1375, 77.2208333333333, 73.9208333333333, 86.4583333333333, 76.9833333333333, 83.3333333333333, 72.7958333333333, 77.7333333333333, 86.0458333333333, 72.2791666666667, 72.5875, 77.625, 86.7166666666667, 71.8583333333333, 72.6208333333333, 78.9458333333333, 76.975, 74.725, 86.7125, 73.7625, 72.25, 94.3291666666667, 76.5375, 91.8708333333333, 74.5583333333333, 93.4833333333333, 72.1541666666667, 70.8083333333333, 73.4625, 79.7583333333333, 76.55, 82.9416666666667, 74.5416666666667, 76.9416666666667, 78.1041666666667, 74.3916666666667, 76.6416666666667, 85.9833333333333, 81.55, 74.9416666666667, 69.925, 79.5208333333333, 82.6416666666667, 74.9083333333333, 82.2083333333333, 95.3916666666667, 79.4541666666667, 67.6666666666667, 77.7958333333333, 85.375, 77.3958333333333, 76.4708333333333, 69.4125, 76.5125, 74.8083333333333, 76.6041666666667, 85.8958333333333, 82.4833333333333, 79.1666666666667, 77.7875, 78.45, 76.125, 82.925, 73.1875, 82.85, 72.7, 92.2208333333333, 93.9875, 88.8958333333333, 85.8333333333333, 74.1208333333333, 84.7833333333333, 80.05, 72.4125, 73.55, 71.5166666666667, 77.3458333333333, 77.9208333333333, 77.65, 81.3625, 74.7125, 81.1166666666667, 79.2333333333333, 81.125, 68.9083333333333, 93.6166666666667, 91.6291666666667, 74.3583333333333, 73.7541666666667, 83.0125, 86.55, 93.6708333333333, 74.9708333333333, 76.8958333333333, 79.2208333333333, 79.0708333333333, 77.9166666666667, 82.4416666666667, 77.3125, 83.2541666666667, 77.5083333333333, 83.6625, 74.5625, 80.4375, 79.6708333333333, 85.4916666666667, 71.6041666666667, 73.225, 76.9041666666667, 76.4, 86.4833333333333, 85.1083333333333, 77.6666666666667, 74.5916666666667, 77.7208333333333, 73.9666666666667, 82.3333333333333, 85.4541666666667, 77.6333333333333, 79.1125, 80.0083333333333, 71.875, 88.925, 85.6208333333333, 84.4416666666667, 86.3625, 80.8, 83.0958333333333, 77.2125, 82.9625, 76.8375, 76.6708333333333, 83.5541666666667, 82.4958333333333, 85.4125, 82.4125, 76.3541666666667, 76.6416666666667]]
    dec = [[-72.0030555555556, -73.3069444444444, -72.9766666666667, -73.2416666666667, -72.3655555555556, -72.3688888888889, -71.8913888888889, -72.2413888888889, -72.9452777777778, -71.2947222222222, -73.4813888888889, -72.8477777777778, -72.9436111111111, -73.3775, -76.0544444444444, -73.9083333333333, -71.1788888888889, -70.9627777777778, -72.1963888888889, -75.4577777777778, -72.0630555555556, -75.5566666666667, -74.1672222222222, -73.2091666666667, -71.1611111111111, -74.3186111111111, -72.0016666666667, -74.1733333333333, -73.4772222222222, -74.0736111111111, -71.1836111111111, -73.5066666666667, -74.2194444444445, -75.1975, -75.0747222222222, -73.4216666666667, -71.9527777777778, -74.3266666666667, -73.3802777777778, -71.2791666666667, -73.0019444444445, -72.1930555555556, -72.5886111111111, -74.0636111111111, -72.8261111111111, -74.4727777777778, -73.755, -75.0016666666667, -73.1194444444444, -73.7283333333333, -73.7486111111111, -72.8908333333333, -72.8744444444444, -73.6697222222222, -72.8841666666667, -71.4608333333333, -74.3561111111111, -73.4783333333333, -74.6191666666667, -73.3980555555556, -72.7919444444444, -73.1516666666667, -71.0202777777778, -73.3955555555556, -72.9327777777778, -73.3488888888889, -73.2569444444444, -74.1561111111111, -73.1197222222222, -73.235, -74.1852777777778, -73.3872222222222, -71.1702777777778, -73.2402777777778, -73.0863888888889, -73.3716666666667, -72.2583333333333, -73.4388888888889, -73.4155555555556, -73.3730555555556, -73.0427777777778, -73.4233333333333, -73.4405555555556, -73.4463888888889, -73.4580555555556, -73.4861111111111, -72.2736111111111, -73.2066666666667, -72.4583333333333], [-68.6394444444445, -68.0022222222222, -67.9716666666667, -68.6811111111111, -68.2083333333333, -71.8583333333333, -72.0566666666667, -68.0263888888889, -68.8825, -71.7077777777778, -66.7980555555556, -68.4441666666667, -67.9755555555556, -68.0836111111111, -67.7833333333333, -69.3802777777778, -67.3577777777778, -68.1522222222222, -67.5347222222222, -67.6266666666667, -69.3333333333333, -67.3416666666667, -72.8275, -68.4005555555556, -69.1897222222222, -67.6597222222222, -67.3258333333333, -69.1919444444445, -67.9288888888889, -67.8469444444445, -69.4197222222222, -67.9644444444444, -72.64, -70.0608333333333, -68.4458333333333, -72.4941666666667, -67.7683333333333, -72.5052777777778, -68.5594444444444, -73.8119444444444, -69.5719444444444, -69.0011111111111, -68.0630555555556, -68.2355555555556, -68.0602777777778, -67.8613888888889, -68.7719444444444, -65.2677777777778, -68.3630555555556, -69.1805555555556, -70.9813888888889, -69.8011111111111, -74.0172222222222, -69.1716666666667, -63.2033333333333, -69.5575, -71.6327777777778, -72.79, -68.4727777777778, -66.9569444444444, -67.6266666666667, -69.185, -67.8105555555556, -67.0494444444445, -66.1994444444444, -68.6283333333333, -67.9083333333333, -68.375, -66.2086111111111, -72.0547222222222, -70.5408333333333, -67.6825, -66.62, -68.3497222222222, -72.1461111111111, -72.5177777777778, -72.0425, -72.5766666666667, -72.3838888888889, -70.0730555555556, -62.3452777777778, -66.2622222222222, -67.6227777777778, -74.8533333333333, -68.9041666666667, -72.2480555555556, -69.8069444444444, -66.9113888888889, -67.7783333333333, -67.5655555555556, -70.4875, -73.5702777777778, -69.9577777777778, -67.7286111111111, -68.6827777777778, -67.6780555555556, -73.7316666666667, -72.6094444444445, -72.2263888888889, -67.6852777777778, -67.7141666666667, -64.2422222222222, -69.0825, -69.8019444444445, -67.9236111111111, -67.4608333333333, -69.15, -69.1541666666667, -68.7266666666667, -71.0005555555556, -67.6997222222222, -67.8491666666667, -66.6991666666667, -68.3055555555556, -68.0491666666667, -68.9172222222222, -69.0794444444444, -69.0475, -72.5683333333333, -72.1725, -68.5419444444444, -68.6286111111111, -69.2719444444444, -69.2486111111111, -68.7536111111111, -69.8030555555556, -67.4711111111111, -69.7058333333333, -70.5794444444445, -68.9208333333333, -66.94, -69.0802777777778, -69.2611111111111, -72.5883333333333, -74.3538888888889, -65.3627777777778, -74.7827777777778, -69.3452777777778, -70.7777777777778, -67.9969444444444, -67.9802777777778, -67.9911111111111, -66.8291666666667, -67.8422222222222, -67.8563888888889, -67.8788888888889, -69.2294444444444, -70.9838888888889, -68.5005555555556, -68.4272222222222]]
    dist_cent = [[1633.20060682, 2320.18925372, 2441.9621715, 1471.75739309, 2376.23276698, 769.377260812, 1217.10160803, 1804.52411396, 5941.22951677, 2131.04017704, 909.015121115, 3069.45508982, 2364.25345255, 2670.26753775, 4890.13489057, 2343.49016847, 3159.23688725, 2711.83412564, 1818.55754823, 4455.35179015, 1123.62228704, 4757.29044595, 3641.00953753, 2648.4875134, 4474.93618064, 3141.57483709, 3560.24707724, 4473.53808798, 1395.16822521, 1590.48289743, 4018.53487003, 1077.40931903, 5110.49822201, 3786.85646257, 4265.8368639, 949.628191291, 3943.07061992, 1767.67615703, 1299.74797416, 4641.567255, 2418.05685486, 897.32557739, 3039.33043949, 1461.11210457, 1901.98375256, 2165.3084445, 2985.15151754, 3722.76959311, 2000.6292342, 4885.31051319, 2625.48765461, 1076.03465615, 3317.63420596, 2482.45477321, 1670.70686978, 1704.67888631, 4679.63812859, 3032.23144391, 2116.49256805, 2422.17347074, 1542.38398274, 3014.72270399, 2369.65890352, 1867.37485958, 688.161984621, 2390.12368317, 1417.8418462, 2759.56312304, 1647.94174103, 2532.76451148, 2666.69527657, 1799.77351614, 1911.04658283, 2378.98229824, 1319.33831563, 646.650576733, 1952.83158151, 2437.45811475, 1828.61711731, 1838.06758207, 1147.07949676, 2459.7173664, 1123.54113051, 1845.94527547, 2373.7888037, 2374.78095982, 1177.87078532, 1347.286255, 724.047475304], [2094.92442356, 2592.17806506, 1926.01374032, 1393.13441938, 2375.35464572, 3301.57725704, 3191.52180512, 1846.0135317, 716.76297379, 3154.98748878, 3076.59112567, 2026.05665608, 2404.83900491, 1555.09594563, 2707.74040478, 2869.88446607, 2126.74995801, 1875.64558761, 3013.12651671, 1938.91000725, 2787.32652108, 3248.71116208, 3895.8901948, 2000.81445911, 2338.12168352, 3264.69637398, 3074.99893799, 2636.47366843, 2026.34617439, 2555.05591583, 1987.45137574, 2962.78454192, 3645.08442741, 4482.58313968, 1839.78891743, 4308.71233216, 2341.82811636, 4655.60143549, 2645.31421226, 4649.97712828, 2659.06169124, 495.958510033, 1712.84745196, 2951.12877171, 2751.29010002, 1869.62116249, 1815.09401991, 4178.76760298, 1500.54520148, 1794.75124939, 2529.48578966, 1640.54629482, 5595.45850131, 776.552314833, 5577.42226175, 2249.44973927, 3280.94442372, 5192.25922245, 2306.05901115, 4799.58913132, 1930.71867668, 1618.96956185, 2738.38874315, 2453.83795634, 4563.65453258, 1390.93463993, 2210.73606753, 2156.22099955, 3489.11690674, 2339.47934036, 1913.24925842, 1765.47926278, 2664.78715134, 1747.9825195, 2559.15521594, 3985.98208869, 2470.42710796, 3449.61772162, 4247.88799427, 4176.76382235, 6981.56554962, 3438.49595369, 2544.39520164, 4846.09325026, 859.207154012, 4248.51583614, 3334.30573069, 3728.40654495, 1759.97165254, 2178.23559988, 1194.19833024, 3858.08899466, 1864.27857663, 1561.90783843, 2706.82540924, 1747.18252909, 4836.99690823, 4641.18832985, 4140.00763706, 2964.26958036, 3078.63612913, 4672.26090061, 2310.20822426, 4180.41248267, 2162.15029365, 2072.33397643, 1238.51032921, 832.357478845, 1488.41649353, 2593.59405836, 1824.3311299, 2309.75855988, 2853.29775447, 2151.58363423, 2434.44194297, 2082.03137409, 1198.78596379, 1997.28834526, 3626.80120077, 3532.99815015, 2924.8945852, 2101.1318558, 2466.86972041, 2506.56843199, 2352.67147924, 2754.38598176, 1932.30496879, 2459.48921061, 1639.6977746, 1695.42812081, 2380.69363485, 2640.55720245, 2121.15585755, 3588.5380412, 4855.28448834, 4857.60200313, 4765.37072491, 1887.71370931, 1942.00303845, 1926.04452804, 2570.73288389, 1574.74408319, 3067.96627803, 1827.77703955, 1778.39064842, 2280.49016348, 1949.46253097, 2136.73909356, 1495.06656829, 1870.72119252]]
    e_d_cent = [[1039.173740299393, 1435.9536198970777, 1778.3059807875993, 833.5600644307985, 944.8627199983013, 9.540124318824205, 823.8633204223586, 1401.4181419730958, 827.5312450161011, 1163.202855888643, 13.237829767759157, 1792.3207837264274, 1560.168479675753, 1148.1705423534424, 755.9245287150616, 892.191969353227, 884.4665374591882, 415.35433225392086, 1875.7404289875428, 442.41685887935273, 797.3852853246149, 893.6966802217544, 715.515865771003, 1468.0163280348056, 180.07593240624257, 310.93678012208176, 1539.3048577844493, 645.5935589311879, 1221.0569425039175, 622.8771608179757, 777.3761195766289, 830.7302823568057, 74.42330047098997, 821.8292306571432, 1103.510336489361, 940.0405957271148, 822.4108254121722, 608.4887378615729, 1148.9388310016825, 587.5570422769338, 1526.421120290217, 1168.2108316701813, 1273.2683945252095, 393.432086570192, 1557.4353087829875, 844.2681890359437, 770.8856764819485, 760.2515378456518, 910.6248838134877, 598.6815415818562, 38.23452198067387, 978.0345113217967, 34.163250596628764, 36.151558902120726, 1015.1877769341709, 494.3376980922027, 68.14875955525937, 1541.2904129789106, 17.571308861202773, 1792.3673345822394, 458.22021965535055, 1283.3482722452554, 825.5956710209496, 1353.1526557220127, 1291.585217521375, 1159.6867266474756, 599.8811014721533, 1580.683975865727, 857.1852247427792, 1537.8170745504758, 1243.6395933750302, 1169.2812348399343, 365.9048145832149, 1640.885818719533, 13.585851453369786, 5.368550393503098, 1899.2919526321182, 1839.912870204379, 1621.6520186557234, 1144.4604372844262, 875.8209539282942, 1765.885395992803, 762.8973515574974, 1847.1410426886998, 1644.595816646902, 1643.8858310928317, 489.95206768810164, 16.70608038529138, 10.544177978144234], [1785.2629840823188, 1041.566284758639, 1467.6674047933768, 1923.8136494672647, 1716.0321102312591, 915.244895489895, 221.1675461512587, 861.867312378096, 840.6042035290191, 844.2414907308923, 1314.2330413768143, 1770.5151998468243, 1633.4007149260217, 368.65492682278284, 1058.123240416804, 1840.4439280460017, 284.2266672489778, 1506.0332777066492, 514.4581837549737, 878.2567435346746, 1708.3601232217165, 472.68214070414155, 997.4875266424194, 1954.127788865912, 1228.248441899879, 856.8220621495436, 221.15248443469196, 2353.736947917064, 1397.013963090411, 1120.4495754771256, 299.74838064309364, 1170.9093515493155, 414.38923232920945, 912.2142994105723, 1534.6450360693136, 645.6729294664426, 264.5515241024804, 345.43316624042524, 247.75738755022581, 653.3364213144663, 1445.3983698014852, 1050.9281645823662, 344.67425076422523, 1995.0077769507548, 1475.4497099997725, 825.5266504486219, 2150.14206872408, 454.5676595110679, 368.51923389480214, 336.38429850370164, 1890.5055794165903, 351.828770708083, 994.9605893924067, 2229.0283582517486, 548.3166216178901, 1742.920966906669, 1984.4888736732746, 566.106410249732, 2434.639259730216, 622.6818281879467, 936.9663122493897, 318.60972914001735, 1931.1872592082586, 165.28422562460204, 89.07313767059404, 423.11626224110586, 243.6303481764325, 1816.3310639967779, 435.7639805821094, 270.7290811126059, 2139.7993528686457, 310.6312946105225, 629.1766510798647, 907.0254731098203, 580.6065658166992, 1309.0795526339405, 681.692149849808, 155.04438853570673, 92.31684124911177, 136.71413051358266, 630.9594773060937, 443.11461492136965, 209.3880413379435, 110.5204345895767, 2013.1230449025077, 1518.1460692179344, 1951.8103229039705, 762.6417741189457, 290.3259190547951, 1400.7507173640981, 442.85762169709034, 1086.835709620925, 914.7013000175225, 367.2400890475026, 2182.627764031268, 1037.4172952821752, 110.51264193361975, 198.1143051051561, 698.273859249986, 1365.9485757880254, 1122.979469109, 413.50630387776516, 1243.5706321936586, 692.719607065467, 229.5825564846332, 261.3031238931711, 2347.1190868041, 2194.0422166832936, 1670.1045987321418, 1841.9016208506491, 299.8714988772865, 1624.557624100974, 1063.2981681602307, 1820.151128560012, 1081.256744265969, 2280.823993202776, 2062.179945562474, 1302.6619137307714, 206.9994601689073, 1212.1724246853282, 2013.6803279587316, 1956.0130420885787, 1593.5985662002408, 1725.1768982775568, 2138.5128286813383, 1919.6897087683137, 315.76110544137276, 1470.8243533267491, 1766.3716811308943, 332.203147610922, 224.8820964437035, 2478.4850127560644, 2377.8401976977834, 207.7094605398353, 132.62699690261672, 1238.1214902654456, 110.50946129093892, 293.3991946533112, 2113.9358140607815, 1403.4104181381058, 1943.7259970491498, 351.0171521341531, 1318.071612339641, 331.1161250625808, 286.95452010783447, 1695.0764266614954, 1387.1609196274565, 1682.152980716172, 395.94052375974496, 1337.190585009628]]
    aarr = [[[9.0, 8.9, 8.75, 9.9, 8.2, 8.1, 9.1, 7.9, 9.6, 9.6, 8.9, 9.6, 9.1, 9.2, 9.2, 9.45, 9.1, 9.3, 7.4, 9.6, 9.85, 9.1, 9.7, 9.55, 9.35, 9.15, 9.4, 9.0, 8.45, 9.7, 8.95, 8.0, 7.95, 9.45, 9.45, 9.65, 9.0, 9.05, 8.8, 9.3, 9.05, 8.0, 8.9, 9.5, 7.8, 9.8, 9.15, 9.8, 9.6, 9.55, 9.6, 9.45, 9.7, 9.75, 8.95, 9.6, 8.2, 8.8, 9.2, 8.7, 8.3, 8.4, 9.1, 8.5, 8.8, 8.35, 8.75, 9.6, 8.8, 8.3, 9.4, 8.3, 8.4, 8.0, 8.0, 8.5, 7.8, 8.7, 8.0, 8.4, 8.0, 8.1, 8.0, 8.0, 8.2, 6.9, 7.0, 7.2, 6.7], [9.20411998265592, 9.14612803567824, 8.79934054945358, 9.82607480270083, 8.04139268515822, 8.09691001300806, 8.79934054945358, 7.90308998699194, 9.77815125038364, 9.73239375982297, 9.17609125905568, 9.36172783601759, 9.30102999566398, 9.39794000867204, 9.20411998265592, 9.32221929473392, 9.10037054511756, 9.44715803134222, 7.39794000867204, 9.81291335664286, 9.77815125038364, 9.04139268515822, 9.47712125471966, 9.73239375982297, 9.07918124604762, 9.06069784035361, 9.32221929473392, 9.04139268515822, 8.39794000867204, 9.96848294855393, 8.77815125038364, 8.20411998265593, 8.04139268515822, 9.57978359661681, 9.61278385671973, 9.49136169383427, 9.14612803567824, 9.30102999566398, 8.5051499783199, 9.25527250510331, 9.07918124604762, 8.14612803567824, 9.0, 9.68124123737559, 7.39794000867204, 9.73239375982297, 9.30102999566398, 9.67209785793572, 9.63346845557959, 9.77815125038364, 9.49136169383427, 9.32221929473392, 9.80617997398389, 9.51851393987789, 9.07918124604762, 9.77815125038364, 8.14612803567824, 8.69897000433602, 9.44715803134222, 8.35, 8.25, 8.4, 9.15, 8.4, 8.65, 8.1, 8.65, 9.45, 8.45, 8.1, 9.25, 8.4, 7.9, 8.1, 8.3, 8.05, 7.7, 8.35, 7.9, 8.1, 8.0, 7.9, 7.8, 8.1, 8.4, 8.34242268082221, 7.9, 8.15, 8.1]], [[8.45, 8.2, 8.2, 8.7, 8.9, 9.3, 8.9, 9.15, 8.4, 8.6, 8.8, 9.1, 8.95, 8.8, 8.3, 8.6, 9.0, 9.5, 8.8, 8.6, 8.7, 9.1, 9.05, 9.1, 9.2, 8.6, 9.0, 8.7, 8.3, 9.05, 9.1, 9.1, 9.0, 9.0, 8.2, 9.2, 8.9, 9.05, 9.05, 9.3, 9.1, 9.15, 8.65, 9.4, 9.2, 9.05, 8.8, 9.5, 7.9, 9.0, 8.4, 8.95, 9.2, 8.3, 9.4, 8.3, 8.9, 9.05, 9.5, 9.2, 9.05, 8.0, 8.2, 9.1, 9.15, 8.3, 8.9, 8.65, 9.1, 8.85, 9.1, 8.7, 8.8, 9.0, 9.3, 9.25, 8.75, 9.05, 9.3, 9.4, 9.3, 9.45, 8.95, 9.3, 8.0, 9.3, 9.1, 8.65, 9.1, 8.6, 8.95, 9.4, 9.15, 9.5, 8.75, 8.8, 9.55, 9.15, 9.2, 9.1, 9.15, 9.2, 9.45, 9.2, 9.15, 9.0, 8.8, 8.0, 8.0, 9.5, 7.9, 9.0, 8.3, 8.8, 8.4, 8.8, 8.85, 8.1, 9.55, 8.3, 9.15, 8.8, 8.2, 7.9, 8.8, 8.3, 8.9, 8.1, 8.3, 8.0, 8.4, 8.7, 8.6, 9.25, 9.2, 9.45, 9.25, 8.8, 9.1, 8.8, 9.55, 7.9, 8.6, 9.1, 7.0, 6.7, 7.7, 7.3, 8.85, 7.5], [8.0, 8.1, 8.25, 8.69897000433602, 9.14612803567824, 9.30102999566398, 9.14612803567824, 9.11394335230684, 8.20411998265593, 8.5051499783199, 8.95424250943932, 9.1, 8.75, 8.60205999132796, 8.35, 8.45, 9.04139268515822, 9.39794000867204, 8.60205999132796, 8.69897000433602, 8.7, 9.09691001300805, 9.11394335230684, 8.89762709129044, 9.30102999566398, 8.55, 8.9, 8.75, 8.25, 9.23044892137827, 9.0, 8.89762709129044, 9.04139268515822, 9.14612803567824, 8.0, 9.25527250510331, 8.69897000433602, 9.0, 8.84509804001426, 9.20411998265592, 9.0, 9.11394335230684, 8.55, 9.23044892137827, 9.14612803567824, 9.17609125905568, 8.69897000433602, 9.34242268082221, 7.9, 9.09691001300805, 8.4, 9.20411998265592, 9.39794000867204, 8.3, 9.41497334797082, 8.05, 8.95424250943932, 8.79934054945358, 9.20411998265592, 9.20411998265592, 9.04139268515822, 8.25, 8.25, 9.14612803567824, 9.17609125905568, 8.09691001300806, 8.69897000433602, 8.5051499783199, 9.17609125905568, 8.84509804001426, 9.17609125905568, 8.60205999132796, 8.79934054945358, 8.7481880270062, 9.27875360095283, 9.20411998265592, 8.60205999132796, 8.95424250943932, 9.23044892137827, 9.14612803567824, 9.38021124171161, 9.36172783601759, 8.85125834871907, 9.25527250510331, 8.2, 9.20411998265592, 9.11394335230684, 8.6, 9.14612803567824, 8.60205999132796, 9.07918124604762, 9.25527250510331, 9.17609125905568, 9.34242268082221, 8.65321251377534, 8.69897000433602, 9.39794000867204, 9.04139268515822, 9.25527250510331, 9.20411998265592, 9.20411998265592, 9.23044892137827, 9.36172783601759, 9.23044892137827, 9.17609125905568, 9.14612803567824, 8.3, 8.3, 7.69897000433602, 9.25527250510331, 8.14612803567824, 9.0, 8.14612803567824, 8.60205999132796, 8.45, 8.55, 8.9, 8.2, 9.30102999566398, 8.39794000867204, 9.07918124604762, 8.60205999132796, 8.1, 8.25, 8.39794000867204, 8.45, 8.65321251377534, 7.95, 8.11394335230684, 8.39794000867204, 8.0, 8.20411998265593, 8.1, 9.11394335230684, 9.04139268515822, 9.47712125471966, 9.20411998265592, 8.60205999132796, 9.20411998265592, 8.79934054945358, 9.25527250510331, 8.15, 8.34242268082221, 9.0, 8.15, 8.3, 8.25, 7.9, 7.69897000433602, 8.35]]]
    darr = [[[18.92, 18.88, 19.04, 18.98, 18.88, 18.96, 18.94, 18.9, 19.06, 19.0, 18.96, 19.06, 19.04, 19.04, 18.88, 18.9, 19.02, 18.98, 18.9, 19.0, 18.98, 18.88, 19.0, 18.88, 18.94, 18.98, 18.86, 19.02, 18.92, 18.94, 19.04, 18.98, 18.96, 19.04, 18.86, 18.98, 18.88, 18.98, 19.0, 19.0, 19.04, 18.98, 19.06, 18.94, 18.9, 19.0, 19.02, 19.02, 19.0, 19.02, 18.96, 18.98, 18.96, 18.96, 18.92, 18.94, 18.96, 19.06, 18.96, 19.04, 18.94, 19.06, 18.92, 18.9, 18.98, 18.88, 18.94, 19.04, 18.92, 18.88, 18.88, 18.9, 18.94, 18.88, 18.96, 18.96, 19.02, 18.88, 18.9, 18.9, 18.94, 19.04, 18.94, 18.9, 18.88, 18.88, 18.94, 18.96, 18.96], [18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9, 18.9]], [[18.42, 18.54, 18.44, 18.44, 18.56, 18.54, 18.48, 18.46, 18.48, 18.54, 18.56, 18.42, 18.42, 18.48, 18.44, 18.58, 18.48, 18.44, 18.52, 18.52, 18.58, 18.52, 18.42, 18.42, 18.54, 18.44, 18.48, 18.6, 18.44, 18.54, 18.48, 18.56, 18.52, 18.42, 18.44, 18.44, 18.48, 18.52, 18.48, 18.44, 18.56, 18.48, 18.5, 18.6, 18.56, 18.46, 18.42, 18.46, 18.48, 18.48, 18.58, 18.48, 18.6, 18.52, 18.44, 18.42, 18.6, 18.44, 18.58, 18.54, 18.52, 18.5, 18.58, 18.5, 18.5, 18.48, 18.5, 18.42, 18.52, 18.48, 18.56, 18.5, 18.52, 18.46, 18.52, 18.58, 18.52, 18.5, 18.5, 18.5, 18.42, 18.52, 18.5, 18.5, 18.52, 18.6, 18.6, 18.44, 18.5, 18.54, 18.5, 18.42, 18.52, 18.48, 18.6, 18.52, 18.5, 18.48, 18.44, 18.56, 18.56, 18.46, 18.54, 18.44, 18.5, 18.5, 18.54, 18.52, 18.44, 18.58, 18.5, 18.42, 18.54, 18.42, 18.44, 18.4, 18.44, 18.44, 18.48, 18.56, 18.6, 18.42, 18.42, 18.4, 18.58, 18.58, 18.48, 18.42, 18.54, 18.48, 18.5, 18.6, 18.58, 18.48, 18.5, 18.6, 18.5, 18.48, 18.42, 18.44, 18.4, 18.5, 18.56, 18.48, 18.5, 18.56, 18.44, 18.42, 18.48, 18.54], [18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, -9999999999.9, 18.5, -9999999999.9, 18.5, -9999999999.9, -9999999999.9, -9999999999.9, 18.5, 18.5, 18.5, -9999999999.9, 18.5, -9999999999.9, -9999999999.9, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, -9999999999.9, 18.5, -9999999999.9, 18.5, -9999999999.9, 18.5, -9999999999.9, 18.5, 18.5, 18.5, -9999999999.9, 18.5, -9999999999.9, -9999999999.9, -9999999999.9, 18.5, -9999999999.9, -9999999999.9, -9999999999.9, -9999999999.9, 18.5, 18.5, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, -9999999999.9, 18.5, -9999999999.9, 18.5, -9999999999.9, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5]]]
    dsigma = [[[0.05, 0.05, 0.06, 0.07, 0.03, 0.05, 0.06, 0.05, 0.05, 0.07, 0.06, 0.06, 0.05, 0.04, 0.06, 0.04, 0.05, 0.06, 0.07, 0.05, 0.05, 0.07, 0.07, 0.06, 0.06, 0.05, 0.07, 0.05, 0.05, 0.06, 0.04, 0.05, 0.06, 0.04, 0.06, 0.05, 0.05, 0.06, 0.04, 0.07, 0.05, 0.06, 0.04, 0.03, 0.06, 0.05, 0.04, 0.05, 0.05, 0.05, 0.06, 0.06, 0.04, 0.06, 0.05, 0.05, 0.06, 0.05, 0.03, 0.06, 0.04, 0.04, 0.06, 0.05, 0.05, 0.04, 0.05, 0.06, 0.04, 0.06, 0.05, 0.04, 0.04, 0.06, 0.04, 0.03, 0.07, 0.07, 0.06, 0.04, 0.06, 0.06, 0.05, 0.07, 0.06, 0.06, 0.03, 0.05, 0.06], [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]], [[0.05, 0.05, 0.06, 0.05, 0.06, 0.07, 0.06, 0.04, 0.07, 0.05, 0.06, 0.04, 0.06, 0.05, 0.06, 0.06, 0.05, 0.06, 0.05, 0.06, 0.04, 0.05, 0.05, 0.06, 0.06, 0.05, 0.04, 0.05, 0.06, 0.06, 0.05, 0.03, 0.05, 0.06, 0.06, 0.03, 0.05, 0.06, 0.06, 0.06, 0.05, 0.04, 0.07, 0.04, 0.06, 0.03, 0.06, 0.06, 0.04, 0.06, 0.04, 0.05, 0.04, 0.06, 0.04, 0.06, 0.06, 0.03, 0.07, 0.07, 0.07, 0.05, 0.06, 0.03, 0.04, 0.06, 0.06, 0.06, 0.05, 0.06, 0.06, 0.06, 0.06, 0.04, 0.04, 0.06, 0.06, 0.06, 0.04, 0.06, 0.05, 0.05, 0.06, 0.05, 0.06, 0.06, 0.06, 0.05, 0.05, 0.07, 0.05, 0.07, 0.06, 0.05, 0.04, 0.07, 0.05, 0.07, 0.05, 0.06, 0.03, 0.05, 0.06, 0.05, 0.05, 0.06, 0.06, 0.07, 0.03, 0.04, 0.06, 0.05, 0.07, 0.06, 0.04, 0.05, 0.03, 0.04, 0.04, 0.07, 0.04, 0.07, 0.06, 0.02, 0.05, 0.06, 0.06, 0.04, 0.06, 0.04, 0.06, 0.06, 0.05, 0.04, 0.06, 0.05, 0.05, 0.02, 0.07, 0.05, 0.06, 0.06, 0.06, 0.06, 0.05, 0.05, 0.05, 0.04, 0.06, 0.03], [0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, -9999999999.9, 0.1, -9999999999.9, 0.1, -9999999999.9, -9999999999.9, -9999999999.9, 0.1, 0.1, 0.1, -9999999999.9, 0.1, -9999999999.9, -9999999999.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, -9999999999.9, 0.1, -9999999999.9, 0.1, -9999999999.9, 0.1, -9999999999.9, 0.1, 0.1, 0.1, -9999999999.9, 0.1, -9999999999.9, -9999999999.9, -9999999999.9, 0.1, -9999999999.9, -9999999999.9, -9999999999.9, -9999999999.9, 0.1, 0.1, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, -9999999999.9, 0.1, -9999999999.9, 0.1, -9999999999.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]]
    gal_names = [['B112', 'B47', 'H86-70', 'HW59', 'L63', 'L62', 'BS265', 'L50', 'AM3', 'HW40', 'B39', 'BS121', 'BS88', 'HW55', 'L106', 'L19', 'H86-197', 'HW67', 'L49', 'L112', 'HW22', 'L111', 'L109', 'HW63', 'HW84', 'L3', 'L100', 'HW86', 'L34', 'HW42', 'L102', 'L30', 'L115', 'HW66', 'L5', 'K38', 'L108', 'L58', 'NGC294', 'HW85', 'B34', 'L45', 'L28', 'HW31', 'L72', 'NGC339', 'L7', 'HW79', 'L91', 'L113', 'L4', 'L27', 'L110', 'L6', 'NGC419', 'HW41', 'L114', 'BS35', 'HW47', 'SOGLE196', 'K63', 'B103', 'B111', 'H86-76', 'H86-174', 'K43', 'K57', 'BS80', 'K55', 'HW52', 'BS75', 'B55', 'HW32', 'B99', 'K61', 'L39', 'H86-190', 'NGC241', 'B48', 'H86-87', 'B124', 'H86-85', 'NGC242', 'H86-97', 'H86-90', 'L35', 'K47', 'B134', 'H86-188'], ['BRHT4B', 'BRHT45A', 'BRHT38B', 'BSDL1035', 'BSDL527', 'BSDL3158', 'KMHK128', 'HS114', 'BSDL1024', 'C11', 'BSDL665', 'HS131', 'KMHK505', 'BSDL716', 'H88-33', 'HS411', 'BSDL675', 'H3', 'HS38', 'HS154', 'H88-331', 'KMHK112', 'KMHK151', 'HS151', 'H88-334', 'KMHK95', 'BSDL77', 'H88-235', 'BSDL677', 'H88-67', 'BSDL2995', 'H88-26', 'KMHK123', 'KMHK1719', 'NGC1838', 'LW397', 'H88-55', 'KMHK1702', 'NGC1697', 'KMHK58', 'KMHK229', 'NGC1917', 'HS116', 'KMHK1023', 'H88-52', 'HS121', 'NGC1865', 'SL133', 'SL230', 'NGC2108', 'KMHK907', 'NGC1795', 'SL13', 'H88-265', 'NGC1997', 'NGC1793', 'SL505', 'LW469', 'SL359', 'HS8', 'HS156', 'HS390', 'SL269', 'BSDL594', 'NGC1644', 'NGC1839', 'SL154', 'SL229', 'SL678', 'LW224', 'SL35', 'SL293', 'HS178', 'KMHK506', 'SL555', 'SL73', 'SL548', 'SL54', 'KMHK1668', 'SL874', 'OHSC28', 'SL674', 'H88-40', 'LW263', 'SL397', 'LW69', 'NGC1751', 'LW54', 'NGC1852', 'SL300', 'SL290', 'LW211', 'SL151', 'SL446A', 'SL351', 'SL444', 'SL5', 'SL870', 'LW393', 'SL132', 'SL96', 'SL549', 'SL707', 'SL869', 'SL162', 'NGC1846', 'OGLE298', 'H88-244', 'NGC1863', 'HS329', 'HS130', 'KMHK1055', 'BSDL761', 'SL588', 'BSDL341', 'HS247', 'H88-269', 'H88-320', 'SL33', 'SL72', 'SL244', 'NGC1836', 'HS412', 'H88-307', 'NGC1860', 'KMHK378', 'H88-188', 'BSDL268', 'SL510', 'NGC2093', 'BSDL779', 'H88-245', 'H88-279', 'SL41', 'NGC2161', 'SL663', 'IC2146', 'H88-333', 'HS264', 'KMHK1045', 'KMHK586', 'SL551', 'BSDL654', 'H88-131', 'SL579', 'KMHK975', 'H88-316', 'KMHK979', 'SL218', 'BSDL631']]

    in_params = {'ra': ra, 'dec': dec, 'dist_cent': dist_cent,
                 'e_d_cent': e_d_cent, 'aarr': aarr, 'darr': darr,
                 'dsigma': dsigma, 'gal_names': gal_names}

    gal_str_pars, rho_plot_pars = gsd(in_params)

    from make_all_plots import make_angles_plot, make_rho_min_plot
    make_angles_plot(gal_str_pars)
    make_rho_min_plot(rho_plot_pars)
