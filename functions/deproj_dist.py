
from astropy.coordinates import Distance, Angle, SkyCoord
from astropy import units as u
import numpy as np


def rho_phi(coord, glx_ctr):
    '''
    Eqs 1, 2 & 3 from van der Marel & Cioni (2001).
    '''
    # # Angular distance between point and center of galaxy.
    # cos_rho = np.cos(coord.dec.radian) * np.cos(glx_ctr.dec.radian) * \
    #     np.cos(coord.ra.radian - glx_ctr.ra.radian) + \
    #     np.sin(coord.dec.radian) * np.sin(glx_ctr.dec.radian)
    # rho = Angle(np.arccos(cos_rho) * 180. / np.pi, unit=u.deg)
    # print 'rho = ', rho.degree

    # # Position angle.
    # cos_phi = (-1. * np.cos(coord.dec.radian) * np.sin(coord.ra.radian -
    #            glx_ctr.ra.radian)) / np.sin(rho.radian)
    # phi = np.arccos(cos_phi) * 180. / np.pi

    # Position angle between center and coordinates. This is the angle between
    # the positive y axis (North) counter-clockwise towards the negative x
    # axis (East).
    Phi = glx_ctr.position_angle(coord)
    delta_ra = coord.ra.degree - glx_ctr.ra.degree
    delta_dec = coord.dec.degree - glx_ctr.dec.degree
    if delta_ra <= 0. and delta_dec >= 0.:
        print 'NW', Phi.degree
    elif delta_ra >= 0. and delta_dec >= 0.:
        print 'NE', Phi.degree
    elif delta_ra >= 0. and delta_dec <= 0.:
        print 'SE', Phi.degree
    elif delta_ra <= 0. and delta_dec <= 0.:
        print 'SW', Phi.degree
    phi = Phi + Angle('90d')

    # Angular separation between center and coordinates.
    rho = coord.separation(glx_ctr)

    return rho, Phi, phi


def gal_theta(glx_PA):
    '''
    PA of the galaxy, rotated 90 deg.
    '''
    theta = glx_PA - Angle('90d')

    return theta


def phi_rotate(coord, glx_ctr, phi):
    '''
    Rotate so that phi is located in the NW quadrant.
    '''
    phi_rot = Angle(90. * int(phi.degree / 90.0001), unit='degree')
    phi = phi - phi_rot
    print phi.degree

    return phi


def vdm_2001_D(glx_incl, glx_dist, rho, phi, theta):
    '''
    Eq 8 from van der Marel & Cioni (2001).
    '''
    D_0 = glx_dist.value
    # Distance to the coordinates passed.
    A = np.cos(glx_incl.radian) * np.cos(rho.radian) - \
        np.sin(glx_incl.radian) * np.sin(rho.radian) * \
        np.sin(phi.radian - theta.radian)
    D = D_0 * np.cos(glx_incl.radian) / (A)

    return D, D_0


def vdm_2001_xy(phi, rho, D):
    '''
    Eq 5 from van der Marel & Cioni (2001).
    '''
    x = D * np.sin(rho.radian) * np.cos(phi.radian)
    y = D * np.sin(rho.radian) * np.sin(phi.radian)
    x, y = Angle(x, unit='degree'), Angle(y, unit='degree')
    # print x.degree, y.degree
    # z = D_0 - D * np.cos(rho.radian)  # <-- not used??

    return x, y


def carrera_2011_xy(rho, phi):
    '''
    Carrera, private communication.
    '''
    y = rho.radian * np.cos(phi.radian)
    x = rho.radian * np.sin(phi.radian)
    x, y = Angle(x, unit=u.rad), Angle(y, unit=u.rad)
    # print x.degree, y.degree

    return x, y


def carrera_2011_xy_inv(rho, phi):
    '''
    '''
    x = rho.radian * np.cos(phi.radian)
    y = rho.radian * np.sin(phi.radian)
    x, y = Angle(x, unit=u.rad), Angle(y, unit=u.rad)
    # print x.degree, y.degree

    return x, y


def vdm_2001_dep_dist(rho, phi, theta, glx_incl, D, D_0):
    '''
    Deprojected angular distance from vdM & Cioni (2001).
    '''
    # Eq 7 from van der Marel & Cioni (2001).
    x_p = D * np.sin(rho.radian) * np.cos(phi.radian - theta.radian)
    y_p = D * (np.sin(rho.radian) * np.cos(glx_incl.radian) *
               np.sin(phi.radian - theta.radian) + np.cos(rho.radian) *
               np.sin(glx_incl.radian)) - D_0 * np.sin(glx_incl.radian)
    x_p, y_p = Angle(x_p, unit='degree'), Angle(y_p, unit='degree')
    print 'M01 x_p,y_p=', x_p.degree, y_p.degree
    # z_p = D * (np.sin(rho.radian) * np.sin(glx_incl.radian) *
    #            np.sin(phi.radian - theta.radian) - np.cos(rho.radian) *
    #            np.cos(glx_incl.radian)) + D_0 * np.cos(glx_incl.radian)
    dep_dist_deg = Angle(np.sqrt(x_p ** 2 + y_p ** 2), unit='degree')

    return dep_dist_deg


def claria_2005_dep_dist(rho, phi, glx_incl, theta):
    '''
    Deprojected distance from Claria et al. 2005.

    This formula is obtained from the Cioni (2009) eqs assuming:
    x = rho*cos(phi)
    y = rho*sin(phi)
    '''
    A = 1 + (np.sin(phi.radian - theta.radian) * np.tan(glx_incl)) ** 2
    dep_dist_deg = Angle(rho * np.sqrt(A), unit='degree')

    return dep_dist_deg


def cioni_2009_dep_dist(glx_incl, theta, x, y):
    '''
    Deprojected angular distance. Eqs 1, 2, 3 & 4 from Cioni (2009).
    '''
    # Rotate the coords system.
    x1 = x * np.cos(theta.radian) + y * np.sin(theta.radian)
    y1 = y * np.cos(theta.radian) - x * np.sin(theta.radian)
    # De-project.
    y2 = y1 / np.cos(glx_incl.radian)
    print 'C09 x1,y2', x1.degree, y2.degree
    # Obtain de-projected distance in decimal degrees.
    dep_dist_deg = Angle(np.sqrt(x1 ** 2 + y2 ** 2), unit='degree')

    return dep_dist_deg


def cioni_2009_dist_kpc(dep_dist_deg, glx_dist):
    '''
    Obtain de-projected distance in the units used for the galaxy center
    distance.
    Eq 5 from Cioni (2009).
    '''
    dep_dist_kpc = Distance(np.tan(dep_dist_deg) * glx_dist,
                            unit=glx_dist.unit)

    return dep_dist_kpc


def deproj_dist(coord,
                glx_ctr=SkyCoord('00h42m44.33s +41d16m07.5s', frame='icrs'),
                glx_PA=Angle('37d42m54s'), glx_incl=Angle('77.5d'),
                glx_dist=Distance(783, unit=u.kpc)):
    """
    Computes deprojected galactocentric distance.

    Based on: https://gist.github.com/jonathansick/9399842

    Parameters
    ----------
    coord : :class:`astropy.coordinates.ICRS`
        Coordinate of points to compute galactocentric distance for.
        Can be either a single coordinate, or array of coordinates.
    glx_ctr : :class:`astropy.coordinates.ICRS`
        Galaxy center.
    glx_PA : :class:`astropy.coordinates.Angle`
        Position angle of galaxy disk.
    glx_incl : :class:`astropy.coordinates.Angle`
        Inclination angle of the galaxy disk.
    glx_dist : :class:`astropy.coordinates.Distance`
        Distance to galaxy.

    Returns
    -------
    dep_dist_deg : class:`astropy.coordinates.Angle`
        Galactocentric distance(s) for coordinate point(s) in decimal degrees.
    dep_dist_kpc : class:`astropy.coordinates.Distance`
        Galactocentric distance(s) for coordinate point(s).
    """

    theta = gal_theta(glx_PA)
    rho, Phi, phi = rho_phi(coord, glx_ctr)

    D, D_0 = vdm_2001_D(glx_incl, glx_dist, rho, phi, theta)
    print 'rho', rho.radian, rho.degree, np.tan(rho)
    dep_dist_deg = vdm_2001_dep_dist(rho, phi, theta, glx_incl, D, D_0)
    print 'd_deg M01 =', dep_dist_deg.degree

    # # phi_r = phi_rotate(coord, glx_ctr, phi)
    dep_dist_deg = claria_2005_dep_dist(rho, phi, glx_incl, theta)
    print 'd_deg C05 =', dep_dist_deg.degree

    x, y = vdm_2001_xy(phi, rho, D)
    print 'vdm, xy=', x.degree, y.degree
    dep_dist_deg = cioni_2009_dep_dist(glx_incl, theta, x, y)
    print 'd_deg C09_1 =', dep_dist_deg.degree

    x, y = carrera_2011_xy(rho, phi)
    # print 'C11, xy=', x.degree, y.degree
    dep_dist_deg = cioni_2009_dep_dist(glx_incl, theta, x, y)
    print 'd_deg C09_2 =', dep_dist_deg.degree

    x, y = carrera_2011_xy_inv(rho, phi)
    # print 'C11, xy=', x.degree, y.degree
    dep_dist_deg = cioni_2009_dep_dist(glx_incl, theta, x, y)
    print 'd_deg C09_3 =', dep_dist_deg.degree

    dep_dist_kpc = cioni_2009_dist_kpc(dep_dist_deg, glx_dist)
    # print 'd_proj =', dep_dist_kpc, '\n'
    print '\n'

    return dep_dist_kpc


def main():
    cent = SkyCoord('05h27.6m', '-69.87d', frame='icrs')
    dist = Distance(51., unit=u.kpc)
    inc, pa = Angle('34.7d'), Angle('189.3d')  # LMC
    ra_lst = ['05h20m', '05h40m', '05h40m', '05h20m']
    dec_lst = ['-68d', '-68d', '-70d', '-70d']

    for i, (ra, dec) in enumerate(zip(*[ra_lst, dec_lst])):
        a = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        deproj_dist(a, cent, pa, inc, dist)

if __name__ == "__main__":
    main()
