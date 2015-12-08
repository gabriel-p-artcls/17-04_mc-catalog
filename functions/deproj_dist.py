
from astropy.coordinates import Distance, Angle, SkyCoord
from astropy import units as u
import numpy as np


def rho_pa(coord, glx_ctr):
    '''
    Equivalent to Eqs 1, 2 & 3 from van der Marel & Cioni (2001).

    Some rotations need to be applied to the 'phi' if using the equations
    from the article instead of astropy.
    '''
    # Angular separation between center and coordinates.
    rho = coord.separation(glx_ctr)
    # Position angle between center and coordinates.
    pa = glx_ctr.position_angle(coord)

    # # Generate angular distance and position angles as those shown in
    # # Table 1 of Carrera et al. (2011).

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

    # # Transform phi so it will match the position of the angles in Table 1.
    # delta_ra = coord.ra.degree - glx_ctr.ra.degree
    # delta_dec = coord.dec.degree - glx_ctr.dec.degree
    # if delta_ra <= 0. and delta_dec >= 0.:
    #     # This is the N-W quadrant.
    #     phi = phi + 360.
    # elif delta_ra >= 0. and delta_dec >= 0.:
    #     # N-E quadrant.
    #     phi = phi
    # elif delta_ra >= 0. and delta_dec <= 0.:
    #     # S-E quadrant.
    #     phi = 360. - phi
    # elif delta_ra <= 0. and delta_dec <= 0.:
    #     # S-W quadrant.
    #     phi = 360. - phi
    # # Obtain the position angle.
    # pa = Angle(phi - 90., unit=u.deg)

    # print 'rho, pa = {:.2f}, {:.1f}'.format(rho.degree * 60., pa.degree)

    return rho, pa


def x_y(rho, pa):
    '''
    Carrera, private communication.
    '''
    y = rho.radian * np.cos(pa.radian)
    x = rho.radian * np.sin(pa.radian)
    x, y = Angle(x, unit=u.rad), Angle(y, unit=u.rad)
    print 'x =', x.radian, 'y =', y.radian

    return x, y


# def claria_dist(rho, phi, glx_PA, glx_incl, glx_dist):
#     '''
#     Claria et al. (2005). This Eq only appears to work if glx_PA=45(deg) so
#     that several terms in the (more general?) expression from Cioni (2009)
#     are simplified.
#     '''
#     theta = glx_PA - Angle('90d')
#     A = 1 + np.sin((phi.radian - theta.radian)) ** 2 * \
#         (np.tan(glx_incl)) ** 2
#     d = float(rho.radian * np.sqrt(A))
#     print 'R_proj =', np.rad2deg(d)

#     # Added by me, not in Claria et al. (2005)
#     dep_dist_kpc = Distance(np.tan(d) * glx_dist, unit=glx_dist.unit)

#     return dep_dist_kpc


def cioni_dist(x, y, glx_PA, glx_incl, glx_dist):
    '''
    Eqs 1, 2, 3 & 4 from Cioni (2009).
    '''

    theta = glx_PA - Angle('90d')

    # Rotate the coords system.
    x1 = x * np.cos(theta.radian) + y * np.sin(theta.radian)
    y1 = y * np.cos(theta.radian) - x * np.sin(theta.radian)
    # print 'x1 =', x1.radian, 'y1 =', y1.radian

    # De-project.
    y2 = y1 / np.cos(glx_incl.radian)
    # print 'y2 =', y2.radian

    # Obtain de-projected distance in decimal degrees.
    dep_dist_rad = np.sqrt(x1.rad ** 2 + y2.rad ** 2)
    dep_dist_deg = Angle(np.rad2deg(dep_dist_rad), unit='degree')
    # print 'R_proj =', dep_dist_deg.degree

    # Obtain de-projected distance in the units used for the galaxy center
    # distance.
    dep_dist_kpc = Distance(np.tan(dep_dist_deg) * glx_dist,
                            unit=glx_dist.unit)

    return dep_dist_kpc


def van_der_marel_dist(rho, phi, coord, glx_ctr, glx_PA, glx_incl, glx_dist):
    '''
    Eqs derived from those given in van der Marel & Cioni (2001).
    '''

    theta = glx_PA - Angle('90d')

    A = np.cos(glx_incl.radian) * np.cos(rho.radian) - \
        np.sin(glx_incl.radian) * np.sin(rho.radian) * \
        np.sin(phi.radian - theta.radian)

    D_0 = glx_dist
    D = D_0 * np.cos(glx_incl.radian) / (A)
    # x = D * np.sin(rho.radian) * np.cos(phi.radian)
    # y = D * np.sin(rho.radian) * np.sin(phi.radian)
    # z = D_0 - D * np.cos(rho.radian)
    # print 'x =', x, 'y =', y

    x_p = D * np.sin(rho.radian) * np.cos(phi.radian - theta.radian)
    y_p = D * (np.sin(rho.radian) * np.cos(glx_incl.radian) *
               np.sin(phi.radian - theta.radian) + np.cos(rho.radian) *
               np.sin(glx_incl.radian)) - D_0 * np.sin(glx_incl.radian)

    # De-projected distance in Kpc.
    dep_dist_kpc = Distance(np.sqrt(x_p ** 2 + y_p ** 2), unit=glx_dist.unit)

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

    rho, pa = rho_pa(coord, glx_ctr)
    x, y = x_y(rho, pa)

    dep_dist_kpc = van_der_marel_dist(rho, pa, coord, glx_ctr, glx_PA,
                                      glx_incl, glx_dist)
    print 'vdM d_proj =', dep_dist_kpc
    dep_dist_kpc = cioni_dist(x, y, glx_PA, glx_incl, glx_dist)
    print 'Cio d_proj =', dep_dist_kpc, '\n'

    return dep_dist_kpc
