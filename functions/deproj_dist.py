
from astropy.coordinates import Distance, Angle, SkyCoord
from astropy import units as u
import numpy as np


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

    # Angular separation between center and coordinates.
    # a = coord.separation(glx_ctr)
    # print a.degree, a.degree * 60.

    # Generate angular distance and position angles as those shown in
    # Table 1 of Carrera et al. (2011).
    # Rho & Phi obtained following Eqs (1) and (2) in
    # van der Marel & Cioni (2001).

    # Angular distances.
    # cos_rho = np.cos(a.dec.radian) * np.cos(c_LMC.dec.radian) * \
    #     np.cos(a.ra.radian - c_LMC.ra.radian) + np.sin(a.dec.radian) * \
    #     np.sin(c_LMC.dec.radian)
    # rho = np.arccos(cos_rho) * 180. / np.pi
    # print 'rho = {:.4f}'.format(rho)

    # Position angles.
    # cos_phi = (-1. * np.cos(a.dec.radian) * np.sin(a.ra.radian -
    #            c_LMC.ra.radian)) / np.sin(np.arccos(cos_rho))
    # phi = np.arccos(cos_phi) * 180. / np.pi
    # # This is so phi matches the values in Table 1. Not sure why.
    # if phi < 45.:
    #     print 'phi = {:.2f}, {:.1f}'.format(phi, 180 + (90 - phi))
    # elif phi > 90.:
    #     print 'phi = {:.2f}, {:.1f}'.format(phi, 360 - (90 + phi))
    # else:
    #     print 'phi = {:.2f}, {:.1f}'.format(phi, 360 - (90 - phi))

    # Convert equatorial coords into angular coords.
    avg_dec = 0.5 * (glx_ctr.dec + coord.dec).radian
    if coord.ra > Angle(180., unit=u.deg):
        ra_inv = Angle((coord.ra.degree - 360.), unit=u.deg)
        x = (glx_ctr.ra - ra_inv) * np.cos(avg_dec)
    else:
        x = (glx_ctr.ra - coord.ra) * np.cos(avg_dec)
    y = glx_ctr.dec - coord.dec

    # Rotate the coords system. Eqs 1, 2, 3 & 4 from Cioni et al. (2009)
    phi = glx_PA - Angle('90d')
    x1 = x * np.cos(phi.radian) + y * np.sin(phi.radian)
    y1 = y * np.cos(phi.radian) - x * np.sin(phi.radian)
    # De-project.
    y2 = y1 / np.cos(glx_incl.radian)

    # Obtain de-projected distance in decimal degrees.
    dep_dist_rad = np.sqrt(x1.rad ** 2 + y2.rad ** 2)
    dep_dist_deg = Angle(np.rad2deg(dep_dist_rad), unit='degree')
    # Obtain de-projected distance in the units used for the galaxy center
    # distance.
    dep_dist_kpc = Distance(np.tan(dep_dist_deg) * glx_dist,
                            unit=glx_dist.unit)

    return dep_dist_deg, dep_dist_kpc
