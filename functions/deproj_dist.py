
from astropy.coordinates import Distance, Angle, SkyCoord
from astropy import units as u
import numpy as np


def rho_phi(coord, glx_ctr):
    '''
    Eqs 1, 2 & 3 from van der Marel & Cioni (2001).

    Parameters
    ----------
    coord : :class:`astropy.coordinates.ICRS`
        Coordinate of points to compute galactocentric distance for.
        Can be either a single coordinate, or array of coordinates.
    glx_ctr : :class:`astropy.coordinates.ICRS`
        Galaxy center.
    '''
    # Angular distance between point and center of galaxy.
    rho = coord.separation(glx_ctr)

    # Position angle between center and coordinates. This is the angle between
    # the positive y axis (North) counter-clockwise towards the negative x
    # axis (East).
    Phi = glx_ctr.position_angle(coord)

    # This is the angle measured counter-clockwise from the x positive axis
    # (West).
    phi = Phi + Angle('90d')

    # x = rho.degree * np.cos(phi)
    # y = rho.degree * np.sin(phi)
    # print x, y
    # if x >= 0. and y >= 0.:
    #     print 'NW', Phi.degree, phi.degree
    # elif x <= 0. and y >= 0.:
    #     print 'NE', Phi.degree, phi.degree
    # elif x <= 0. and y <= 0.:
    #     print 'SE', Phi.degree, phi.degree
    # elif x >= 0. and y <= 0.:
    #     print 'SW', Phi.degree, phi.degree

    return rho, Phi, phi


def gal_theta(glx_PA):
    '''
    PA of the galaxy, rotated 90 deg as defined in vdM et al. (2002).
    '''
    theta = glx_PA + Angle('90d')

    return theta


def vdm_2001_xy(rho, phi):
    '''
    Eq 4 from van der Marel & Cioni (2001).
    '''
    # x,y in degrees.
    x = rho.degree * np.cos(phi.radian)
    y = rho.degree * np.sin(phi.radian)

    return x, y


# def carrera_2011_xy(rho, Phi):
#     '''
#     Carrera, private communication said:
#     x = rho * np.sin(Phi)
#     y = rho * np.cos(Phi)
#     but did not specify if they used Phi or phi.
#     '''
#     # phi = Phi.radian + np.pi / 2.
#     # x = -1. * rho.degree * np.cos(phi)
#     # y = rho.degree * np.sin(phi)
#     x = rho.degree * np.sin(Phi.radian)
#     y = rho.degree * np.cos(Phi.radian)

#     return x, y

def vdm_2001_D(glx_incl, D_0, rho, phi, theta):
    '''
    Eq 8 from van der Marel & Cioni (2001).
    '''
    # Distance to the coordinates passed.
    # A = np.cos(glx_incl.radian) * np.cos(rho.radian) - \
    #     np.sin(glx_incl.radian) * np.sin(rho.radian) * \
    #     np.sin(phi.radian - theta.radian)
    # A = np.cos(a) * np.cos(b) - np.sin(a) * np.sin(b) * np.sin(c - d)
    s = np.sin(phi.radian - theta.radian)
    A = 0.5 * ((1-s) * np.cos(glx_incl.radian - rho.radian) +
               (1+s) * np.cos(glx_incl.radian + rho.radian))
    D = D_0 * np.cos(glx_incl.radian) / A

    return D


def vdm_2001_dep_dist(rho, phi, theta, glx_incl, D, D_0):
    '''
    Deprojected angular distance from vdM & Cioni (2001).

    Eq 7 from van der Marel & Cioni (2001) / Eq 2 from van der Marel (2001).
    '''
    # x_p = D * np.sin(rho.radian) * np.cos(phi.radian - theta.radian)

    # Both expressions are equivalent.
    # y_p = D * (np.sin(rho.radian) * np.cos(glx_incl.radian) *
    #            np.sin(phi.radian - theta.radian) + np.cos(rho.radian) *
    #            np.sin(glx_incl.radian)) - D_0 * np.sin(glx_incl.radian)
    # y_p = D * np.sin(rho.radian) * np.sin(phi.radian - theta.radian) / \
    #     np.cos(glx_incl.radian)

    # z_p = 0 since the source is located *on* the inclined disk.
    # z_p = D * (np.sin(rho.radian) * np.sin(glx_incl.radian) *
    #            np.sin(phi.radian - theta.radian) - np.cos(rho.radian) *
    #            np.cos(glx_incl.radian)) + D_0 * np.cos(glx_incl.radian)
    # d_kpc = np.sqrt(x_p ** 2 + y_p ** 2)

    # The above is equivalent to using the cosine law.
    d_kpc = np.sqrt(D_0 ** 2 + D ** 2 - 2 * D_0 * D * np.cos(rho.radian))

    return d_kpc


def vdm_2001_dep_dist_kpc(rho, phi, theta, glx_incl, D_0):
    '''
    Deprojected angular distance from vdM & Cioni (2001).

    Eq 7 from van der Marel & Cioni (2001) / Eq 2 from van der Marel (2001).
    '''
    s = np.sin(phi.radian - theta.radian)
    A = 0.5 * ((1-s)*np.cos(glx_incl.radian - rho.radian) +
               (1+s)*np.cos(glx_incl.radian + rho.radian))
    D = np.cos(glx_incl.radian) / A

    # The above is equivalent to using the cosine law.
    d_kpc = D_0*np.sqrt(1. + D**2 - 2*D*np.cos(rho.radian))

    return d_kpc


def cioni_2009_dep_dist(glx_incl, theta, x, y):
    '''
    Deprojected angular distance. Eqs 1, 2, 3 & 4 from Cioni (2009).
    '''
    # Rotate the coords system.
    x1 = x*np.cos(theta.radian) + y*np.sin(theta.radian)
    y1 = y*np.cos(theta.radian) - x*np.sin(theta.radian)
    # De-project.
    y2 = y1 / np.cos(glx_incl.radian)
    # Obtain de-projected distance in decimal degrees.
    dep_dist_deg = Angle(np.sqrt(x1**2 + y2**2), unit='degree')

    return dep_dist_deg


def cioni_2009_dist_kpc(dep_dist_deg, D_0):
    '''
    Obtain de-projected distance in the units used for the galaxy center
    distance.
    Eq 5 from Cioni (2009).
    '''
    dep_dist_kpc = Distance(np.tan(dep_dist_deg) * D_0, unit=D_0.unit)

    return dep_dist_kpc


def deproj_dist(glx_PA, glx_incl, glx_dist, rho, phi):
    """
    Computes deprojected galactocentric distance.

    Based on: https://gist.github.com/jonathansick/9399842

    Parameters
    ----------
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

    # Obtain 'theta' position angle for the galaxy.
    theta = gal_theta(glx_PA)
    # Convert distance to galaxy to kpc.
    D_0 = Distance(glx_dist.kpc, unit=u.kpc)

    # Obtain angular projected distance and position angle for the cluster.
    # rho, Phi, phi = rho_phi(coord, glx_ctr)

    # D = vdm_2001_D(glx_incl, D_0, rho, phi, theta)
    # dep_dist_kpc_M01 = vdm_2001_dep_dist(rho, phi, theta, glx_incl, D, D_0)
    dep_dist_kpc_M01 = vdm_2001_dep_dist_kpc(rho, phi, theta, glx_incl, D_0)
    # print 'd_kpc M01 = {:.12f}'.format(dep_dist_kpc_M01)

    # This gives the values for the deprojected angular distance in
    # Carrera et al. (2011).
    # x, y = carrera_2011_xy(rho, Phi)

    # vdm&C01 (x,y values) + Cioni 2009 (C09).
    # x, y = vdm_2001_xy(rho, phi)
    # dep_dist_deg = cioni_2009_dep_dist(glx_incl, theta, x, y)
    # print 'd_deg C09 = {:.8f}'.format(dep_dist_deg)

    # Convert distance value to kpc using C09 equation.
    # dep_dist_kpc = cioni_2009_dist_kpc(dep_dist_deg, D_0)
    # print 'd_kpc C09 = {:.12f}'.format(dep_dist_kpc)

    # # TEST.
    # s = np.sin(phi.radian - theta.radian)
    # A = 0.5 * ((1-s) * np.cos(glx_incl.radian - rho.radian) +
    #            (1+s) * np.cos(glx_incl.radian + rho.radian))
    # D = D_0 * np.cos(glx_incl.radian) / A
    # B = 1 + (np.sin(phi.radian - theta.radian) * np.tan(glx_incl)) ** 2
    # # Replace rho by D*sin(rho) --> Equivalent to vdMC01.
    # print 'd_kpc XXX = {:.12f}'.format(D * np.sin(rho.radian) * np.sqrt(B))
    # # TEST

    # print dep_dist_kpc_M01, dep_dist_kpc_M01 - dep_dist_kpc

    # return dep_dist_deg, dep_dist_deg - dep_dist_deg_p
    return dep_dist_kpc_M01


if __name__ == "__main__":

    # # Values in all 4 quadrants.
    # cent = SkyCoord('05h27.6m', '-69.87d', frame='icrs')
    # dist = Distance(51., unit=u.kpc)
    # inc, pa = Angle('34.7d'), Angle('189.3d')  # LMC
    # ra_lst = ['05h10m', '05h50m', '05h50m', '05h10m']
    # dec_lst = ['-66d', '-66d', '-72d', '-72d']

    # # Noel SMC data.
    # cent = SkyCoord('0h52m42s', '-72d49m', frame='icrs')
    # dist = Distance(51., unit=u.kpc)
    # inc, pa = Angle('34.7d'), Angle('189.3d')  # LMC
    # ra_lst = ['0h57m', '0h37m', '0h36m', '1h11m', '1h12m', '0h35m', '1h16m',
    #           '1h0m', '0h47m', '0h33m', '0h49m', '1h2m', '0h53m']
    # dec_lst = ['-73d53m', '-72d18m', '-72d25m', '-72d49m', '-72d36m',
    #            '-72d1m', '-72d59m', '-74d57m', '-75d30m', '-70d28m',
    #            '-75d44m', '-75d46m', '-76d46m']

    def claria_2005_dep_dist(rho, phi, theta, glx_incl):
        '''
        Deprojected distance from Claria et al. 2005. This formula is obtained
        from the vdM&C01 eqs:

        x = rho.cos(phi)
        y = rho.sin(phi)

        and the Cioni (2009) eqs assuming:

        p = phi & p' = theta (= PA - 90)

        or

        p = Phi & p' = PA (= theta + 90)
        '''
        A = 1 + (np.sin(phi.radian - theta.radian) * np.tan(glx_incl)) ** 2
        dep_dist_deg = Angle(rho * np.sqrt(A), unit='degree')

        return dep_dist_deg


    def phi_palma_func(coord, glx_ctr):
        '''
        '''
        # Angular separation between center and coordinates.
        cos_rho = np.cos(coord.dec.radian) * np.cos(glx_ctr.dec.radian) * \
            np.cos(coord.ra.radian - glx_ctr.ra.radian) + \
            np.sin(coord.dec.radian) * np.sin(glx_ctr.dec.radian)
        rho = Angle(np.arccos(cos_rho) * 180. / np.pi, unit=u.deg)

        # Position angle.
        cos_phi = (-1. * np.cos(coord.dec.radian) * np.sin(coord.ra.radian -
                   glx_ctr.ra.radian)) / np.sin(rho.radian)
        phi_palma = Angle(np.arccos(cos_phi) * 180. / np.pi, unit=u.deg)

        # sin_phi = (np.sin(coord.dec.radian) * np.cos(glx_ctr.dec.radian) -
        #            np.cos(coord.dec.radian) * np.sin(glx_ctr.dec.radian) *
        #            np.cos(coord.ra.radian - glx_ctr.ra.radian)) / \
        #     np.sin(rho.radian)
        # phi_palma = Angle(np.arcsin(sin_phi) * 180. / np.pi, unit=u.deg)

        Phi_palma = phi_palma - Angle('90d')

        return Phi_palma, phi_palma


    theta = gal_theta(pa)
    # Claria et al. 2005 equation. Equivalent to vdM&C01 + C09.
    dep_dist_deg = claria_2005_dep_dist(rho, phi, theta, inc)
    # print dep_dist_deg.degree

    # This gives the Palma et al. values. Notice the use of the galaxy's PA
    # with the 'phi_p' value.
    Phi_p, phi_p = phi_palma_func(ra_dec, cent)
    dep_dist_deg_p = claria_2005_dep_dist(rho, phi_p, pa, inc)
    # print dep_dist_deg_p.degree
