
from astropy.coordinates import Distance, Angle, SkyCoord
from astropy import units as u
import numpy as np


def rho_phi(coord, glx_ctr):
    '''
    Eqs 1, 2 & 3 from van der Marel & Cioni (2001).
    '''
    # Angular distance between point and center of galaxy.
    rho = coord.separation(glx_ctr)

    # Position angle between center and coordinates. This is the angle between
    # the positive y axis (North) counter-clockwise towards the negative x
    # axis (East).
    Phi = glx_ctr.position_angle(coord)
    phi = Phi + Angle('90d')

    # x = rho.degree * np.cos(phi)
    # y = rho.degree * np.sin(phi)
    # if x >= 0. and y >= 0.:
    #     print 'NW', phi.degree, phi_palma.degree
    # elif x <= 0. and y >= 0.:
    #     print 'NE', phi.degree, phi_palma.degree
    # elif x <= 0. and y <= 0.:
    #     print 'SE', phi.degree, phi_palma.degree
    # elif x >= 0. and y <= 0.:
    #     print 'SW', phi.degree, phi_palma.degree

    return rho, Phi, phi


def phi_palma(coord, glx_ctr):
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
    # sin_phi = (np.sin(coord.dec.radian) * np.cos(glx_ctr.dec.radian) -
    #            np.cos(coord.dec.radian) * np.sin(glx_ctr.dec.radian) *
    #            np.cos(coord.ra.radian - glx_ctr.ra.radian)) / \
    #     np.sin(rho.radian)
    phi_palma = Angle(np.arccos(cos_phi) * 180. / np.pi, unit=u.deg)
    Phi_palma = phi_palma - Angle('90d')

    return Phi_palma, phi_palma


def gal_theta(glx_PA):
    '''
    PA of the galaxy, rotated 90 deg.
    '''
    # theta = glx_PA - Angle('90d')  # Acc to Cioni (2009)
    theta = glx_PA + Angle('90d')  # Acc to vdM et al. (2002)
    # ???????

    return theta


def vdm_2001_D(glx_incl, D_0, rho, phi, theta):
    '''
    Eq 8 from van der Marel & Cioni (2001).
    '''
    # Distance to the coordinates passed.
    A = np.cos(glx_incl.radian) * np.cos(rho.radian) - \
        np.sin(glx_incl.radian) * np.sin(rho.radian) * \
        np.sin(phi.radian - theta.radian)
    D = D_0 * np.cos(glx_incl.radian) / (A)

    return D


def vdm_2001_xy(rho, phi):
    '''
    Eq 4 from van der Marel & Cioni (2001).
    '''
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


def vdm_2001_dep_dist(rho, phi, theta, glx_incl, D, D_0):
    '''
    Deprojected angular distance from vdM & Cioni (2001).
    '''
    # Eq 7 from van der Marel & Cioni (2001).
    x_p = D * np.sin(rho.radian) * np.cos(phi.radian - theta.radian)
    y_p = D * (np.sin(rho.radian) * np.cos(glx_incl.radian) *
               np.sin(phi.radian - theta.radian) + np.cos(rho.radian) *
               np.sin(glx_incl.radian)) - D_0 * np.sin(glx_incl.radian)
    # z_p = 0
    # z_p = D * (np.sin(rho.radian) * np.sin(glx_incl.radian) *
    #            np.sin(phi.radian - theta.radian) - np.cos(rho.radian) *
    #            np.cos(glx_incl.radian)) + D_0 * np.cos(glx_incl.radian)
    d_kpc = np.sqrt(x_p ** 2 + y_p ** 2)

    # # The above is equivalent to using the cosine law.
    # d_kpc = np.sqrt(D_0**2+D**2-2*D_0*D*np.cos(rho.radian))

    return d_kpc


def claria_2005_dep_dist(rho, phi, glx_incl, theta):
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
    rho_p = np.sqrt(A)
    dep_dist_deg = Angle(rho * rho_p, unit='degree')

    return rho_p, dep_dist_deg


def cioni_2009_dep_dist(glx_incl, theta, x, y):
    '''
    Deprojected angular distance. Eqs 1, 2, 3 & 4 from Cioni (2009).
    '''
    # Rotate the coords system.
    x1 = x * np.cos(theta.radian) + y * np.sin(theta.radian)
    y1 = y * np.cos(theta.radian) - x * np.sin(theta.radian)
    # De-project.
    y2 = y1 / np.cos(glx_incl.radian)
    # Obtain de-projected distance in decimal degrees.
    dep_dist_deg = Angle(np.sqrt(x1 ** 2 + y2 ** 2), unit='degree')

    return dep_dist_deg


def cioni_2009_dist_kpc(dep_dist_deg, D_0):
    '''
    Obtain de-projected distance in the units used for the galaxy center
    distance.
    Eq 5 from Cioni (2009).
    '''
    dep_dist_kpc = Distance(np.tan(dep_dist_deg) * D_0, unit=D_0.unit)

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

    D_0 = glx_dist
    D = vdm_2001_D(glx_incl, D_0, rho, phi, theta)
    dep_dist_kpc_M01 = vdm_2001_dep_dist(rho, phi, theta, glx_incl, D, D_0)
    print 'd_kpc M01 = {:.5f}'.format(dep_dist_kpc_M01)
    import pdb; pdb.set_trace()  # breakpoint 7a6db5c5 //

    # This gives the values for the deprojected angular distance in
    # Carrera et al. (2011).
    # x, y = carrera_2011_xy(rho, Phi)

    # # vdm&C01 (x,y values) + Cioni 2009 (C09).
    # x, y = vdm_2001_xy(rho, phi)
    # # theta = glx_PA - Angle('90d')  # Acc to Cioni (2009), but NOT to vdM&C01
    # dep_dist_deg = cioni_2009_dep_dist(glx_incl, theta, x, y)
    # print 'd_deg C09 = {:.5f}'.format(dep_dist_deg.degree)

    # Claria et al. 2005 equation. Equivalent to vdM&C01 + C09.
    rho_p, dep_dist_deg = claria_2005_dep_dist(rho, Phi, glx_incl, glx_PA)

    # TEST
    X = D_0 * np.tan(rho) * rho_p
    print 'd_kpc XXX = {:.5f}'.format(X)
    # TEST

    # # This gives the Palma et al. values. Notice the use of the galaxy's PA
    # # with the 'phi_p' value.
    # Phi_p, phi_p = phi_palma(coord, glx_ctr)
    # dep_dist_deg = claria_2005_dep_dist(rho, phi_p, glx_incl, glx_PA)

    # print 'd_deg C05 = {:.5f}'.format(dep_dist_deg.degree)
    dep_dist_kpc = cioni_2009_dist_kpc(dep_dist_deg, D_0)
    print 'd_kpc C09 = {:.5f}'.format(dep_dist_kpc)

    dep_dist_kpc = Angle(0., unit='degree')
    print ''
    return dep_dist_kpc


def main():
    # Values in all 4 quadrants.
    cent = SkyCoord('05h27.6m', '-69.87d', frame='icrs')
    dist = Distance(51., unit=u.kpc)
    inc, pa = Angle('34.7d'), Angle('189.3d')  # LMC
    ra_lst = ['05h20m', '05h40m', '05h40m', '05h20m']
    dec_lst = ['-68d', '-68d', '-70d', '-70d']

    # # Palma et al. 2015
    # cent = SkyCoord('05h20m47s', '-69d28m41s', frame='icrs')
    # dist = Distance(50.12, unit=u.kpc)
    # inc, pa = Angle('35.8d'), Angle('145.d')
    # ra_lst = ['7h7m39s']
    # dec_lst = ['-69d59m2s']

    # # Random values.
    # ra_lst, dec_lst = [], []
    # N = 25
    # for _ in range(N):
    #     h = int(np.random.uniform(4., 7., 1)[0])
    #     m = np.random.uniform(0., 59., 1)[0]
    #     d = np.random.uniform(-75., -65., 1)[0]
    #     ra_lst.append(str(h) + 'h' + str(m) + 'm')
    #     dec_lst.append(str(d) + 'd')
    # print ra_lst
    # print dec_lst, '\n'

    # # Carrera et al 2011, first two sources.
    # # cent = SkyCoord('05h27.6m', '-69.87d', frame='icrs')
    # cent = SkyCoord('82.25d', '-69.5d', frame='icrs')
    # dist = Distance(51., unit=u.kpc)
    # ra_lst = ['05h12m', '05h14m']
    # dec_lst = ['-66d48m', '-65d03m']
    # inc, pa = Angle('34.7d'), Angle('189.3d')

    for i, (ra, dec) in enumerate(zip(*[ra_lst, dec_lst])):
        a = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
        deproj_dist(a, cent, pa, inc, dist)

if __name__ == "__main__":
    main()
