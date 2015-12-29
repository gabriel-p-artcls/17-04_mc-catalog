
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
import numpy as np
try:
    from functions.MCs_data import MCs_data
except:
    from MCs_data import MCs_data


def dist_err_mag_2_pc(d_pc, e_dm, e_E):
    '''
    Error associated to the distance in pc, obtained via the distance modulus
    relation.

    d_pc = 10 ** (0.2*(dm + 5 - 3.1*E))

    e_d_pc = sqrt((diff(d,dm)*e_dm)**2 + (diff(d,E)*e_E)**2)

    where: dm is distance modulus, E is E(B-V) extinction, diff(x,y) is the
    partial derivative of x with respect to y, and e_X is the standard
    deviation of parameter X.
    '''

    diff_dm = d_pc * np.log(10) * 0.2
    diff_E = d_pc * np.log(10) * 0.2 * (-3.1)
    e_d_pc = np.sqrt((diff_dm * e_dm) ** 2 + (diff_E * e_E) ** 2)

    return e_d_pc


def dist_err_2_pts(d_pc, c1, e_r_0, c2, e_r):
    '''
    Error associated to the distance between two points in spherical
    coordinates:

    c1 = (r_0, ra_0, dec_0)
    c2 = (r, ra, dec)

    Distance between two points in spherical coordinates
    (http://math.stackexchange.com/a/833110/37846):

    A = cos(dec)*cos(dec_0)*cos(ra-ra_0)+sin(dec)*sin(dec_0)
    d_pc = sqrt(r**2 + r_0**2 - 2*r*r_0*A)

    Error in distance between points, assuming only the distances of the
    points contain error, not the alpha,delta coordinates.

    e_d_pc = (1/d_pc) * sqrt([(r-r_0*A)e_r]**2 + [(r_0-r*A)e_r_0]**2)

    where:

    e_r, e_r_0

    are the errors in the distance for each point.
    '''

    ra_0, dec_0, r_0 = c1.ra.radian, c1.dec.radian, c1.distance
    ra, dec, r = c2.ra.radian, c2.dec.radian, c2.distance

    alpha_delta_par = np.cos(dec) * np.cos(dec_0) * np.cos(ra - ra_0) +\
        np.sin(dec) * np.sin(dec_0)

    A = np.sqrt(((r - r_0 * alpha_delta_par) * e_r) ** 2 +
                ((r_0 - r * alpha_delta_par) * e_r_0) ** 2)
    e_d_pc = A / d_pc

    return e_d_pc


def dist_2_cloud_center(gal, ra_deg, dec_deg, dist_mod, e_dm):
    '''
    Obtain the 3D distance in parsecs between the center of a cluster and
    the center of the Magellanic cloud it belongs to.
    '''

    # Retrieve the galaxy's center coords and distance modulus (and its error).
    gal_center, gal_dist, gal_e_dm = MCs_data(gal)
    # Obtain error in distance in parsecs. Use 0. error in extinction since
    # this is the true distance modulus.
    e_gal_dist_pc = dist_err_mag_2_pc(gal_dist, gal_e_dm, 0.)

    # *Individual* distance (ASteCA) for each cluster (in parsecs).
    # ASteCA gives the distance modulus as dm = -5 + 5*log(d), so to transform
    # that into the distance in parsecs 'd', we do:
    d_clust = Distance(10 ** (0.2 * (float(dist_mod) + 5)), unit=u.pc)
    # Obtain error in distance in parsecs.
    e_cl_dist_pc = dist_err_mag_2_pc(d_clust, e_dm, 0.)

    # Galaxy center coordinate.
    c1 = SkyCoord(ra=gal_center.ra, dec=gal_center.dec, distance=gal_dist)
    # Cluster coordinate.
    c2 = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, distance=d_clust)

    # 3D distance between cluster and center of galaxy,  in parsecs.
    dist_pc = c1.separation_3d(c2)
    # Error for the 3D distance.
    e_d_pc = dist_err_2_pts(dist_pc, c1, e_gal_dist_pc, c2, e_cl_dist_pc)
    print dist_pc, e_d_pc

    return dist_pc


if __name__ == "__main__":
    dist_2_cloud_center(0, 12.75, -72.5, 18.52, 0.02)
