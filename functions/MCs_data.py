
from astropy.coordinates import Distance, SkyCoord
from astropy import units as u


def MCs_data(gal):
    '''
    Store center and distance data for the MCs.

    S/LMC distance stored in parsecs.
    de Grijs et al. (2015): SMC (m-M)o = 18.96 +/- 0.02 mag
    de Grijs et al. (2014): LMC (m-M)o = 18.49 +/- 0.09 mag

    '''

    # SMC central coords stored in degrees.
    c_SMC = SkyCoord('00h52m45s', '-72d49m43s', frame='icrs')
    # ^ (13.1875, -72.82861111)
    # d_SMC ~ 61944.11 pc (18.96 mag)
    d_SMC = Distance(10 ** (0.2 * (18.96 + 5)), unit=u.pc)
    e_dm_SMC = 0.02

    # LMC central coords stored in degrees.
    c_LMC = SkyCoord('05h20m57s', '-69d28m41s', frame='icrs')
    # ^ (80.2375, -69.47805556)
    # d_LMC ~ 49888.45 pc (18.49 mag)
    d_LMC = Distance(10 ** (0.2 * (18.49 + 5)), unit=u.pc)
    e_dm_LMC = 0.09

    if gal == 0:  # SMC
        gal_center, gal_dist, e_gal_dist = c_SMC, d_SMC, e_dm_SMC
    else:  # LMC
        gal_center, gal_dist, e_gal_dist = c_LMC, d_LMC, e_dm_LMC

    return gal_center, gal_dist, e_gal_dist
