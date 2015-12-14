
import os
import numpy as np
from astropy.coordinates import Distance, Angle, SkyCoord
from astropy import units as u
import sys
# Define main path.
r_path = os.path.realpath(__file__)[:-57]
sys.path.insert(0, r_path + 'github/mc-catalog/functions/')
from deproj_dist import deproj_dist
from lin_fit_conf_bands import linear_fit

"""
Data taken from Carrera et al. (2011); Table 1
"""
# cent = SkyCoord('05h27.6m', '-69.87d', frame='icrs')
# dist = Distance(51., unit=u.kpc)
# ra_lst = ['05h12m', '05h14m', '05h13m', '04h22m', '05h13m', '06h49m',
#           '05h12m', '03h56m', '05h11m', '07h14m']
# dec_lst = ['-66d48m', '-65d03m', '-63d33m', '-71d00m', '-61d59m', '-71d10m',
#            '-61d00m', '-71d00m', '-60d00m', '-71d00m']
# inc, pa = Angle('34.7d'), Angle('189.3d')  # LMC
# fe_h = np.array([-0.58, -0.64, -0.54, -0.43, -1.03, -0.74, -1.04, -0.87,
#                  -0.99, -0.85])
# sigma_feh = np.array([0.47, 0.54, 0.48, 0.39, 0.64, 0.62, 0.70, 0.56, 0.20,
#                       0.35])

"""
Data taken from Carrera et al. (2008), "The Chemical Enrichment History of
the Small Magellanic Cloud and Its Gradients"

Inclination and position angle data taken from Cioni et al. (2009), "The
Chemical Enrichment History of the Small Magellanic Cloud and Its Gradients"

Center taken from Noel et al. (2007)
(https://ui.adsabs.harvard.edu/#abs/2007AJ....133.2037N/abstract) which is
'Paper I' according to Carrera et al. (2008).
"""
# # cent = SkyCoord('00h52m45s', '-72d49m43s', frame='icrs')
# cent = SkyCoord('00h52m42s', '-72d49m', frame='icrs')
# dist = Distance(60.26, unit=u.kpc)
# inc, pa = Angle('65.5d'), Angle('45.d')
# ra_lst = ['00h57m', '00h37m', '00h36m', '01h11m', '01h12m', '00h35m', '01h16m',
#           '01h00m', '00h47m', '00h33m', '00h49m', '01h02m', '00h53m']
# dec_lst = ['-73d53m', '-72d18m', '-72d25m', '-72d49m', '-72d36m', '-72d01m',
#            '-72d59m', '-74d57m', '-75d30m', '-70d28m', '-75d44m', '-75d46m',
#            '-76d46m']
# fe_h = np.array([-1.01, -0.95, -0.98, -1.08, -1.16, -1.09, -0.96, -1.07,
#                  -1.15, -1.58, -1.00, -1.29, -1.64])
# sigma_feh = np.array([0.33, 0.17, 0.25, 0.21, 0.32, 0.24, 0.26, 0.28, 0.27,
#                       0.57, 0.28, 0.42, 0.50])

'''
Data taken from Palma et al. (2015),

'''
# cent = SkyCoord('05h20m47s', '-69d28m41s', frame='icrs')
# dist = Distance(50.12, unit=u.kpc)
# inc, pa = Angle('35.8d'), Angle('145.d')
# ra_lst = ['4h30m40s', '4h35m38s', '4h37m39s', '4h37m52s', '4h39m42s',
# '4h43m14s', '4h46m5s', '4h46m25s', '4h46m40s', '4h47m26s', '4h47m30s',
# '4h48m37s', '4h49m0s', '4h49m7s', '4h49m14s', '4h49m27s', '4h49m41s',
# '4h50m21s', '4h50m29s', '4h50m48s', '4h50m58s', '4h51m11s', '4h51m30s',
# '4h51m41s', '4h52m45s', '4h52m54s', '4h53m52s', '4h54m5s', '4h54m5s',
# '4h54m12s', '4h55m1s', '4h55m3s', '4h55m39s', '4h55m39s', '4h55m42s',
# '4h55m52s', '4h56m26s', '4h56m28s', '4h56m29s', '4h56m29s', '4h56m32s',
# '4h56m51s', '4h56m52s', '4h56m54s', '4h57m8s', '4h57m22s', '4h57m26s',
# '4h57m34s', '4h58m10s', '4h58m15s', '4h58m15s', '4h58m51s', '4h58m54s',
# '4h59m15s', '4h59m38s', '4h59m46s', '4h59m53s', '5h0m4s', '5h4m29s', '5h4m34s']
# dec_lst = ['-66d57m25s', '-73d43m54s', '-66d11m58s', '-69d1m42s', '-74d1m2s',
# '-73d48m43s', '-66d54m41s', '-72d34m6s', '-67d41m7s', '-67d39m36s',
# '-72d35m18s', '-68d33m31s', '-72d38m24s', '-67d20m30s', '-72d3m24s',
# '-72d46m54s', '-72d14m50s', '-72d49m36s', '-67d19m36s', '-72d34m36s',
# '-67d36m36s', '-67d32m1s', '-67d27m15s', '-72d13m13s', '-72d31m5s',
# '-72d10m23s', '-69d34m14s', '-69d40m54s', '-69d45m30s', '-69d48m25s',
# '-67d42m51s', '-67d57m52s', '-67d43m34s', '-67d49m19s', '-67d46m54s',
# '-69d42m21s', '-67d56m19s', '-67d41m46s', '-67d37m22s', '-69d59m0s',
# '-69d58m54s', '-70d6m3s', '-68d0m20s', '-68d0m8s', '-70d6m42s', '-62d32m5s',
# '-67d41m7s', '-65d16m0s', '-68d3m37s', '-67d46m2s', '-68d2m57s', '-69d57m28s',
# '-67d50m49s', '-67d54m32s', '-69d33m22s', '-69d48m4s', '-67d55m25s',
# '-67d48m2s', '-68d20m55s', '-68d12m30s']
# fe_h = [0.]*len(ra_lst)
# sigma_feh = [0.]*len(ra_lst)

#
dist_kpc = []
for i, (ra, dec) in enumerate(zip(*[ra_lst, dec_lst])):
    a = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    d_kpc = deproj_dist(a, cent, pa, inc, dist)
    dist_kpc.append(d_kpc.value)

print '\nNon-weighted fit.'
fit = np.polyfit(dist_kpc, fe_h, 1)
fit_nw = np.poly1d(fit)
print fit_nw
print '\nWeighted fit.'
a, b, sa, sb, rchi2, dof = linear_fit(np.array(dist_kpc), fe_h, sigma_feh)
print a, sa, b, sb, rchi2, dof

# import matplotlib.pyplot as plt
# plt.scatter(dist_kpc, fe_h)
# plt.show()
