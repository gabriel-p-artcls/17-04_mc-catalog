
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
cent = SkyCoord('05h20m47s', '-69d28m41s', frame='icrs')
dist = Distance(50.12, unit=u.kpc)
inc, pa = Angle('35.8d'), Angle('145.d')
ra_lst = ['4h30m40s', '4h35m38s', '4h37m39s', '4h37m52s', '4h39m42s', '4h43m14s', '4h46m5s', '4h46m25s', '4h46m40s', '4h47m26s', '4h47m30s', '4h48m37s', '4h49m0s', '4h49m7s', '4h49m14s', '4h49m27s', '4h49m41s', '4h50m21s', '4h50m29s', '4h50m48s', '4h50m58s', '4h51m11s', '4h51m30s', '4h51m41s', '4h52m45s', '4h52m54s', '4h53m52s', '4h54m5s', '4h54m5s', '4h54m12s', '4h55m1s', '4h55m3s', '4h55m39s', '4h55m39s', '4h55m42s', '4h55m52s', '4h56m26s', '4h56m28s', '4h56m29s', '4h56m29s', '4h56m32s', '4h56m51s', '4h56m52s', '4h56m54s', '4h57m8s', '4h57m22s', '4h57m26s', '4h57m34s', '4h58m10s', '4h58m15s', '4h58m15s', '4h58m51s', '4h58m54s', '4h59m15s', '4h59m38s', '4h59m46s', '4h59m53s', '5h0m4s', '5h4m29s', '5h4m34s', '5h5m25s', '5h5m36s', '5h5m40s', '5h5m54s', '5h6m3s', '5h6m2s', '5h6m9s', '5h6m13s', '5h6m25s', '5h6m32s', '5h6m34s', '5h6m41s', '5h7m19s', '5h7m21s', '5h7m32s', '5h7m35s', '5h7m37s', '5h7m46s', '5h7m47s', '5h7m56s', '5h8m28s', '5h8m51s', '5h8m53s', '5h8m54s', '5h9m0s', '5h9m12s', '5h9m15s', '5h9m21s', '5h9m23s', '5h10m2s', '5h10m18s', '5h10m30s', '5h10m32s', '5h10m33s', '5h10m36s', '5h10m39s', '5h10m40s', '5h10m46s', '5h10m54s', '5h10m56s', '5h11m9s', '5h11m11s', '5h11m40s', '5h11m41s', '5h12m25s', '5h12m30s', '5h12m40s', '5h13m43s', '5h13m48s', '5h15m7s', '5h15m15s', '5h15m21s', '5h15m22s', '5h15m25s', '5h15m26s', '5h15m36s', '5h15m39s', '5h15m47s', '5h16m4s', '5h16m27s', '5h16m31s', '5h16m32s', '5h16m43s', '5h16m50s', '5h16m52s', '5h16m56s', '5h17m20s', '5h17m20s', '5h17m26s', '5h17m49s', '5h18m5s', '5h18m41s', '5h18m47s', '5h18m47s', '5h19m2s', '5h19m5s', '5h19m9s', '5h19m24s', '5h19m54s', '5h19m55s', '5h20m4s', '5h20m5s', '5h20m12s', '5h20m21s', '5h21m3s', '5h21m5s', '5h21m7s', '5h21m9s', '5h21m14s', '5h21m45s', '5h21m46s', '5h22m3s', '5h22m17s', '5h23m6s', '5h23m12s', '5h24m13s', '5h24m28s', '5h24m30s', '5h25m27s', '5h25m28s', '5h25m57s', '5h26m12s', '5h26m24s', '5h26m34s', '5h26m46s', '5h26m49s', '5h26m53s', '5h27m18s', '5h28m39s', '5h28m50s', '5h29m20s', '5h29m39s', '5h29m46s', '5h29m48s', '5h29m56s', '5h29m59s', '5h30m26s', '5h30m34s', '5h30m34s', '5h31m24s', '5h31m42s', '5h31m46s', '5h31m51s', '5h31m57s', '5h31m58s', '5h32m3s', '5h32m23s', '5h33m2s', '5h33m20s', '5h33m21s', '5h34m13s', '5h34m39s', '5h37m46s', '5h39m8s', '5h40m24s', '5h41m21s', '5h41m30s', '5h41m38s', '5h41m49s', '5h41m58s', '5h42m8s', '5h42m29s', '5h43m15s', '5h43m20s', '5h43m29s', '5h43m35s', '5h43m38s', '5h43m43s', '5h43m56s', '5h44m11s', '5h44m14s', '5h44m15s', '5h44m42s', '5h44m47s', '5h44m57s', '5h45m1s', '5h45m20s', '5h45m25s', '5h45m27s', '5h45m32s', '5h45m46s', '5h45m50s', '5h45m56s', '5h45m59s', '5h46m11s', '5h46m12s', '5h46m37s', '5h46m41s', '5h46m47s', '5h46m48s', '5h46m51s', '5h46m51s', '5h46m52s', '5h46m52s', '5h48m0s', '5h48m12s', '5h48m12s', '5h48m28s', '5h48m33s', '5h48m35s', '5h48m46s', '5h49m17s', '5h49m36s', '5h50m3s', '5h50m15s', '5h50m17s', '5h50m28s', '5h50m45s', '5h50m48s', '5h52m11s', '5h53m15s', '5h53m23s', '5h53m27s', '5h53m27s', '5h55m35s', '5h55m42s', '5h57m51s', '5h58m33s', '6h0m38s', '6h1m52s', '6h2m2s', '6h6m31s', '6h7m29s', '6h8m15s', '6h8m53s', '6h10m42s', '6h13m27s', '6h14m28s', '6h14m41s', '6h14m54s', '6h15m17s', '6h15m57s', '6h17m19s', '6h21m34s', '6h29m58s', '7h7m39s']
dec_lst = ['-66d57m25s', '-73d43m54s', '-66d11m58s', '-69d1m42s', '-74d1m2s', '-73d48m43s', '-66d54m41s', '-72d34m6s', '-67d41m7s', '-67d39m36s', '-72d35m18s', '-68d33m31s', '-72d38m24s', '-67d20m30s', '-72d3m24s', '-72d46m54s', '-72d14m50s', '-72d49m36s', '-67d19m36s', '-72d34m36s', '-67d36m36s', '-67d32m1s', '-67d27m15s', '-72d13m13s', '-72d31m5s', '-72d10m23s', '-69d34m14s', '-69d40m54s', '-69d45m30s', '-69d48m25s', '-67d42m51s', '-67d57m52s', '-67d43m34s', '-67d49m19s', '-67d46m54s', '-69d42m21s', '-67d56m19s', '-67d41m46s', '-67d37m22s', '-69d59m0s', '-69d58m54s', '-70d6m3s', '-68d0m20s', '-68d0m8s', '-70d6m42s', '-62d32m5s', '-67d41m7s', '-65d16m0s', '-68d3m37s', '-67d46m2s', '-68d2m57s', '-69d57m28s', '-67d50m49s', '-67d54m32s', '-69d33m22s', '-69d48m4s', '-67d55m25s', '-67d48m2s', '-68d20m55s', '-68d12m30s', '-68d30m4s', '-68d37m46s', '-68d38m12s', '-67d2m58s', '-68d37m37s', '-68d1m35s', '-68d26m45s', '-68d3m53s', '-68d22m22s', '-68d21m44s', '-68d25m38s', '-67d50m28s', '-68d20m54s', '-66d49m45s', '-67d34m13s', '-67d27m39s', '-68d32m31s', '-67d51m41s', '-66d47m53s', '-67d21m28s', '-66d46m14s', '-67d58m49s', '-68d5m1s', '-66d47m8s', '-67d59m0s', '-68d26m39s', '-67d42m0s', '-62d22m46s', '-67d46m42s', '-66d42m0s', '-67d51m0s', '-68d24m2s', '-66d56m24s', '-67d7m39s', '-70d29m15s', '-66d43m45s', '-68d45m13s', '-67d29m6s', '-67d28m16s', '-67d37m36s', '-67d40m57s', '-67d37m37s', '-68d43m36s', '-67d33m56s', '-68d46m19s', '-67d17m28s', '-67d37m24s', '-67d24m10s', '-66d37m12s', '-68d58m43s', '-68d52m57s', '-69d6m27s', '-69d2m32s', '-68d40m52s', '-69d3m2s', '-69d8m21s', '-68d54m31s', '-69d14m39s', '-69d6m9s', '-69d4m49s', '-69d10m58s', '-68d55m7s', '-69d12m13s', '-69d3m35s', '-69d4m13s', '-68d40m58s', '-69d9m25s', '-69d12m49s', '-69d6m55s', '-68d28m22s', '-69d10m18s', '-69d4m46s', '-69d13m32s', '-69d16m37s', '-69d0m4s', '-68d52m14s', '-69d15m36s', '-68d52m52s', '-68d57m53s', '-68d48m7s', '-69d15m55s', '-63d28m49s', '-68d54m15s', '-69d14m48s', '-69d5m51s', '-69d4m16s', '-69d8m9s', '-69d7m2s', '-68d47m0s', '-68d55m2s', '-68d43m53s', '-70d2m44s', '-70d2m0s', '-75d26m48s', '-70d46m40s', '-75d34m0s', '-67d43m43s', '-67d40m41s', '-73d34m13s', '-69d46m32s', '-69d45m4s', '-70d58m53s', '-69d43m51s', '-69d50m27s', '-69d51m3s', '-69d50m17s', '-69d48m54s', '-73d40m48s', '-73d37m49s', '-71d37m58s', '-70d34m46s', '-70d59m2s', '-71d0m2s', '-63d38m58s', '-72d3m17s', '-67d52m44s', '-75d20m57s', '-63d12m12s', '-68d9m27s', '-72d2m33s', '-72d8m46s', '-68d14m8s', '-67d59m28s', '-67d52m43s', '-67d58m18s', '-64d14m32s', '-67d59m49s', '-67d50m56s', '-68d9m8s', '-75d22m35s', '-67d51m23s', '-68d18m20s', '-74d46m58s', '-74d51m12s', '-69d15m10s', '-69d3m46s', '-69d11m6s', '-69d18m48s', '-68d55m15s', '-69d2m51s', '-69d22m0s', '-65d21m46s', '-69d2m3s', '-66d15m44s', '-69d9m44s', '-66d12m31s', '-69d15m51s', '-69d13m23s', '-69d10m50s', '-69d20m0s', '-70d39m20s', '-70d40m10s', '-70d25m31s', '-70d24m22s', '-70d19m59s', '-70d32m34s', '-70d36m6s', '-70d24m5s', '-69d20m43s', '-70d45m34s', '-70d43m9s', '-69d22m49s', '-69d16m19s', '-70d43m46s', '-70d43m12s', '-69d4m57s', '-70d46m33s', '-70d50m52s', '-70d49m58s', '-70d35m21s', '-69d25m11s', '-70d30m40s', '-70d48m21s', '-69d11m23s', '-70d28m30s', '-70d28m0s', '-70d33m24s', '-70d32m52s', '-70d29m0s', '-70d18m39s', '-70d28m23s', '-70d47m54s', '-70d41m35s', '-67d43m5s', '-70d25m40s', '-70d36m56s', '-70d32m33s', '-70d34m34s', '-71d42m28s', '-71d51m30s', '-71d53m32s', '-70d4m16s', '-71d41m10s', '-71d42m57s', '-62d20m43s', '-74d21m14s', '-66d24m2s', '-65d28m37s', '-70d4m10s', '-72d21m19s', '-60d31m24s', '-72d13m35s', '-72d29m39s', '-62d59m15s', '-72d23m2s', '-71d31m44s', '-70d41m45s', '-72d36m34s', '-69d48m7s', '-72d30m19s', '-73d47m7s', '-70d4m23s', '-70d3m39s', '-72d47m24s', '-69d20m1s', '-69d59m2s']
fe_h = [0.]*len(ra_lst)
sigma_feh = [0.]*len(ra_lst)

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
