
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import Distance, Angle, SkyCoord
from astropy import units as u
import sys
sys.path.insert(0, '/media/rest/github/mc-catalog/functions')
from deproj_dist import deproj_dist
from lin_fit_conf_bands import linear_fit


"""
Data taken from Carrera et al. (2011); Table 1
"""
cent = SkyCoord('05h27.6m', '-69.87d', frame='icrs')
dist = Distance(51., unit=u.kpc)
ra_lst = ['05h12m', '05h14m', '05h13m', '04h22m', '05h13m', '06h49m', '05h12m',
          '03h56m', '05h11m', '07h14m']
dec_lst = ['-66d48m', '-65d03m', '-63d33m', '-71d00m', '-61d59m', '-71d10m',
           '-61d00m', '-71d00m', '-60d00m', '-71d00m']
inc, pa = Angle('34.7d'), Angle('189.3d')  # LMC
fe_h = np.array([-0.58, -0.64, -0.54, -0.43, -1.03, -0.74, -1.04, -0.87,
                 -0.99, -0.85])
sigma_feh = np.array([0.47, 0.54, 0.48, 0.39, 0.64, 0.62, 0.70, 0.56, 0.20,
                      0.35])

"""
Data taken from Carrera et al. (2008), "The Chemical Enrichment History of
the Small Magellanic Cloud and Its Gradients"

Inclination and position angle data taken from Cioni et al. (2009), "The
Chemical Enrichment History of the Small Magellanic Cloud and Its Gradients"
"""
# cent = SkyCoord('00h52m45s', '-72d49m43s', frame='icrs')
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

#
dist_kpc = []
for i, (ra, dec) in enumerate(zip(*[ra_lst, dec_lst])):
    a = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    d_kpc = deproj_dist(a, cent, pa, inc, dist)[1]
    print d_kpc
    dist_kpc.append(d_kpc.value)

print '\nNon-weighted fit.'
fit = np.polyfit(dist_kpc, fe_h, 1)
fit_nw = np.poly1d(fit)
print fit_nw
print '\nWeighted fit.'
a, b, sa, sb, rchi2, dof = linear_fit(np.array(dist_kpc), fe_h, sigma_feh)
print a, sa, b, sb, rchi2, dof

plt.scatter(dist_kpc, fe_h)
plt.show()
