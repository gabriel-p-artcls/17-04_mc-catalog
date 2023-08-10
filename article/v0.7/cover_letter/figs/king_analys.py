
import numpy as np
from scipy import integrate as integrate
import matplotlib.pyplot as plt


def kprof(x, rc, rt):
    '''
    Standard King profile, minus the central (and field) density.
    '''
    return (1. / np.sqrt(1. + (x / rc) ** 2) -
            1. / np.sqrt(1. + (rt / rc) ** 2)) ** 2


def kprof_crowded(x, cr, rc, rt):
    '''
    King profile affected by a crowding photometry function, minus the central
    (and field) density.
    '''
    crow_f = cr + x * (1. - cr) / rt
    return crow_f *\
        (1. / np.sqrt(1. + (x / rc) ** 2) -
            1. / np.sqrt(1. + (rt / rc) ** 2)) ** 2


def kprof_int(x, rc, rt):
    '''
    Standard King profile with a 2*Pi*r term added, for integration (as stated
    in King 1962), minus the central density term (k) that is assumed constant.
    '''
    return 2. * np.pi * x *\
        (1. / np.sqrt(1. + (x / rc) ** 2) -
            1. / np.sqrt(1. + (rt / rc) ** 2)) ** 2


def kprof_crowded_int(x, cr, rc, rt):
    '''
    King profile affected by a crowding photometry function., and with the
    2*Pi*r term added for integration, minus the central density term (k)
    that is assumed constant.
    '''
    crow_f = cr + x * (1. - cr) / rt
    return crow_f * 2. * np.pi * x *\
        (1. / np.sqrt(1. + (x / rc) ** 2) -
            1. / np.sqrt(1. + (rt / rc) ** 2)) ** 2


def num_memb_conc_param(x, x_t):
    '''
    Number of members according to a King profile, minus a term of Pi*k*rc**2,
    that is assumed constant.
    '''
    return np.log(1. + x) - 4. *\
        ((np.sqrt(1. + x) - 1.) / np.sqrt(1. + x_t)) + (x / (1. + x_t))


# Some rc, rt values for testing.
rc, rt = 18.8, 130.
x_t = (rt / rc) ** 2
r = np.arange(1., 130., 1.)
x = (r / rc) ** 2
# This parameter controls how the crowding function affects the King profile.
cr = 0.5

# Standards King profile.
plt.plot(r, kprof(r, rc, rt), c='b')
# King profile affected by crowding.
plt.plot(r, kprof_crowded(r, cr, rc, rt), c='r')
plt.show()

# Number of members for radius r using King's equation (18).
y0 = np.pi * (rc ** 2) * num_memb_conc_param(x, x_t)
# Number of members for radius r integrating King's density profile. This
# should be equal to the y0 curve.
y1 = []
for _ in r:
    y1.append(integrate.quad(kprof_int, 0., _, args=(rc, rt))[0])
# Number of members for radius r integrating King's density profile affected
# by crowding.
y2 = []
for _ in r:
    y2.append(integrate.quad(kprof_crowded_int, 0., _, args=(cr, rc, rt))[0])

membs_t = np.pi * (rc ** 2) * num_memb_conc_param(x_t, x_t)
membs_t_crow = integrate.quad(kprof_crowded_int, 0., rt, args=(cr, rc, rt))[0]
print('Tot members: {:.0f}'.format(membs_t))
print('Tot members crowded: {:.0f}'.format(membs_t_crow))
plt.plot(np.sqrt(x) * rc, y0, c='r')
plt.plot(np.sqrt(x) * rc, y1, c='b')
plt.plot(np.sqrt(x) * rc, y2, c='g')
plt.show()

# Fraction of cluster members lost using a radius that is smaller than rt.
rad_f = [0.3, 0.4, 0.5, 0.6, 0.7]
curves1, curves2 = [], []
for frac in rad_f:
    pts1, pts2 = [[], []], [[], []]
    for x_t in np.arange(1., 100., 5.):
        # x = (r/r_c)**2
        # x_t = (r_t/r_c)**2
        # Total number of members up to the tidal radius.
        x = x_t
        a = num_memb_conc_param(x, x_t)
        x = x_t * (frac ** 2)
        # Total number of members up to a fraction of the tidal radius.
        b = num_memb_conc_param(x, x_t)
        pts1[0].append(np.sqrt(x_t))
        pts1[1].append(1. - b / a)

        # Same as above, except the profile is here affected by crowding.
        rc = 1.
        rt = np.sqrt(x_t) * rc
        a = integrate.quad(kprof_crowded_int, 0., rt, args=(cr, rc, rt))[0]
        b = integrate.quad(
            kprof_crowded_int, 0., rt * frac, args=(cr, rc, rt))[0]
        pts2[0].append(np.sqrt(x_t))
        pts2[1].append(1. - b / a)

    curves1.append(pts1)
    curves2.append(pts2)

cols = ['b', 'g', 'r', 'c', 'm']
for i, c in enumerate(curves1):
    plt.plot(c[0], c[1], label='r = {} r_t'.format(rad_f[i]), lw=2.5,
             color=cols[i])
for i, c in enumerate(curves2):
    plt.plot(c[0], c[1], ls='--', lw=1.5, color=cols[i])
plt.legend()
plt.grid(b=True, which='major', color='gray', linestyle='--', lw=1)
plt.xlabel('c = (r_t / r_c)')
plt.ylabel('Fraction of members lost')
plt.title('Cluster members lost as a function of radius and concentration')
plt.xlim(1., 9.)
plt.ylim(0., 0.7)
plt.show()
