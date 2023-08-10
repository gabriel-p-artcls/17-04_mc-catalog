
import numpy as np
import matplotlib.pyplot as plt


def num_memb_conc_param(x, x_t):
    # Obtain approximate number of members.
    n_c_k = np.log(1. + x) - 4. *\
        ((np.sqrt(1. + x) - 1.) / np.sqrt(1. + x_t)) + (x / (1. + x_t))
    return n_c_k


def king_tot_membs(f_0, r_c, x_t):
    # Obtain approximate number of members.
    t_m = np.pi * f_0 * (r_c ** 2) * (
        np.log(1. + x_t) - ((3. * np.sqrt(1. + x_t) - 1.) *
                            (np.sqrt(1. + x_t) - 1.)) / (1. + x_t))
    return t_m


f_0 = 9.
r_c, x_t = 18.8, 6.9 ** 2
t_m = king_tot_membs(f_0, r_c, x_t)
print('Total members Godfrooij: {}'.format(int(t_m)))

r_c, r_t = 18.8, 130.
c = r_t / r_c
x_t = c ** 2
max_d, bck = 2.52, 0.089
print((max_d - bck) / (1. - 1. / np.sqrt(1. + c ** 2)) ** 2)
xy = []
for _ in np.arange(0., 0.5, 0.05):
    f_0 = (max_d - bck) / (1. - 1. / np.sqrt(1. + c ** 2)) ** 2
    tot_memb = np.pi * (r_c ** 2) * f_0 * num_memb_conc_param(x_t, x_t)
    xy.append([max_d, tot_memb])
    print(f_0, int(tot_memb))
plt.plot(*zip(*xy))
plt.show()
