
from os import walk
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.offsetbox as offsetbox
import numpy as np


def get_data(to_min, to_max, c_min, c_max):
    '''
    '''

    # Path to data file.
    path = '../OCs_data/MC_all/'

    # Go through each cluster data file.
    t1_e, c_e, jt, kt, jc, kc = [], [], 0., 0., 0., 0.
    for root, dirs, files in walk(path):
        for f in files:
            with open(root + f) as f1:
                # if f == 'SL551.OUT':
                #     break
                print f
                for line in f1:
                    if not line.strip().startswith('#'):
                        l = line.split()
                        if l:
                            if to_min <= float(l[3]) <= to_max:
                                t1_e.append([float(l[3]), float(l[4])])
                                jt += 1
                                if float(l[4]) > 0.1:
                                    kt += 1
                                    # print l
                            if c_min <= float(l[5]) <= c_max:
                                c_e.append([float(l[5]), float(l[6])])
                                jc += 1
                                if float(l[6]) > 0.1:
                                    kc += 1

    print jt, kt, jc, kc
    perc = [jt, kt, jc, kc]
    return t1_e, c_e, perc


def make_plot(minv, maxv, data, perc, avrg, e_avrg, id_tc):
    """
    """
    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(1, 1)

    print('Mean {} error in [{}, {}]: {}'.format(id_tc, minv, maxv,
          np.mean(zip(*data)[1])))

    ax = plt.subplot(gs[0])
    plt.xlabel(r'${}\, (mag)$'.format(id_tc))
    plt.ylabel(r'$\sigma_{{{}}}\, (mag)$'.format(id_tc))
    plt.xlim(minv, maxv)
    plt.ylim(-0.005, 0.5)
    plt.minorticks_on()
    plt.axhline(y=np.mean(zip(*data)[1]), lw=3, ls='--', color='r')
    plt.scatter(zip(*data)[0], zip(*data)[1], lw=0.5)
    # Average per 0.5 magnitude interval.
    plt.scatter(avrg, e_avrg, c='g', s=50)
    # for _ in np.arange(minv, maxv, 0.5):
    # Text box.
    text1 = r'$N={}$'.format(int(perc[0]))
    text2 = r'$\sigma_{{{}}}>0.1\;mag:\;{:.2f}\%$'.format(
        id_tc, perc[1] / perc[0] * 100.)
    text3 = r'$\overline{{\sigma_{{{}}}}}={:.4f}\;mag$'.format(
        id_tc, np.mean(zip(*data)[1]))
    text = text1 + '\n' + text2 + '\n' + text3
    ob = offsetbox.AnchoredText(text, loc=9, prop=dict(size=12))
    ob.patch.set(alpha=0.95)
    ax.add_artist(ob)

    # Output png file.
    fig.tight_layout()
    plt.savefig('{}_photom_errors.png'.format(id_tc), dpi=300,
                bbox_inches='tight')


def average_errors(data, minv, maxv, step):
    """
    """

    lims = np.arange(minv, maxv + step, step)
    avrg, e_avrg = [], []
    for _ in np.arange(minv + step / 2., maxv, step):
        avrg.append(_)

    for i, mi in enumerate(lims[:-1]):
        e_int = []
        for st in data:
            if mi <= st[0] < lims[i + 1]:
                e_int.append(st[1])
        print mi
        e_avrg.append(np.mean(e_int))

    return avrg, e_avrg


def main():
    '''
    Obtain photometric error for the T1 magnitude, and (C-T1) color.
    '''

    to_min, to_max = 17., 21.
    c_min, c_max = -0.5, 0.5

    # Read data from ASteCA output file.
    t1_e, c_e, perc = get_data(to_min, to_max, c_min, c_max)

    t1_avrg, e_t1_avrg = average_errors(t1_e, to_min, to_max, 0.5)
    # Plot for T1
    make_plot(to_min, to_max, t1_e, [perc[0], perc[1]], t1_avrg, e_t1_avrg,
              'T_1')

    c_avrg, e_c_avrg = average_errors(c_e, c_min, c_max, 0.1)
    # Plot for (C-T1)
    make_plot(c_min, c_max, c_e, [perc[2], perc[3]], c_avrg, e_c_avrg,
              '(C-T_1)')


if __name__ == "__main__":
    main()
