
from os import walk
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.offsetbox as offsetbox
import numpy as np


def get_t1_data():
    '''
    '''

    # Path to data file.
    path = '../OCs_data/MC_all/'

    # Go through each cluster data file.
    t1_e, j, k = [], 0., 0.
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
                            if 18.4 <= float(l[3]) <= 19.1:
                                t1_e.append([float(l[3]), float(l[4])])
                                j += 1
                                if float(l[4]) > 0.1:
                                    k += 1
                                    # print l

    print j, k
    perc = [j, k]
    return t1_e, perc


def make_plot(t1_e, perc):
    """
    """
    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(1, 1)

    print 'Mean T1 error in [18.4, 19.1]:', np.mean(zip(*t1_e)[1])

    ax = plt.subplot(gs[0])
    plt.xlabel(r'$T_1\, (mag)$')
    plt.ylabel(r'$\sigma_{T_1}\, (mag)$')
    plt.xlim(18.37, 19.13)
    plt.ylim(-0.005, 0.5)
    plt.minorticks_on()
    plt.axhline(y=np.mean(zip(*t1_e)[1]), lw=3, ls='--', color='r')
    plt.scatter(zip(*t1_e)[0], zip(*t1_e)[1], lw=0.5)
    # Text box.
    text1 = r'$N={}$'.format(int(perc[0]))
    text2 = r'$\sigma_{{T_1}}>0.1\;mag:\;{:.2f}\%$'.format(
        perc[1] / perc[0] * 100.)
    text3 = r'$\overline{{\sigma_{{T_1}}}}={:.4f}\;mag$'.format(
        np.mean(zip(*t1_e)[1]))
    text = text1 + '\n' + text2 + '\n' + text3
    ob = offsetbox.AnchoredText(text, loc=9, prop=dict(size=12))
    ob.patch.set(alpha=0.95)
    ax.add_artist(ob)

    # Output png file.
    fig.tight_layout()
    plt.savefig('t1_photom_errors.png', dpi=300, bbox_inches='tight')


def main():
    '''
    Obtain photometric error for the T1 magnitude, in the range 18.4-19.1
    '''

    # Read data from ASteCA output file.
    t1_e, perc = get_t1_data()

    make_plot(t1_e, perc)


if __name__ == "__main__":
    main()
