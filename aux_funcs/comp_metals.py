
from os import walk
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def read_final_data():
    '''
    Read data from 'asteca_output_final.dat'.
    '''

    run, in_data = 'none', []
    with open('../asteca_output_final.dat', "r") as f:
        for line in f:
            if line.strip().startswith('#'):
                if line.strip().startswith('#>>>'):
                    run = line.split()[1]
            else:
                # Store run, name, and z
                in_data.append([run, line.split()[0], line.split()[20]])

    return in_data


def search_cluster(out_path, cl):
    """
    Extract NOT rounded metallicity values from the 'pantalla_analysis.out'
    files.
    """
    cl_found, mets = False, []
    # Store subdir names and file names.
    for root, dirs, files in walk(out_path):
        for f in files:
            # Store paths to output images.
            if f.startswith('pantalla'):
                # Go through 'pantalla_analysis' file.
                with open(root + f) as f1:
                    f_store, n_store = [], []
                    for i, line in enumerate(f1):
                        f_store.append(line)
                        if line.startswith('Analyzing') or\
                                line.startswith('Analizing'):
                                n_store.append(line.split()[2])
                        if line.startswith('Best fit param'):
                            met = f_store[i - 1].split()[2][1:-1]
                            if cl[1] == n_store[0]:
                                print cl[0], cl[1], cl[2], f, n_store[0],\
                                    met
                                cl_found = True
                                mets = [cl[2], met]
                                break
                            else:
                                n_store = []

    return cl_found, mets


def get_met_data(in_data):
    '''
    For each cluster in 'asteca_output_final.dat', store its rounded and
    its NOT rounded metallicity value.

    Clusters from the 18th_ru/2nd_run can't be matched since I did not store
    their 'pantalla_analysis.out' files.
    '''

    mets = [[], []]
    for cl in in_data:
        # Path to data file.
        out_path = '../runs/' + str(cl[0]) + '_run/output/err_out/'
        # Search cluster in 'err_out'
        cl_found, mets0 = search_cluster(out_path, cl)
        if mets0:
            mets[0].append(float(mets0[0]))
            mets[1].append(float(mets0[1]))
        # I cluster was not found, try in second folder.
        if not cl_found:
            # print cl, 'not found'
            out_path = '../runs/' + str(cl[0]) + '_run/output/err_out_2/'
            cl_found, dummy = search_cluster(out_path, cl)
            if not cl_found:
                print cl, 'not found a second time'

    return mets


def make_plot(mets):
    """
    """
    fig = plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(2, 2)

    # Rounded minus real values
    met_diffs = np.array(mets[0]) - np.array(mets[1])
    print 'Mean z diffs:', np.mean(met_diffs)

    ax0 = plt.subplot(gs[0])
    plt.xlabel('z (rounded - real)')
    ax0.hist(met_diffs, 20)
    plt.axvline(x=np.mean(met_diffs), lw=3, ls='--', color='r')
    plt.annotate(r'$\overline{{\Delta z}}={:.6f}$'.format(np.mean(met_diffs)),
                 xy=(0.75, 0.95), xycoords='axes fraction')

    ax1 = plt.subplot(gs[1])
    xmin, xmax = -0.0001, 0.02
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.xlabel('z round')
    plt.ylabel('z real')
    plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')
    ax1.scatter(mets[0], mets[1], s=50.)

    # [Fe/H] values.
    feh_round, feh_real = np.log10(np.array(mets[0]) / 0.0152), \
        np.log10(np.array(mets[1]) / 0.0152)
    # Rounded minus real values
    feh_diffs = feh_round - feh_real
    print 'Mean [Fe/H] diffs:', np.mean(feh_diffs)

    ax2 = plt.subplot(gs[2])
    plt.xlabel('[Fe/H] (rounded - real)')
    ax2.hist(feh_diffs, 20, color='g')
    plt.axvline(x=np.mean(feh_diffs), lw=3, ls='--', color='r')
    plt.annotate(r'$\overline{{\Delta [Fe/H]}}={:.4f}$'.format(
        np.mean(feh_diffs)), xy=(0.05, 0.95), xycoords='axes fraction')

    ax3 = plt.subplot(gs[3])
    xmin, xmax = -2.5, 0.5
    plt.xlim(xmin, xmax)
    plt.ylim(xmin, xmax)
    plt.xlabel('[Fe/H] round')
    plt.ylabel('[Fe/H] real')
    plt.plot([xmin, xmax], [xmin, xmax], 'k', ls='--')
    ax3.scatter(feh_round, feh_real, s=50., c='g')

    # Output png file.
    fig.tight_layout()
    plt.savefig('mets_round_comp.png', dpi=300, bbox_inches='tight')


def main():
    '''
    Compare metallicity values, reals (obtained with the GA) versus
    rounded.
    This is to assess the impact of issue #248.
    '''

    # Read data from ASteCA output file.
    in_data = read_final_data()

    mets = get_met_data(in_data)

    make_plot(mets)


if __name__ == "__main__":
    main()
