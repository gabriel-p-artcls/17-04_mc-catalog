
import numpy as np
import random
from astroML.plotting import hist
import operator
from ..inp import input_params as g
from bayesian_da import sort_members


def bin_edges_f(bin_method, mag_col_cl):
    '''
    Obtain bin edges for each photometric dimension using the cluster region
    diagram.
    '''
    bin_edges = []
    if bin_method in ['sturges', 'sqrt']:
        if bin_method == 'sturges':
            b_num = 1 + np.log2(len(mag_col_cl[0]))
        else:
            b_num = np.sqrt(len(mag_col_cl[0]))

        for mag_col in mag_col_cl:
            bin_edges.append(np.histogram(mag_col, bins=int(b_num))[1])

    elif bin_method == 'bb':
        # Based on Bonatto & Bica (2007) 377, 3, 1301-1323. Fixed bin width
        # of 0.25 for colors and 0.5 for magnitudes.
        b_num = [(max(mag_col_cl[0]) - min(mag_col_cl[0])) / 0.25,
                 (max(mag_col_cl[1]) - min(mag_col_cl[1])) / 0.5]

        for i, mag_col in enumerate(mag_col_cl):
            bin_edges.append(np.histogram(mag_col, bins=int(b_num[i]))[1])

    else:
        for mag_col in mag_col_cl:
            bin_edges.append(hist(mag_col, bins=bin_method)[1])

    return bin_edges


def get_clust_histo(memb_prob_avrg_sort, mag_col_cl, bin_edges):
    '''
    Generate the N-dimensional cluster region histogram, with each star
    positioned in its corresponding cell.
    '''

    # Cluster region N-dimensional histogram.
    cl_hist = np.histogramdd(np.array(zip(*mag_col_cl)), bins=bin_edges)[0]
    # np.shape(cl_hist) gives the tuple containing one element per dimension,
    # indicating how many cells that dimension was divided into.

    # Add a very small amount to each outer-most edge so the 'np.digitize'
    # function will position the stars on the edges correctly.
    for i, b_e in enumerate(bin_edges):
        bin_edges[i][0] = b_e[0] - (abs(b_e[0]) / 100.)
        bin_edges[i][-1] = b_e[-1] + (b_e[-1] / 100.)

    # Position each cluster region star in its corresponding N-dimensional
    # cell/bin.
    cl_st_indx = []
    # Store indexes for each dimension.
    for i, mag_col in enumerate(mag_col_cl):
        # Set correct indexes for array substracting 1, since 'np.digitize'
        # counts one more bin to the right by default.
        cl_st_indx.append(np.digitize(mag_col, bin_edges[i]) - 1)

    # Create empty list with the same dimension as the photometric N-histogram.
    cl_hist_p = np.empty(shape=cl_hist.shape + (0,)).tolist()
    # Position stars in their corresponding N-histogram cells. Since the stars
    # are already sorted by their MPs, they will be correctly sorted in the
    # final list here too.
    for i, h_indx in enumerate(zip(*cl_st_indx)):
        # Store stars.
        reduce(operator.getitem, list(h_indx), cl_hist_p).append(
            memb_prob_avrg_sort[i])

    return cl_hist_p, cl_hist


def get_fl_reg_hist(field_region, bin_edges, cl_hist):
    '''
    Obtain the average number of field region stars in each cell defined for
    the N-dimensional cluster region photometric diagram.
    '''

    # Empty field region array shaped like the cluster region array.
    f_hist = np.zeros(shape=np.shape(cl_hist))
    # Add stars in all the defined field regions.
    for freg in field_region:
        Q = np.array(zip(*freg)[1:])
        # Create list with all magnitudes and colors defined.
        mag_col_fl = [Q[4], Q[2]]
        f_hist = f_hist + np.histogramdd(np.array(zip(*mag_col_fl)),
                                         bins=bin_edges)[0]

    # Average number of stars in each cell/bin and round to integer.
    f_hist = np.around(f_hist / len(field_region), 0)

    return f_hist


def get_fit_stars(cl_hist_p, f_hist, flag_decont_skip):
    '''
    Iterate through each N-dimensional cell of the cluster region array and
    remove the excess of field stars in each one, selecting those with the
    lowest assigned MPs if the DA was applied. Otherwise select random stars.
    '''

    # Only flatten list if more than 1 cell was defined.
    if len(cl_hist_p) > 1:
        cl_hist_p_flat = np.asarray(cl_hist_p).flatten()
    else:
        cl_hist_p_flat = cl_hist_p[0]

    # Flatten arrays to access all of its elements.
    f_hist_flat = f_hist.flatten()

    red_memb_fit, red_memb_no_fit = [], []
    # For each cell defined.
    for i, cl_cell in enumerate(cl_hist_p_flat):

        # Get average number of field regions in this cell.
        N_fl_reg = f_hist_flat[i]

        if N_fl_reg > 0.:
            # Discard the excess of N_reg_fl stars from this cluster region.

            # If the DA was not applied, discard N_fl_reg *random* stars in
            # the cell.
            if flag_decont_skip:

                if int(N_fl_reg) < len(cl_cell):
                    # Generate list with randomized cell indexes.
                    ran_indx = random.sample(
                        xrange(len(cl_cell)), len(cl_cell))

                    # Store len(cl_cell) - N_fl_reg stars
                    red_memb_fit.append([cl_cell[i] for i in
                                         ran_indx[:-int(N_fl_reg)]])
                    # Discard N_fl_reg stars.
                    red_memb_no_fit.append([cl_cell[i] for i in
                                            ran_indx[-int(N_fl_reg):]])
                else:
                    # Discard *all* stars in the cell.
                    red_memb_no_fit.append(cl_cell)
            else:
                # Discard those N_fl_reg with the smallest MPs, keep the rest.
                red_memb_fit.append(cl_cell[:-int(N_fl_reg)])
                red_memb_no_fit.append(cl_cell[-int(N_fl_reg):])
        else:
            # No field region stars in this cell, keep all stars.
            red_memb_fit.append(cl_cell)

    # Flatten lists of stars and re-sort according to highest MPs.
    red_memb_fit = sort_members([i for sublst in red_memb_fit for i in sublst])
    red_memb_no_fit = sort_members([i for sublst in red_memb_no_fit for i in
                                    sublst])

    # Minimum probability of  selected stars.
    min_prob = red_memb_fit[-1][-1]

    return red_memb_fit, red_memb_no_fit, min_prob


def main(decont_algor_return, field_region):
    '''
    Takes the photometric diagram (CMD, CCD, etc.) of the cluster region with
    assigned MPs, divides it into sub-regions (cells) according to the
    density within it, and removes in each sub-region a number of stars
    equal to the average excess due to field star contamination.
    '''

    memb_prob_avrg_sort, flag_decont_skip = decont_algor_return
    local_bin = g.rm_params[1]

    # Remove ID's and zip.
    P = np.array(zip(*memb_prob_avrg_sort)[1:], dtype='float')
    # Create list with all magnitudes and colors defined.
    mag_col_cl = [P[4], P[2]]

    # Obtain bin edges.
    bin_edges = bin_edges_f(local_bin, mag_col_cl)

    # Obtain N-dimensional cluster region histogram.
    cl_hist_p, cl_hist = get_clust_histo(memb_prob_avrg_sort, mag_col_cl,
                                         bin_edges)

    # Obtain field regions histogram (only number of stars in each cell).
    f_hist = get_fl_reg_hist(field_region, bin_edges, cl_hist)

    # Obtain stars separated in list to be used by the BF func and list of
    # those discarded stars.
    red_memb_fit, red_memb_no_fit, min_prob = get_fit_stars(
        cl_hist_p, f_hist, flag_decont_skip)

    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # plt.scatter(P[4], P[2])
    # ax.set_xticks(bin_edges[0], minor=False)
    # ax.set_yticks(bin_edges[1], minor=False)
    # ax.xaxis.grid(True, which='major')
    # ax.yaxis.grid(True, which='major')
    # ax.invert_yaxis()
    # plt.show()

    # Store and pass for plotting purposes.
    red_plot_pars = [min_prob, bin_edges]

    # Check the number of stars selected.
    if len(red_memb_fit) < 10:
        print ("  WARNING: less than 10 stars left after reducing\n"
               "  by 'local' method. Using full list.")
        red_memb_fit, red_memb_no_fit, red_plot_pars = memb_prob_avrg_sort, \
            [], [0.]

    return red_memb_fit, red_memb_no_fit, red_plot_pars
