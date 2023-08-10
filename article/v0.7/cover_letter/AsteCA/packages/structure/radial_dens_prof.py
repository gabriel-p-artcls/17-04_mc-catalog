
import math


def main(hist_lst, cent_bin):
    """
    Calculate the density profile by counting the number of stars in the center
    bin first (r aprox width_bin/2 px), then moving to the 8 adjacent bins
    (r aprox width_bin + (width_bin/2) px), then the following 16 (r aprox
    2*width_bin + (width_bin/2) px), then the next 24 (r aprox 3*width_bin +
    (width_bin/2) px) and so forth. These "square rings" have consecutive
    multiples of 8 bins each (if no bin falls outside the range of the x,y
    chart) and an area of N*width_bin^2 [px^2], where N is the number of bins.
    Dividing the number of stars in each "square ring" by this area, we get all
    the "ring densities" for those approximate values of r.
    """

    hist_2d, bin_width = hist_lst[0], hist_lst[-1]
    x_c_b, y_c_b = cent_bin

    # Initialize lists.
    radii, rdp_points, poisson_error = [], [], []

    # Percentage of the frame where the RDP will be calculated.
    rdp_perc = 0.95
    # Total length in both axis.
    # x_length, y_length = len(hist_2d[0]), len(hist_2d[1])
    # Minimum length in either axis.
    # min_length = min(x_length, y_length)
    # Length where the RDP will be defined, in bin units.
    # rdp_length = min_length * rdp_perc

    # Number of bins that define the length of the largest square ring
    # around the center coordinates.
    # bins = int(rdp_length)
    # Number of "square rings" to generate.
    # sq_rings = int(bins * 0.5)

    # Use max x,y length for the RDP.
    rdp_length = int(max(len(hist_2d), len(hist_2d[0])) * rdp_perc)

    square_rings = []
    # Iterate through all the square rings.
    # for i in range(sq_rings):
    for i in range(rdp_length):
        # Store here the coordinates of the bins.
        bins_coords = []

        # Initialize bin_count for this square ring.
        ring_count, bin_count = 0, 0

        # Iterate through bins in the x dimension for the 2D hist.
        for xindex, xitem in enumerate(hist_2d):
            # Iterate through bins in the y dimension for the 2D hist.
            for yindex, st_in_bin in enumerate(xitem):

                # Bin is in the top row
                if yindex == (y_c_b + i) and abs(xindex - x_c_b) <= i:
                    # Add stars in bin to corresponding ring.
                    ring_count = ring_count + st_in_bin
                    # Add 1 more bin to the "square ring".
                    bin_count += 1
                    bins_coords.append([xindex - x_c_b, i])
                # Bin is in the bottom row
                elif yindex == (y_c_b - i) and abs(xindex - x_c_b) <= i:
                    ring_count = ring_count + st_in_bin
                    bin_count += 1
                    bins_coords.append([xindex - x_c_b, -i])
                # Bin is in the left column
                elif xindex == (x_c_b - i) and abs(yindex - y_c_b) <= (i - 1):
                    ring_count = ring_count + st_in_bin
                    bin_count += 1
                    bins_coords.append([-i, yindex - y_c_b])
                # Bin is in the right column
                elif xindex == (x_c_b + i) and abs(yindex - y_c_b) <= (i - 1):
                    ring_count = ring_count + st_in_bin
                    bin_count += 1
                    bins_coords.append([i, yindex - y_c_b])

        if bin_count == 0:
            break

        print('Ring {}, bins total: {} bins: {}, diff: {}'.format(
            i, i*8, bin_count, i*8 - bin_count))
        # Store bin coordinates in each square ring.
        square_rings.append(bins_coords)
        # If no stars are inside this square ring, set value to 1 to avoid a
        # division by zero.
        bin_count = 1 if bin_count == 0 else bin_count
        # The number of bins times the area of each bin gives the area of
        # this square ring.
        area = bin_count * (bin_width ** 2)

        # Calculate density corresponding to "square ring" i
        rdp_points.append(ring_count / area)
        # Obtain the Poisson error bar for each value
        poisson_error.append(math.sqrt(ring_count) / area)

        # Store values for radii to go with the densities obtained above
        # and stored in 'rdp_points'
        radii.append(bin_width / 2. + (bin_width * i))

    # Transform from bin units to coordinate units before passing.
    rdp_length = rdp_length * bin_width
    print('rdp length: {}'.format(rdp_length))

    print 'Radial density profile (RDP) calculated.'

    return radii, rdp_points, poisson_error, square_rings, rdp_length
