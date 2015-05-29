

import numpy as np
from os.path import exists, join, dirname, realpath, isfile, getsize
from os import makedirs, getcwd, listdir
import shutil
from itertools import cycle, chain


def done_move(dst_dir, data_files, folds_num):
    '''
    Move 'data_file' to 'dst_dir' directory.
    '''

    # Cycle through sub-folders reversing when end of list is reached.
    folds_num_rev = reversed(folds_num)
    c = cycle(chain(folds_num, folds_num_rev))

    # For each cluster file.
    for i, cl_f in enumerate(data_files):

        # Name of sub-folder where file is moved.
        sub_dir = dst_dir + c.next()

        # If the sub-dir doesn't exist, create it before moving the file.
        if not exists(sub_dir):
            makedirs(sub_dir)
        try:
            shutil.copy(cl_f, sub_dir)
            print '{} data file copied.'.format(cl_f.split('/')[-1])
        except:
            print "Data file already exists in destination folder."


def get_files_by_file_size(mypath, reverse=False):
    """
    Return list of file paths in directory sorted by file size
    Use reverse=True to return larger files first in list.
    """

    # Get list of files
    filepaths = []
    for basename in listdir(mypath):
        filename = join(mypath, basename)
        if isfile(filename):
            filepaths.append(filename)

    # Re-populate list with filename, size tuples
    for i in xrange(len(filepaths)):
        filepaths[i] = (filepaths[i], getsize(filepaths[i]))

    # Sort list by file size
    # If reverse=True sort from largest to smallest
    # If reverse=False sort from smallest to largest
    filepaths.sort(key=lambda filename: filename[1], reverse=reverse)

    # Re-populate list with just filenames
    for i in xrange(len(filepaths)):
        filepaths[i] = filepaths[i][0]

    return filepaths


def main():

    # Get path where the code is running
    # mypath = realpath(join(getcwd(), dirname(__file__)))
    mypath = '/media/rest/github/asteca-project/asteca/input/dont_read/MC_all/'

    # Read clusters file names and get their sizes. Larger files are located
    # first in the list.
    clusts_names_sz = get_files_by_file_size(mypath, reverse=True)

    # Number of folders to create.
    folds_num = [str(_ + 1).zfill(2) for _ in np.arange(30)]
    dst_dir = mypath + 'input_'

    # Move files into folders.
    done_move(dst_dir, clusts_names_sz, folds_num)

    print 'End.'


if __name__ == "__main__":
    main()
