
from os.path import join, realpath, dirname
from os import getcwd
from os import walk
import shutil


def skip_comments(f):
    '''
    Read lines that DO NOT start with a # symbol.
    '''
    for line in f:
        if not line.strip().startswith('#'):
            yield line


def get_clusts_names(part):
    '''
    Read asteca_output_part_xx.dat and store the names of the files.
    '''

    # Path to data file.
    out_file = 'asteca_used_' + part + '.dat'

    # Read data file
    with open(out_file) as f:
        as_names = []

        for line in skip_comments(f):
            as_names.append(line.split()[0])

    return as_names


def get_files_names(clusts_names, mypath):
    """
    Return list of file paths in directory sorted by file size
    Use reverse=True to return larger files first in list.
    """

    # Get list of files
    filepaths = []

    # Store subdir names and file names in cl_files.
    for root, dirs, files in walk(mypath):
        for f in files:
            # Store paths to output images.
            if f.endswith('.png') and not f.endswith('_tiers.png'):
                if f[:-4] in clusts_names:
                    filepaths.append(join(root, f))

    return filepaths


def done_move(dst_dir, files_names):
    '''
    Go through all sub-folders and search for each cluster in the partial
    asteca_output_part.dat file. Copy each match to the given output folder.
    '''

    # For each cluster file.
    for i, cl_f in enumerate(files_names):
        try:
            shutil.copy(cl_f, dst_dir)
            print '{}, {} data file copied.'.format(i + 1, cl_f.split('/')[-1])
        except:
            print "Data file already exists in destination folder."


def main():

    # Get path where the code is running
    mypath = realpath(join(getcwd(), dirname(__file__)))
    print mypath

    # Partial output file.
    num = '3rd'
    part = mypath[:-9]
    print part

    # Read clusters file names.
    clusts_names = get_clusts_names(num)
    print clusts_names

    # Read files names.
    path = part + '/runs/' + num + '_run/output/'
    print path
    files_names = get_files_names(clusts_names, path)
    print files_names

    # Number of folders to create.
    dst_dir = '/home/gabriel/Descargas/mc_asteca_out/'

    # Copy files into destination folder.
    done_move(dst_dir, files_names)

    print 'End.'


if __name__ == "__main__":
    main()
