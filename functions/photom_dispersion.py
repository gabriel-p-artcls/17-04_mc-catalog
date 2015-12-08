
import numpy as np
import CMD_obs_vs_asteca as cmd


def get_disp(cl_reg_fit):
    '''
    Calculate the 2D photometric dispersion for the CMD.
    '''
    p_mean_col = np.mean(cl_reg_fit[0])
    p_mean_mag = np.mean(cl_reg_fit[1])
    col_disp = (np.array(cl_reg_fit[0]) - p_mean_col) ** 2
    mag_disp = (np.array(cl_reg_fit[1]) - p_mean_mag) ** 2
    # print p_mean_col, p_mean_mag
    p_disp = sum(np.sqrt(col_disp + mag_disp)) / len(cl_reg_fit[0])
    return p_disp


def main(r_path, cl):
    '''
    '''
    # Fetch which run holds this cluster's membership data.
    run = cmd.get_cl_run(cl)
    # Fetch what 'input_XX' folder in the above run contains the
    # membership file.
    inpt = cmd.get_input_folder(r_path, cl, run)
    # Membership data for cluster.
    cl_reg_fit, cl_reg_no_fit = cmd.get_memb_data(r_path, run, inpt, cl)

    # Obtain photometric dispersion
    p_disp = get_disp(cl_reg_fit)
    # print cl, p_disp
    return p_disp


if __name__ == "__main__":
    main()
