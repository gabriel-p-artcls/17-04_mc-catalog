
from astropy.coordinates import Angle, SkyCoord
from astropy import units as u
import numpy as np
from deproj_dist import deproj_dist as dd
from functions.MCs_data import MCs_data
import scipy.interpolate


def ccc(l1, l2):
    '''
    Concordance correlation coefficient.
    See: https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
    '''
    return 2 * np.cov(l1, l2)[0, 1] / (np.var(l1) + np.var(l2) +
                                       (np.mean(l1) - np.mean(l2)) ** 2)


def xy_interp(rang, N):
    return np.linspace(rang[0], rang[1], N)


def get_deproj_dist(ra, dec, cent, dist, inc, pa):
    '''
    Obtain deprojected distance between cluster and the center of the MC
    in kpc.
    '''
    dist_kpc = []
    for r, d in zip(*[ra, dec]):
        c = SkyCoord(ra=r*u.degree, dec=d*u.degree)
        d_kpc = dd(c, cent, pa, inc, dist)[1].value
        dist_kpc.append(d_kpc)

    return dist_kpc


def gsd(in_params):
    '''
    Calculate the best match for the inclination and position angles of the
    MCs, based on the distance assigned to each cluster by ASteCA.
    '''
    ra, dec, dist_cent, aarr = [in_params[_] for _ in [
        'ra', 'dec', 'dist_cent', 'aarr']]
    dist_cent_kpc = [np.asarray(dist_cent[0]) / 1000.,
                     np.asarray(dist_cent[1]) / 1000.]

    # # S/LMC central coords stored in degrees.
    # c_SMC = SkyCoord('00h52m45s', '-72d49m43s', frame='icrs')
    # # ^ (13.1875, -72.82861111)
    # c_LMC = SkyCoord('05h20m57s', '-69d28m41s', frame='icrs')
    # # ^ (80.2375, -69.47805556)
    # cent = [c_SMC, c_LMC]

    # # S/LMC distance stored in parsecs.
    # # ~ 61944.11 pc (18.96 mag)
    # d_SMC = Distance(10 ** (0.2 * (18.96 + 5)) / 1000., unit=u.kpc)
    # # ~ 49888.45 pc (18.49 mag)
    # d_LMC = Distance(10 ** (0.2 * (18.49 + 5)) / 1000., unit=u.kpc)
    # dist = [d_SMC, d_LMC]

    # Inclination range, S/LMC.
    inc_rang = [[20., 65.], [20., 65.]]
    # Position angle range S/LMC.
    pa_rang = [[100., 300.], [100., 300.]]

    # Grid size: N x N
    N = 35

    gal_str_pars, best_inc_pa_ang = [], []
    # SMC, LMC = 0, 1
    for j in [0, 1]:
        # Center and distance to galaxy.
        cent, dist, e_dm_dist = MCs_data(j)
        # Inclination and position angles grid.
        x_slmc, y_slmc = xy_interp(inc_rang[j], N), xy_interp(pa_rang[j], N)
        # Inclination and position angles grid in degrees.
        inc_lst = [Angle(_, unit=u.degree) for _ in x_slmc]
        pa_lst = [Angle(_, unit=u.degree) for _ in y_slmc]

        # Obtain deprojected distances for each cluster.
        x, y, z = [], [], []
        for inc in inc_lst:
            for pa in pa_lst:
                dist_kpc = get_deproj_dist(ra[j], dec[j], cent, dist, inc, pa)
                c = ccc(dist_kpc, dist_cent_kpc[j])
                x.append(inc.degree)
                y.append(pa.degree)
                z.append(c)
                print inc.degree, pa.degree, c
        print '\n'

        # Use data from file to speed things up.
        # x, y, z = [], [], []
        # gal_name = ['SMC', 'LMC']
        # with open('functions/temp_' + gal_name[j] + '.dat') as f:
        # # with open('temp_' + gal_name[j] + '.dat') as f:
        #     for line in f:
        #         l = line.split()
        #         x.append(Angle(l[0], unit=u.degree).degree)
        #         y.append(Angle(l[1], unit=u.degree).degree)
        #         z.append(float(l[2]))

        # As arrays.
        x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
        # Limits.
        xmin, xmax = inc_rang[j][0] - 0.1, inc_rang[j][1] + 0.1
        ymin, ymax = pa_rang[j][0] - 0.1, pa_rang[j][1] + 0.1
        # Interpolating function.
        z = z.reshape(N, N)
        rbs = scipy.interpolate.RectBivariateSpline(x_slmc, y_slmc, z)
        # Finer grid to interpolate on.
        xi, yi = xy_interp(inc_rang[j], 200), xy_interp(pa_rang[j], 200)
        # Get values on grid.
        zi = rbs(xi, yi)

        # Obtain maximum density value.
        max_idx = np.unravel_index(zi.argmax(), zi.shape)
        best_inc_pa_ang.append([xi[max_idx[0]], yi[max_idx[1]]])
        # print zi.max()

        # 10% max values.
        ten_perc = int(len(zi.flatten()) * 0.1)
        # Get indices.
        indices = np.argpartition(zi.flatten(), -ten_perc)[-ten_perc:]
        # Get coordinates.
        max_coords = np.vstack(np.unravel_index(indices, zi.shape)).T
        # Store all coordinates.
        x_lst, y_lst = [], []
        for c in max_coords:
            x_lst.append(xi[c[0]])
            y_lst.append(yi[c[1]])
        # calculate standard deviation.
        x_std, y_std = np.std(np.asarray(x_lst)), np.std(np.asarray(y_lst))

        gal_str_pars.append([xmin, xmax, ymin, ymax, xi, yi, zi, max_idx,
                             x_std, y_std])

    # Append deprojected distances obtained using the max density angles.
    for j in [0, 1]:
        dist_kpc = get_deproj_dist(ra[j], dec[j], cent, dist,
                                   Angle(best_inc_pa_ang[j][0], unit=u.degree),
                                   Angle(best_inc_pa_ang[j][1], unit=u.degree))
        # import matplotlib.pyplot as plt
        # plt.scatter(dist_cent_kpc[j], dist_kpc)
        # plt.show()
        c = ccc(dist_kpc, dist_cent_kpc[j])
        gal_str_pars.append([-0.01, 8.05, -0.01, 8.05, dist_cent_kpc[j],
                            dist_kpc, aarr[j][0], c, '', ''])

    return gal_str_pars


# if __name__ == "__main__":
#     ra = [[15.5958333333333, 12.1375, 10.9333333333333, 17.225, 15.1416666666667, 15.0041666666667, 14.3416666666667, 13.5625, 357.245833333333, 15.1041666666667, 11.3583333333333, 16.0916666666667, 14.4583333333333, 16.8333333333333, 22.6583333333333, 9.425, 18.875, 18.2583333333333, 13.3541666666667, 24.0041666666667, 11.6833333333333, 23.75, 23.3083333333333, 17.5541666666667, 25.4291666666667, 4.60416666666667, 19.5666666666667, 25.5916666666667, 11.8, 15.2833333333333, 21.2333333333333, 11.475, 29.1833333333333, 18.0166666666667, 5.66666666666667, 14.45, 22.8833333333333, 14.4458333333333, 13.275, 25.6166666666667, 11.2166666666667, 13.1458333333333, 10.7458333333333, 13.8875, 15.9708333333333, 14.425, 6.17916666666667, 20.7, 18.2125, 27.3666666666667, 5.3625, 10.35, 23.6083333333333, 5.76666666666667, 17.0791666666667, 15.1458333333333, 27.5791666666667, 11.9583333333333, 16.0166666666667, 12.3625, 17.6958333333333, 15.2333333333333, 15.4916666666667, 11.5041666666667, 14.325, 15.2041666666667, 17.0583333333333, 14.0583333333333, 16.8791666666667, 16.7375, 13.6291666666667, 12.5875, 14.3333333333333, 15.35, 17.2583333333333, 12.325, 15.1375, 10.8875, 12.1541666666667, 11.775, 16.2583333333333, 11.7291666666667, 10.9083333333333, 12.05, 11.8541666666667, 12.0041666666667, 15.7958333333333, 17.2541666666667, 15.0583333333333], [76.4166666666667, 74.225, 82.9916666666667, 78.8541666666667, 76.1416666666667, 88.0458333333333, 72.3083333333333, 76.5083333333333, 78.8125, 87.7, 76.9458333333333, 77.3, 76.1375, 77.2208333333333, 73.9208333333333, 86.4583333333333, 76.9833333333333, 83.3333333333333, 72.7958333333333, 77.7333333333333, 86.0458333333333, 72.2791666666667, 72.5875, 77.625, 86.7166666666667, 71.8583333333333, 72.6208333333333, 78.9458333333333, 76.975, 74.725, 86.7125, 73.7625, 72.25, 94.3291666666667, 76.5375, 91.8708333333333, 74.5583333333333, 93.4833333333333, 72.1541666666667, 70.8083333333333, 73.4625, 79.7583333333333, 76.55, 82.9416666666667, 74.5416666666667, 76.9416666666667, 78.1041666666667, 74.3916666666667, 76.6416666666667, 85.9833333333333, 81.55, 74.9416666666667, 69.925, 79.5208333333333, 82.6416666666667, 74.9083333333333, 82.2083333333333, 95.3916666666667, 79.4541666666667, 67.6666666666667, 77.7958333333333, 85.375, 77.3958333333333, 76.4708333333333, 69.4125, 76.5125, 74.8083333333333, 76.6041666666667, 85.8958333333333, 82.4833333333333, 79.1666666666667, 77.7875, 78.45, 76.125, 82.925, 73.1875, 82.85, 72.7, 92.2208333333333, 93.9875, 88.8958333333333, 85.8333333333333, 74.1208333333333, 84.7833333333333, 80.05, 72.4125, 73.55, 71.5166666666667, 77.3458333333333, 77.9208333333333, 77.65, 81.3625, 74.7125, 81.1166666666667, 79.2333333333333, 81.125, 68.9083333333333, 93.6166666666667, 91.6291666666667, 74.3583333333333, 73.7541666666667, 83.0125, 86.55, 93.6708333333333, 74.9708333333333, 76.8958333333333, 79.2208333333333, 79.0708333333333, 77.9166666666667, 82.4416666666667, 77.3125, 83.2541666666667, 77.5083333333333, 83.6625, 74.5625, 80.4375, 79.6708333333333, 85.4916666666667, 71.6041666666667, 73.225, 76.9041666666667, 76.4, 86.4833333333333, 85.1083333333333, 77.6666666666667, 74.5916666666667, 77.7208333333333, 73.9666666666667, 82.3333333333333, 85.4541666666667, 77.6333333333333, 79.1125, 80.0083333333333, 71.875, 88.925, 85.6208333333333, 84.4416666666667, 86.3625, 80.8, 83.0958333333333, 77.2125, 82.9625, 76.8375, 76.6708333333333, 83.5541666666667, 82.4958333333333, 85.4125, 82.4125, 76.3541666666667, 76.6416666666667]] 
#     dec = [[-72.0030555555556, -73.3069444444444, -72.9766666666667, -73.2416666666667, -72.3655555555556, -72.3688888888889, -71.8913888888889, -72.2413888888889, -72.9452777777778, -71.2947222222222, -73.4813888888889, -72.8477777777778, -72.9436111111111, -73.3775, -76.0544444444444, -73.9083333333333, -71.1788888888889, -70.9627777777778, -72.1963888888889, -75.4577777777778, -72.0630555555556, -75.5566666666667, -74.1672222222222, -73.2091666666667, -71.1611111111111, -74.3186111111111, -72.0016666666667, -74.1733333333333, -73.4772222222222, -74.0736111111111, -71.1836111111111, -73.5066666666667, -74.2194444444445, -75.1975, -75.0747222222222, -73.4216666666667, -71.9527777777778, -74.3266666666667, -73.3802777777778, -71.2791666666667, -73.0019444444445, -72.1930555555556, -72.5886111111111, -74.0636111111111, -72.8261111111111, -74.4727777777778, -73.755, -75.0016666666667, -73.1194444444444, -73.7283333333333, -73.7486111111111, -72.8908333333333, -72.8744444444444, -73.6697222222222, -72.8841666666667, -71.4608333333333, -74.3561111111111, -73.4783333333333, -74.6191666666667, -73.3980555555556, -72.7919444444444, -73.1516666666667, -71.0202777777778, -73.3955555555556, -72.9327777777778, -73.3488888888889, -73.2569444444444, -74.1561111111111, -73.1197222222222, -73.235, -74.1852777777778, -73.3872222222222, -71.1702777777778, -73.2402777777778, -73.0863888888889, -73.3716666666667, -72.2583333333333, -73.4388888888889, -73.4155555555556, -73.3730555555556, -73.0427777777778, -73.4233333333333, -73.4405555555556, -73.4463888888889, -73.4580555555556, -73.4861111111111, -72.2736111111111, -73.2066666666667, -72.4583333333333], [-68.6394444444445, -68.0022222222222, -67.9716666666667, -68.6811111111111, -68.2083333333333, -71.8583333333333, -72.0566666666667, -68.0263888888889, -68.8825, -71.7077777777778, -66.7980555555556, -68.4441666666667, -67.9755555555556, -68.0836111111111, -67.7833333333333, -69.3802777777778, -67.3577777777778, -68.1522222222222, -67.5347222222222, -67.6266666666667, -69.3333333333333, -67.3416666666667, -72.8275, -68.4005555555556, -69.1897222222222, -67.6597222222222, -67.3258333333333, -69.1919444444445, -67.9288888888889, -67.8469444444445, -69.4197222222222, -67.9644444444444, -72.64, -70.0608333333333, -68.4458333333333, -72.4941666666667, -67.7683333333333, -72.5052777777778, -68.5594444444444, -73.8119444444444, -69.5719444444444, -69.0011111111111, -68.0630555555556, -68.2355555555556, -68.0602777777778, -67.8613888888889, -68.7719444444444, -65.2677777777778, -68.3630555555556, -69.1805555555556, -70.9813888888889, -69.8011111111111, -74.0172222222222, -69.1716666666667, -63.2033333333333, -69.5575, -71.6327777777778, -72.79, -68.4727777777778, -66.9569444444444, -67.6266666666667, -69.185, -67.8105555555556, -67.0494444444445, -66.1994444444444, -68.6283333333333, -67.9083333333333, -68.375, -66.2086111111111, -72.0547222222222, -70.5408333333333, -67.6825, -66.62, -68.3497222222222, -72.1461111111111, -72.5177777777778, -72.0425, -72.5766666666667, -72.3838888888889, -70.0730555555556, -62.3452777777778, -66.2622222222222, -67.6227777777778, -74.8533333333333, -68.9041666666667, -72.2480555555556, -69.8069444444444, -66.9113888888889, -67.7783333333333, -67.5655555555556, -70.4875, -73.5702777777778, -69.9577777777778, -67.7286111111111, -68.6827777777778, -67.6780555555556, -73.7316666666667, -72.6094444444445, -72.2263888888889, -67.6852777777778, -67.7141666666667, -64.2422222222222, -69.0825, -69.8019444444445, -67.9236111111111, -67.4608333333333, -69.15, -69.1541666666667, -68.7266666666667, -71.0005555555556, -67.6997222222222, -67.8491666666667, -66.6991666666667, -68.3055555555556, -68.0491666666667, -68.9172222222222, -69.0794444444444, -69.0475, -72.5683333333333, -72.1725, -68.5419444444444, -68.6286111111111, -69.2719444444444, -69.2486111111111, -68.7536111111111, -69.8030555555556, -67.4711111111111, -69.7058333333333, -70.5794444444445, -68.9208333333333, -66.94, -69.0802777777778, -69.2611111111111, -72.5883333333333, -74.3538888888889, -65.3627777777778, -74.7827777777778, -69.3452777777778, -70.7777777777778, -67.9969444444444, -67.9802777777778, -67.9911111111111, -66.8291666666667, -67.8422222222222, -67.8563888888889, -67.8788888888889, -69.2294444444444, -70.9838888888889, -68.5005555555556, -68.4272222222222]]
#     dist_cent = [[1633.20060682, 2320.18925372, 2441.9621715, 1471.75739309, 2376.23276698, 769.377260812, 1217.10160803, 1804.52411396, 5941.22951677, 2131.04017704, 909.015121115, 3069.45508982, 2364.25345255, 2670.26753775, 4890.13489057, 2343.49016847, 3159.23688725, 2711.83412564, 1818.55754823, 4455.35179015, 1123.62228704, 4757.29044595, 3641.00953753, 2648.4875134, 4474.93618064, 3141.57483709, 3560.24707724, 4473.53808798, 1395.16822521, 1590.48289743, 4018.53487003, 1077.40931903, 5110.49822201, 3786.85646257, 4265.8368639, 949.628191291, 3943.07061992, 1767.67615703, 1299.74797416, 4641.567255, 2418.05685486, 897.32557739, 3039.33043949, 1461.11210457, 1901.98375256, 2165.3084445, 2985.15151754, 3722.76959311, 2000.6292342, 4885.31051319, 2625.48765461, 1076.03465615, 3317.63420596, 2482.45477321, 1670.70686978, 1704.67888631, 4679.63812859, 3032.23144391, 2116.49256805, 2422.17347074, 1542.38398274, 3014.72270399, 2369.65890352, 1867.37485958, 688.161984621, 2390.12368317, 1417.8418462, 2759.56312304, 1647.94174103, 2532.76451148, 2666.69527657, 1799.77351614, 1911.04658283, 2378.98229824, 1319.33831563, 646.650576733, 1952.83158151, 2437.45811475, 1828.61711731, 1838.06758207, 1147.07949676, 2459.7173664, 1123.54113051, 1845.94527547, 2373.7888037, 2374.78095982, 1177.87078532, 1347.286255, 724.047475304], [2094.92442356, 2592.17806506, 1926.01374032, 1393.13441938, 2375.35464572, 3301.57725704, 3191.52180512, 1846.0135317, 716.76297379, 3154.98748878, 3076.59112567, 2026.05665608, 2404.83900491, 1555.09594563, 2707.74040478, 2869.88446607, 2126.74995801, 1875.64558761, 3013.12651671, 1938.91000725, 2787.32652108, 3248.71116208, 3895.8901948, 2000.81445911, 2338.12168352, 3264.69637398, 3074.99893799, 2636.47366843, 2026.34617439, 2555.05591583, 1987.45137574, 2962.78454192, 3645.08442741, 4482.58313968, 1839.78891743, 4308.71233216, 2341.82811636, 4655.60143549, 2645.31421226, 4649.97712828, 2659.06169124, 495.958510033, 1712.84745196, 2951.12877171, 2751.29010002, 1869.62116249, 1815.09401991, 4178.76760298, 1500.54520148, 1794.75124939, 2529.48578966, 1640.54629482, 5595.45850131, 776.552314833, 5577.42226175, 2249.44973927, 3280.94442372, 5192.25922245, 2306.05901115, 4799.58913132, 1930.71867668, 1618.96956185, 2738.38874315, 2453.83795634, 4563.65453258, 1390.93463993, 2210.73606753, 2156.22099955, 3489.11690674, 2339.47934036, 1913.24925842, 1765.47926278, 2664.78715134, 1747.9825195, 2559.15521594, 3985.98208869, 2470.42710796, 3449.61772162, 4247.88799427, 4176.76382235, 6981.56554962, 3438.49595369, 2544.39520164, 4846.09325026, 859.207154012, 4248.51583614, 3334.30573069, 3728.40654495, 1759.97165254, 2178.23559988, 1194.19833024, 3858.08899466, 1864.27857663, 1561.90783843, 2706.82540924, 1747.18252909, 4836.99690823, 4641.18832985, 4140.00763706, 2964.26958036, 3078.63612913, 4672.26090061, 2310.20822426, 4180.41248267, 2162.15029365, 2072.33397643, 1238.51032921, 832.357478845, 1488.41649353, 2593.59405836, 1824.3311299, 2309.75855988, 2853.29775447, 2151.58363423, 2434.44194297, 2082.03137409, 1198.78596379, 1997.28834526, 3626.80120077, 3532.99815015, 2924.8945852, 2101.1318558, 2466.86972041, 2506.56843199, 2352.67147924, 2754.38598176, 1932.30496879, 2459.48921061, 1639.6977746, 1695.42812081, 2380.69363485, 2640.55720245, 2121.15585755, 3588.5380412, 4855.28448834, 4857.60200313, 4765.37072491, 1887.71370931, 1942.00303845, 1926.04452804, 2570.73288389, 1574.74408319, 3067.96627803, 1827.77703955, 1778.39064842, 2280.49016348, 1949.46253097, 2136.73909356, 1495.06656829, 1870.72119252]]
#     aarr = [[[9.0, 8.9, 8.75, 9.9, 8.2, 8.1, 9.1, 7.9, 9.6, 9.6, 8.9, 9.6, 9.1, 9.2, 9.2, 9.45, 9.1, 9.3, 7.4, 9.6, 9.85, 9.1, 9.7, 9.55, 9.35, 9.15, 9.4, 9.0, 8.45, 9.7, 8.95, 8.0, 7.95, 9.45, 9.45, 9.65, 9.0, 9.05, 8.8, 9.3, 9.05, 8.0, 8.9, 9.5, 7.8, 9.8, 9.15, 9.8, 9.6, 9.55, 9.6, 9.45, 9.7, 9.75, 8.95, 9.6, 8.2, 8.8, 9.2, 8.7, 8.3, 8.4, 9.1, 8.5, 8.8, 8.35, 8.75, 9.6, 8.8, 8.3, 9.4, 8.3, 8.4, 8.0, 8.0, 8.5, 7.8, 8.7, 8.0, 8.4, 8.0, 8.1, 8.0, 8.0, 8.2, 6.9, 7.0, 7.2, 6.7], [9.20411998265592, 9.14612803567824, 8.79934054945358, 9.82607480270083, 8.04139268515822, 8.09691001300806, 8.79934054945358, 7.90308998699194, 9.77815125038364, 9.73239375982297, 9.17609125905568, 9.36172783601759, 9.30102999566398, 9.39794000867204, 9.20411998265592, 9.32221929473392, 9.10037054511756, 9.44715803134222, 7.39794000867204, 9.81291335664286, 9.77815125038364, 9.04139268515822, 9.47712125471966, 9.73239375982297, 9.07918124604762, 9.06069784035361, 9.32221929473392, 9.04139268515822, 8.39794000867204, 9.96848294855393, 8.77815125038364, 8.20411998265593, 8.04139268515822, 9.57978359661681, 9.61278385671973, 9.49136169383427, 9.14612803567824, 9.30102999566398, 8.5051499783199, 9.25527250510331, 9.07918124604762, 8.14612803567824, 9.0, 9.68124123737559, 7.39794000867204, 9.73239375982297, 9.30102999566398, 9.67209785793572, 9.63346845557959, 9.77815125038364, 9.49136169383427, 9.32221929473392, 9.80617997398389, 9.51851393987789, 9.07918124604762, 9.77815125038364, 8.14612803567824, 8.69897000433602, 9.44715803134222, 8.35, 8.25, 8.4, 9.15, 8.4, 8.65, 8.1, 8.65, 9.45, 8.45, 8.1, 9.25, 8.4, 7.9, 8.1, 8.3, 8.05, 7.7, 8.35, 7.9, 8.1, 8.0, 7.9, 7.8, 8.1, 8.4, 8.34242268082221, 7.9, 8.15, 8.1]], [[8.45, 8.2, 8.2, 8.7, 8.9, 9.3, 8.9, 9.15, 8.4, 8.6, 8.8, 9.1, 8.95, 8.8, 8.3, 8.6, 9.0, 9.5, 8.8, 8.6, 8.7, 9.1, 9.05, 9.1, 9.2, 8.6, 9.0, 8.7, 8.3, 9.05, 9.1, 9.1, 9.0, 9.0, 8.2, 9.2, 8.9, 9.05, 9.05, 9.3, 9.1, 9.15, 8.65, 9.4, 9.2, 9.05, 8.8, 9.5, 7.9, 9.0, 8.4, 8.95, 9.2, 8.3, 9.4, 8.3, 8.9, 9.05, 9.5, 9.2, 9.05, 8.0, 8.2, 9.1, 9.15, 8.3, 8.9, 8.65, 9.1, 8.85, 9.1, 8.7, 8.8, 9.0, 9.3, 9.25, 8.75, 9.05, 9.3, 9.4, 9.3, 9.45, 8.95, 9.3, 8.0, 9.3, 9.1, 8.65, 9.1, 8.6, 8.95, 9.4, 9.15, 9.5, 8.75, 8.8, 9.55, 9.15, 9.2, 9.1, 9.15, 9.2, 9.45, 9.2, 9.15, 9.0, 8.8, 8.0, 8.0, 9.5, 7.9, 9.0, 8.3, 8.8, 8.4, 8.8, 8.85, 8.1, 9.55, 8.3, 9.15, 8.8, 8.2, 7.9, 8.8, 8.3, 8.9, 8.1, 8.3, 8.0, 8.4, 8.7, 8.6, 9.25, 9.2, 9.45, 9.25, 8.8, 9.1, 8.8, 9.55, 7.9, 8.6, 9.1, 7.0, 6.7, 7.7, 7.3, 8.85, 7.5], [8.0, 8.1, 8.25, 8.69897000433602, 9.14612803567824, 9.30102999566398, 9.14612803567824, 9.11394335230684, 8.20411998265593, 8.5051499783199, 8.95424250943932, 9.1, 8.75, 8.60205999132796, 8.35, 8.45, 9.04139268515822, 9.39794000867204, 8.60205999132796, 8.69897000433602, 8.7, 9.09691001300805, 9.11394335230684, 8.89762709129044, 9.30102999566398, 8.55, 8.9, 8.75, 8.25, 9.23044892137827, 9.0, 8.89762709129044, 9.04139268515822, 9.14612803567824, 8.0, 9.25527250510331, 8.69897000433602, 9.0, 8.84509804001426, 9.20411998265592, 9.0, 9.11394335230684, 8.55, 9.23044892137827, 9.14612803567824, 9.17609125905568, 8.69897000433602, 9.34242268082221, 7.9, 9.09691001300805, 8.4, 9.20411998265592, 9.39794000867204, 8.3, 9.41497334797082, 8.05, 8.95424250943932, 8.79934054945358, 9.20411998265592, 9.20411998265592, 9.04139268515822, 8.25, 8.25, 9.14612803567824, 9.17609125905568, 8.09691001300806, 8.69897000433602, 8.5051499783199, 9.17609125905568, 8.84509804001426, 9.17609125905568, 8.60205999132796, 8.79934054945358, 8.7481880270062, 9.27875360095283, 9.20411998265592, 8.60205999132796, 8.95424250943932, 9.23044892137827, 9.14612803567824, 9.38021124171161, 9.36172783601759, 8.85125834871907, 9.25527250510331, 8.2, 9.20411998265592, 9.11394335230684, 8.6, 9.14612803567824, 8.60205999132796, 9.07918124604762, 9.25527250510331, 9.17609125905568, 9.34242268082221, 8.65321251377534, 8.69897000433602, 9.39794000867204, 9.04139268515822, 9.25527250510331, 9.20411998265592, 9.20411998265592, 9.23044892137827, 9.36172783601759, 9.23044892137827, 9.17609125905568, 9.14612803567824, 8.3, 8.3, 7.69897000433602, 9.25527250510331, 8.14612803567824, 9.0, 8.14612803567824, 8.60205999132796, 8.45, 8.55, 8.9, 8.2, 9.30102999566398, 8.39794000867204, 9.07918124604762, 8.60205999132796, 8.1, 8.25, 8.39794000867204, 8.45, 8.65321251377534, 7.95, 8.11394335230684, 8.39794000867204, 8.0, 8.20411998265593, 8.1, 9.11394335230684, 9.04139268515822, 9.47712125471966, 9.20411998265592, 8.60205999132796, 9.20411998265592, 8.79934054945358, 9.25527250510331, 8.15, 8.34242268082221, 9.0, 8.15, 8.3, 8.25, 7.9, 7.69897000433602, 8.35]]]
#     in_params = {'ra': ra, 'dec': dec, 'dist_cent': dist_cent, 'aarr': aarr}

#     gsd(in_params)
