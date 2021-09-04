from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import split


def dNdE(N, en):
    res = []
    for i in range(len(en) - 1):
        res.append(N[i] / (en[i + 1] - en[i]))
    return np.array(res)


def rect_window(dots, q):
    vertices = ((dots[0][0], dots[0][0], dots[0][1], dots[0][1], dots[0][0]),
                (dots[1][0], dots[1][1], dots[1][1], dots[1][0], dots[1][0]))
    return np.array([(vertices[0][0], j) for j in np.linspace(vertices[1][0], vertices[1][1], q)] +
                    [(i, vertices[1][1]) for i in np.linspace(vertices[0][1], vertices[0][2], q)] +
                    [(vertices[0][2], j) for j in np.linspace(vertices[1][2], vertices[1][3], q)] +
                    [(i, vertices[1][3]) for i in np.linspace(vertices[0][3], vertices[0][4], q)])


def selector_rect(data, ra_bounds=(0, 360), dec_bounds=(-90, 90), deepcopy_req=True,
                  lon_name='GaLon', lat_name='GaLat'):
    if deepcopy_req:
        return deepcopy(data.loc[(ra_bounds[0] < data[lon_name]) & (data[lon_name] < ra_bounds[1]) &
                                 (dec_bounds[0] < data[lat_name]) & (data[lat_name] < dec_bounds[1])])
    return data.loc[(ra_bounds[0] < data[lon_name]) & (data[lon_name] < ra_bounds[1]) &
                    (dec_bounds[0] < data[lat_name]) & (data[lat_name] < dec_bounds[1])]


def selector(data, bounds=((0, -90), (0, 90), (360, 90), (360, -90)), en_th=1, deepcopy_req=True,
             lon_name='GaLon', lat_name='GaLat'):
    from matplotlib import pyplot as plt
    area = Polygon(bounds)
    x, y = area.exterior.xy
    plt.plot(x, y, 'r-')
    res_rows = []
    for n, i in data.iterrows():
        if Point(i[lon_name], i[lat_name]).within(area):
            res_rows.append(n)
    if deepcopy_req:
        return deepcopy(data.loc[res_rows])
    return data.loc[res_rows]


def Gal2RaDec(points):
    res = []
    for i in points:
        point = SkyCoord(i[0] * u.deg, i[1] * u.deg, frame='galactic').transform_to(ICRS())
        res.append((point.ra.degree, point.dec.degree))
    return np.array(res)


def RaDec2Gal(points, out_of_bounds_fix):
    if not out_of_bounds_fix:
        res = []
        for i in points:
            point = SkyCoord(i[0] * u.deg, i[1] * u.deg, frame='icrs').transform_to(Galactic())
            res.append((point.l.degree, point.b.degree))
        return np.array(res)
    else:
        res1 = []
        res2 = []
        gal_points = []
        flag = 1
        old_point = SkyCoord(points[0][0] * u.deg, points[0][1] * u.deg, frame='icrs').transform_to(Galactic())
        for i in points:
            point = SkyCoord(i[0] * u.deg, i[1] * u.deg, frame='icrs').transform_to(Galactic())
            if abs(old_point.l.degree - point.l.degree) > 300:
                flag = (flag * 2) % 3
            if flag == 1:
                res1.append((point.l.degree, point.b.degree))
            if flag == 2:
                res2.append((point.l.degree, point.b.degree))
            old_point = point
        # gal_points= np.array(gal_points)
        # span = max(gal_points[:,0])- min(gal_points[:,0])
        # for i in gal_points:
        #     if point.l.degree > 360:
        #         res1.append((point.l.degree, point.b.degree))
        #     else:
        #         res2.append((point.l.degree, point.b.degree))
        if len(res2) > 0:
            res2.append(res2[-1])
            return [np.array(res1), np.array(res2)]
        else:
            return [np.array(res1)]


def window_shift(win, q):
    shift0 = max(win[:, 0]) - min(win[:, 0])
    pol = Polygon(win)
    cuts = []
    for i in np.linspace(min(win[:, 1]), max(win[:, 1]), 150):
        line = LineString([(0, i), (360, i)])
        parts = split(line, pol)
        if len(parts) > 2:
            cuts.append(parts[1].length)
    shift = max(max(cuts) + 1, (360 - shift0) / q)
    print(shift, cuts)
    return [np.array([((i + shift * k) % 360, j) for i, j in win]) for k in range(1, q + 1)]


def background_calculator(data, window_bounds, off_zones_number, E=1, rect=True, rect_interpolation_rate=100,
                          en_name='log10(E/GeV)', **kwargs):
    '''
    :param data: data with points to be analized
    :param window_bounds: if rect is True:bounds of rectangular window. Must be array of 2 tuples ((x1,x2), (y1,y2)).
    If rect if False: array of tooples each is verticle of polygon [(x1,y1),(x2,y2),...(xn,yn)]
    :param off_zones_number: number of zones to be used for background calculation
    :param E: energy treshhold. Any event below E will not be considered
    :param rect: if True zone is rectangle, else polygon. Consider proper number of points in polygon, because it will
    undergoo non-linear trasformation
    :param rect_interpolation_rate: numper of dots each side of rectangle will be split on
    :param kwargs: key ford arguments for selector.  Check selector function arguments
    :return: array of integer. First is number of events in window, others - number of events in off zones.
    '''

    data = deepcopy(data.loc[data[en_name] > np.log10(E)])
    window = rect_window(window_bounds, rect_interpolation_rate)
    windowRaDec = Gal2RaDec(window)
    off_zones_RaDec = window_shift(windowRaDec, off_zones_number)
    off_zones = []

    plt.figure(1)
    plt.plot(windowRaDec[:, 0], windowRaDec[:, 1], 'b-', label='window')
    number_of_parts = []
    for i in off_zones_RaDec:
        plt.plot(i[:, 0], i[:, 1], 'r.', label='off-zones')
        gal_zone = RaDec2Gal(i, True)
        number_of_parts.extend([len(gal_zone)] + [0] * (len(gal_zone) - 1))
        off_zones.extend(gal_zone)
    plt.title('Window and off-zones in isrc')
    plt.legend(loc=3)
    plt.show()

    plt.figure(2)
    plt.plot(window[:, 0], window[:, 1], 'b-', label='window')
    sempl = selector(data, window, en_th=E, **kwargs)
    numbers = [sempl.shape[0]]
    plt.plot(sempl['GaLon'], sempl['GaLat'], 'g+')
    plt.text(window[0, 0], window[0, 1], str(numbers[-1]))
    sempls = []
    for i in off_zones:
        sempls.append(selector(data, i, en_th=E, **kwargs))

    for n, i in enumerate(off_zones):
        N = 0
        for j in range(number_of_parts[n]):
            N = N + sempls[n + j].shape[0]
        if N:
            plt.text(i[0, 0], i[0, 1], str(N))
            numbers.append(N)
        plt.plot(sempls[n]['GaLon'], sempls[n]['GaLat'], 'g+')
        plt.plot(i[:, 0], i[:, 1], 'r-', label='off-zones')
    plt.title('Window and off-zones in galactic coordinates')
    plt.legend(loc=4)
    plt.show()
    numbers = np.array(numbers)
    background = np.average(numbers[1:])
    print('number of events in zones ', numbers)
    print('average background is', background)
    print('non background events in window', numbers[0] - background)
    return numbers


def bounds_solver(zone):
    # finction is obsolet and should not be used
    pass
    # line = LineString([(360,-5), (360,5)])
    # new_areas = []
    # out_of_bounds_index = []
    # area = Polygon(zone)
    # areas = split(area, line)
    # new_areas = [np.array(list(zip(i.exterior.xy[0],i.exterior.xy[1]))) for i in areas]
    # for i, pol in zip(range(len(new_areas)),  new_areas):
    #     if min(zone[:,0])> 360:
    #         new_areas[i][:,0] = new_areas[i][:,0]%360
    # return new_areas
