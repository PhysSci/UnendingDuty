from copy import deepcopy
import numpy as np
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import split
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic

def dNdE(N, en):
    res = []
    for i in range(len(en)-1):
        res.append(N[i]/(en[i+1]-en[i]))
    return np.array(res)


def rect_window(dots, q):
    vertices = ((dots[0][0], dots[0][0], dots[0][1], dots[0][1], dots[0][0]),
                (dots[1][0], dots[1][1], dots[1][1], dots[1][0], dots[1][0]))
    return np.array([(vertices[0][0], j) for j in np.linspace(vertices[1][0], vertices[1][1], q)] +
                    [(i, vertices[1][1]) for i in np.linspace(vertices[0][1], vertices[0][2], q)] +
                    [(vertices[0][2], j) for j in np.linspace(vertices[1][2], vertices[1][3], q)] +
                    [(i, vertices[1][3]) for i in np.linspace(vertices[0][3], vertices[0][4], q)])

def selector_rect(data, ra_bounds = (0, 360), dec_bounds = (-90,90), deepcopy_req = True,
                  lon_name = 'GaLon', lat_name = 'GaLat'):
    if deepcopy_req:
        return deepcopy(data.loc[(ra_bounds[0] < data[lon_name]) & (data[lon_name] < ra_bounds[1]) &
                                 (dec_bounds[0] < data[lat_name]) & (data[lat_name] < dec_bounds[1])])
    return data.loc[(ra_bounds[0] < data[lon_name]) & (data[lon_name] < ra_bounds[1]) &
                    (dec_bounds[0] < data[lat_name]) & (data[lat_name] < dec_bounds[1])]


def selector(data, bounds = ((0,-90),(0, 90),(360, 90),(360, -90)),en_th = 1, deepcopy_req = True,
                  lon_name = 'GaLon', lat_name = 'GaLat', en_name = 'log10(E/GeV)'):
    from matplotlib import pyplot as plt
    area = Polygon(bounds)
    x,y = area.exterior.xy
    plt.plot(x,y, 'r-')
    res_rows = []
    for n,i in data.iterrows():
        if Point(i[lon_name], i[lat_name]).within(area) and i[en_name] > np.log10(en_th):
            res_rows.append(n)
    if deepcopy_req:
        return deepcopy(data.loc[res_rows])
    return  data.loc[res_rows]


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
            if abs(old_point.l.degree-point.l.degree) > 300:
                flag = (flag*2)%3
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
        if len(res2) >0:
            res2.append(res2[-1])
            return [np.array(res1), np.array(res2)]
        else:
            return [np.array(res1)]


def window_shift(win,q):
    shift = max(win[:, 0]) - min(win[:,0])
    return [np.array([((i+shift*k)%360,j) for i,j in win]) for k in range(1, q+1)]


def bounds_solver(zone):
        line = LineString([(360,-5), (360,5)])
        new_areas = []
        out_of_bounds_index = []
        area = Polygon(zone)
        areas = split(area, line)
        new_areas = [np.array(list(zip(i.exterior.xy[0],i.exterior.xy[1]))) for i in areas]
        for i, pol in zip(range(len(new_areas)),  new_areas):
            if min(zone[:,0])> 360:
                new_areas[i][:,0] = new_areas[i][:,0]%360
        return new_areas
