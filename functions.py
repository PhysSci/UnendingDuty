from copy import deepcopy

def selector(data, ra_bounds = (0, 360), dec_bounds = (-90,90), deepcopy_req = True):
    if deepcopy_req:
        return deepcopy(data.loc[(ra_bounds[0] < data['GaLon']) & (data['GaLon'] < ra_bounds[1]) &
                                 (dec_bounds[0] < data['GaLat']) & (data['GaLat'] < dec_bounds[1])])
    return data.loc[(ra_bounds[0] < data['GaLon']) & (data['GaLon'] < ra_bounds[1]) &
                    (dec_bounds[0] < data['GaLat']) & (data['GaLat'] < dec_bounds[1])]