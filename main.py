import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
from copy import deepcopy

#data = pd.read_csv('icecube_10year_ps/events/IC86_VII_exp.csv', sep = '\s+')
#data = pd.read_json('/events/IC_86_VII_with_gal.json')
data = pd.read_json('IC_86_VII_galactic.json')
en_data = data['log10(E/GeV)']
N, en = np.histogram(en_data, bins = 500)
en1 = [(en[i+1]+en[i])/2 for i in range(len(en)-1)]

def dNdE(N, en):
    res = []
    for i in range(len(en)-1):
        res.append(N[i]/(en[i+1]-en[i]))
    return np.array(res)
diff_N = dNdE(N, en)

Ra = data['RA[deg]']
Dec = data['Dec[deg]']

# for i in data.iterrows():
#     if i[1]['log10(E/GeV)']> 8:
#         print(i[1]['RA[deg]'], i[1]['Dec[deg]'],)
#         plt.plot(i[1]['RA[deg]'], i[1]['Dec[deg]'], 'b+')


def selector(data, ra_bounds = (0, 360), dec_bounds = (-90,90), deepcopy_req = True):
    if deepcopy_req:
        return deepcopy(data.loc[(ra_bounds[0] < data['GaLon']) & (data['GaLon'] < ra_bounds[1]) &
                                 (dec_bounds[0] < data['GaLat']) & (data['GaLat'] < dec_bounds[1])])
    return data.loc[(ra_bounds[0] < data['GaLon']) & (data['GaLon'] < ra_bounds[1]) &
                    (dec_bounds[0] < data['GaLat']) & (data['GaLat'] < dec_bounds[1])]

data_new = selector(data, (25,100), (-5,5))
# for i in range(len(data)):
#     x = data['RA[deg]'][i]
#     y = data['Dec[deg]'][i]
#     point = SkyCoord(x * u.deg, y * u.deg, frame='icrs').transform_to(Galactic())
#     data.at[i,'GaLon'] = point.l.degree
#     data.at[i,'GaLat'] = point.b.degree
#     # try:
#     #     data.at[i,'dL (Mpc)'] = float(data['dL (Mpc)'][i].split(" ")[0])
#     # except ValueError:
#     #     data = data.drop(i)
#
# data.to_excel("IC_86_VII_galactic.xlsx")
# data.to_json("IC_86_VII_galactic.json")
# data.to_csv("IC_86_VII_galactic.csv")
# np.savetxt("IC_86_VII_galactic.txt", data.values, fmt = '%s', delimiter= '        ', header='      '.join(data.columns))

#plt.show()