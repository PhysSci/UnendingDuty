import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
from functions import *



data = pd.read_json('IC_86_VII_galactic.json')
background_calculator(data, ((25,100), (-5,5)),9, E = 2e3)

# en_data = data['log10(E/GeV)']
# N, en = np.histogram(en_data, bins = 500)
# en1 = [(en[i+1]+en[i])/2 for i in range(len(en)-1)]
#
# diff_N = dNdE(N, en)

# window = rect_window(((25,100), (-5,5)), 100)
# windowRaDec = Gal2RaDec(window)
# off_zones_RaDec = window_shift(windowRaDec, 8)
# off_zones = []
#
# plt.figure(1)
# plt.plot(windowRaDec[:,0], windowRaDec[:, 1], 'b-', label = 'window')
# number_of_parts = []
# for i in off_zones_RaDec:
#     plt.plot(i[:,0], i[:, 1], 'r.', label = 'off-zones')
#     gal_zone = RaDec2Gal(i, True)
#     number_of_parts.extend([len(gal_zone)] + [0]*(len(gal_zone)-1))
#     off_zones.extend(gal_zone)
# plt.title('Window and off-zones in isrc')
# plt.legend(loc = 3)
# plt.show()
#
# plt.figure(2)
# E = 1e3
# plt.plot(window[:,0], window[:, 1], 'b-', label = 'window')
# sempl = selector(data, window, en_th=E)
# numbers= [sempl.shape[0]]
# plt.plot(sempl['GaLon'], sempl['GaLat'], 'g+')
# plt.text(window[0,0], window[0,1], str(numbers[-1]))
# #off_zones = np.array(off_zones)
# #for i in off_zones: off_zones2.extend(bounds_solver(i))
# sempls = []
# for i in off_zones:
#     sempls.append(selector(data, i, en_th=E))
#
# for n,i in enumerate(off_zones):
#     N = 0
#     for j in range(number_of_parts[n]):
#         N = N + sempls[n+j].shape[0]
#     if N:  plt.text(i[0,0], i[0,1], str(N))
#     numbers.append(N)
#     plt.plot(sempls[n]['GaLon'], sempls[n]['GaLat'], 'g+')
#     plt.plot(i[:,0], i[:, 1], 'r-', label = 'off-zones')
# plt.title('Window and off-zones in galactic coordinates')
# plt.legend(loc = 4)
# plt.show()
# print(numbers)
#plt.plot(en[1:], N)
#plt.show()

#dataN = selector(data, (25,100), (-5,5))

# for i in dataN.iterrows():
#     # if i[1]['log10(E/GeV)']> 8:
#     #     print(i[1]['RA[deg]'], i[1]['Dec[deg]'],)
#     plt.plot(i[1]['RA[deg]'], i[1]['Dec[deg]'], 'b+')
# plt.show()


