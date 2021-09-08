import pandas as pd
from functions import *



data = pd.read_json('IC_86_VII_galactic.json')
#background_calculator(data, ((25,100), (-5,5)),9, E = 2e3)

bounds = (1e3, 4e4)
data_x = data.loc[(data['log10(E/GeV)'] > 3) & (data['log10(E/GeV)'] < 4.4)]
en_data = 10**(selector(data_x, rect_window(((25,100), (-5,5)), 10))['log10(E/GeV)'])

#N, en = np.histogram(en_data, bins = 500)
#en1 = [(en[i+1]+en[i])/2 for i in range(len(en)-1)]
N = []
en = []
N_clear = []
for i in np.linspace(bounds[0], bounds[1], 20):
    N.append(en_data.loc[en_data >i].shape[0])
    N_clear.append(background_calculator(data, ((25,100), (-5,5)),9, E = i))
    en.append(i)
plt.figure(3)
plt.plot(en, N, 'g-')
plt.plot(en, N_clear, 'b-')
# for i, E in enumerate(en1):
#     if N[i]> 10:
#         plt.plot(E, N[i], 'k.')
# diff_N = dNdE(N, en)
plt.show()

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


