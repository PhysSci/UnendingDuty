import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic
from functions import *



data = pd.read_json('IC_86_VII_galactic.json')
en_data = data['log10(E/GeV)']
N, en = np.histogram(en_data, bins = 500)
en1 = [(en[i+1]+en[i])/2 for i in range(len(en)-1)]

diff_N = dNdE(N, en)


dataN = selector(data, (25,100), (-5,5))

for i in dataN.iterrows():
    # if i[1]['log10(E/GeV)']> 8:
    #     print(i[1]['RA[deg]'], i[1]['Dec[deg]'],)
    plt.plot(i[1]['RA[deg]'], i[1]['Dec[deg]'], 'b+')
plt.show()

