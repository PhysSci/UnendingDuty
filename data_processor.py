
import pandas as pd
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic



def to_galactic(in_path, out_path):
    data = pd.read_csv(in_path)
    for i in range(len(data)):
        x = data['RA[deg]'][i]
        y = data['Dec[deg]'][i]
        point = SkyCoord(x * u.deg, y * u.deg, frame='icrs').transform_to(Galactic())
        data.at[i,'GaLon'] = point.l.degree
        data.at[i,'GaLat'] = point.b.degree
        # try:
        #     data.at[i,'dL (Mpc)'] = float(data['dL (Mpc)'][i].split(" ")[0])
        # except ValueError:
        #     data = data.drop(i)

    data.to_excel(out_path + ".xlsx")
    data.to_json(out_path + ".json")
    data.to_csv(out_path + ".csv")
    np.savetxt(out_path + ".txt", data.values, fmt = '%s', delimiter= '        ', header='      '.join(data.columns))

#plt.show()