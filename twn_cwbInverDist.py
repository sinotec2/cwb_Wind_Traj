#!/cluster/miniconda/envs/py37/bin/python
from pandas import *
import twd97, sys
import numpy as np
from scipy.io import FortranFile

dir = '/home/backup/data/cwb/e-service/read_web/'
dfS = read_csv(dir + 'stats_tab.csv')
# drop the closed station
dfss = dfS.loc[dfS.END.map(lambda x: x == '\u3000')]
# drop the precipitation stations
dfS = dfss.loc[dfss.stno.map(lambda x: x[:2] != 'C1')]
no_data=['466850', '467550', '467790']
dfS = dfS.loc[dfS.stno.map(lambda x: x not in no_data)].reset_index(drop=True)
# coordinate transformation of stations
lat, lon = np.array(dfS.LAT), np.array(dfS.LON)
xy = np.array([twd97.fromwgs84(i, j) for i, j in zip(lat, lon)])
x, y = (xy[:, i] for i in [0, 1])
dfS['xy'] = x // 1000 * 10000 + y // 1000
dfS['twd97_x'] = x
dfS['twd97_y'] = y
col='stno,stat_name,twd97_x,twd97_y'.split(',')
dfS[col].set_index('stno').to_csv('stat_wnd.csv')
# sys.exit('OK')
# grid_generation
Longitude_Pole = 120.9900
Latitude_Pole = 23.61000
nx, ny, delta_xy = 83, 137, 3000
x0_LCP, y_LCP = -124500, -205500
xcent, ycent = twd97.fromwgs84(Latitude_Pole, Longitude_Pole)
x0, y0 = np.array([x0_LCP, y_LCP]) + np.array([xcent, ycent])
xe, ye = np.array([x0, y0]) + [nx * delta_xy, ny * delta_xy]
x_mesh = [i for i in range(int(x0 / 1000) - 1, int(xe / 1000) + 2)]
y_mesh = [i for i in range(int(y0 / 1000) - 1, int(ye / 1000) + 2)]
x_g, y_g = np.meshgrid(x_mesh, y_mesh)

# distance square inverse as a weighting
R2 = np.zeros(shape=(len(y_mesh), len(x_mesh), len(dfS)))
for s in range(len(dfS)):
  xm = x_g - x[s]/1000.
  ym = y_g - y[s]/1000.
  R2[:, :, s] = 1. / (xm * xm + ym * ym)
for j in range(len(y_mesh)):
  for i in range(len(x_mesh)):
    sR2 = sum(R2[j, i, :])
    R2[j, i, :] = R2[j, i, :] / sR2

# store the matrix
fnameO = 'R414_252_431.bin'
with FortranFile(fnameO, 'w') as f:
  f.write_record(R2)
with FortranFile('x_mesh.bin', 'w') as f:
  f.write_record(x_mesh)
with FortranFile('y_mesh.bin', 'w') as f:
  f.write_record(y_mesh)
