#!/opt/anaconda3/envs/py27/bin/python
import numpy as np
from pandas import *
import os, sys, subprocess, time, json
from scipy.io import FortranFile
from datetime import datetime, timedelta
import twd97


def getarg():
  """ read time period and station name from argument(std input)
  traj2kml.py -t daliao -d 20171231 """
  import argparse
  ap = argparse.ArgumentParser()
  ap.add_argument("-t", "--STNAM", required=True, type=str, help="station name,sep by ,or Lat,Lon")
  ap.add_argument("-d", "--DATE", required=True, type=str, help="yyyymmddhh")
  ap.add_argument("-b", "--BACK", required=True, type=str, help="True or False")
  args = vars(ap.parse_args())
  return [args['STNAM'], args['DATE'],args['BACK']]

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def nstnam():
  import json
  fn = open(path+'sta_list.json')
  d_nstnam = json.load(fn)
  d_namnst = {v: k for k, v in d_nstnam.items()}
  return (d_nstnam, d_namnst)


def ws_uv(ws, wd):
  PAI = np.pi
  RAD = (270. - wd) * PAI / 180.0
  u = ws * np.cos(RAD)
  v = ws * np.sin(RAD)
  return u, v


def beyond(xpp, ypp):
  xp_km, yp_km = int(xpp // 1000), int(ypp // 1000)
  boo = not ((xp_km - x_mesh[0]) * (xp_km - x_mesh[-1]) < 0 and (yp_km - y_mesh[0]) * (yp_km - y_mesh[-1]) < 0)
  return [boo, (xp_km, yp_km)]


def opendf(pdate):
  ymd = pdate.strftime('%Y%m%d')
  fname= path+'../' + ymd[:4] + '/cwb' + ymd + '.csv'
  if not os.path.isfile(fname):
    dfT=DataFrame({})
    sys.exit('no file for '+fname)
  else:
    try:
      dfT = read_csv(fname)
      dfT['stno'] = [i[:6] for i in dfT.stno_name]
      dfT = dfT.loc[dfT.stno.map(lambda x: x in stno)].reset_index(drop=True)
      dfT = dfT.fillna(0)
      ws, wd = np.array(dfT.WS), np.array(dfT.WD)
      uv = np.array([ws_uv(i, j) for i, j in zip(ws, wd)])
      dfT['u'], dfT['v'] = (uv[:, i] for i in [0, 1])
      dfT.ObsTime = [int(i)-1 for i in dfT.ObsTime]
    except:
      dfT=DataFrame({})
  return dfT, ymd

def uvb(r,u,v):
  dfuv=DataFrame({'R':r,'u':u,'v':v})
  dfuv=dfuv.sort_values('R',ascending=False).reset_index(drop=True)
  rr,uu,vv=np.array(dfuv.R)[:ns3],np.array(dfuv.u)[:ns3],np.array(dfuv.v)[:ns3]
  rr=rr/sum(rr)
  ub,vb=sum(rr*uu),sum(rr*vv)
  return ub,vb

path='/Users/Data/cwb/e-service/surf_trj/'
# restore the matrix
nx, ny, ns = 252, 414, 431
fnameO = path+'R%d_%d_%d.bin' % (ny, nx, ns)
with FortranFile(fnameO, 'r') as f:
  R2 = f.read_record(dtype=np.float64)
R2 = R2.reshape(ny, nx, ns)
with FortranFile(path+'x_mesh.bin', 'r') as f:
  x_mesh = list(f.read_record(dtype=np.int64))
with FortranFile(path+'y_mesh.bin', 'r') as f:
  y_mesh = list(f.read_record(dtype=np.int64))

(d_nstnam, d_namnst) = nstnam()
stnam, DATE, BACK = getarg()
BACK=str2bool(BACK)
BF=-1
if not BACK:BF=1
bdate = datetime(int(DATE[:4]), int(DATE[4:6]), int(DATE[6:8]), int(DATE[8:]))
nam = [i for i in stnam.split(',')]
if len(nam) > 1:
  try:
    lat = float(nam[0])
    lon = float(nam[1])
  except:
    sys.exit('more than two station, suggest executing iteratively')
  else:
    # in case of lat,lon
    if lat < 90.:
      xy0 = twd97.fromwgs84(lat,lon)
      x0, y0 =([xy0[i]] for i in [0,1])
      nam[0] = str(round(lat,2))+'_'+str(round(lon,2))+'_'
    #   in case of twd97_x,y
    else:
      # test the coordinate unit
      if lat>1000.:
        x0, y0 = [lat],[lon]
        nam[0] = str(int(lat/1000))+'+'+str(int(lon/1000))+'_'
      else:
        x0, y0 = [lat*1000],[lon*1000]
        nam[0] = str(int(lat))+'_'+str(int(lon))+'_'

# len(nam)==1, read the location from csv files
else:
  for stnam in nam:
    if stnam not in d_namnst: sys.exit("station name not right: " + stnam)
  nst = [int(d_namnst[i]) for i in nam]
  # locations of air quality stations
  # read from the EPA web.sprx
  fname = path+'sta_ll.csv'
  sta_list = read_csv(fname)
  x0, y0 = [], []
  for s in nst:
    sta1 = sta_list.loc[sta_list.ID == s].reset_index(drop=True)
    x0.append(list(sta1['twd_x'])[0])
    y0.append(list(sta1['twd_y'])[0])

xp, yp = x0, y0
dfS = read_csv(path+'stat_wnd.csv')
if len(dfS) != ns: sys.exit('ns not right')
stno = list(dfS.stno)
pdate = bdate
df, ymd0 = opendf(pdate)
if len(df)==0:sys.exit('no cwb data for date of:'+ymd0)
delt = 15
s = 0
o_ymdh,o_time,o_xp,o_yp,l_xp,l_yp=[],[],[],[],[],[]
itime=0
ymdh=int(DATE)
o_ymdh.append('ymd='+DATE)
o_time.append('hour='+str(itime))
o_xp.append(xp[s])
o_yp.append(yp[s])
l_xp.append(xp[s])
l_yp.append(yp[s])
ns3=int(ns)
while not beyond(xp[s], yp[s])[0] and len(df)!=0:
  df1 = df.loc[(df.ObsTime == ymdh) & (df.stno.map(lambda x:x in stno))].reset_index(drop=True)
  df1 = df1.drop_duplicates()
  ldf1=len(df1)
  if ldf1 < ns:
    next_date= bdate + timedelta(hours=24*BF)
    boo=pdate>next_date	
    if not BACK:boo=pdate<next_date
    if boo:
      ns2 = set(df1.stno)
      miss = set(stno) - set(ns2)
      if len(miss)!=0:
        for m in miss:
          df2= DataFrame({'stno_name':[m],'ObsTime':[ymdh]})
          df1=df1.append(df2,ignore_index=True, sort=False)
    else:
      print 'df1 not right' + str(ymdh)
      break
  df1=df1.sort_values('stno_name').reset_index(drop=True)
  df1 = df1.fillna(0)
  u, v = np.array(df1.u), np.array(list(df1.v))
  for sec in range(0, 3601, delt):
    boo, (xp_km, yp_km) = beyond(xp[s], yp[s])
    if boo: break
    ix, iy = x_mesh.index(xp_km), y_mesh.index(yp_km)
    if sec == 0:
      ix0, iy0 = ix, iy
      ub, vb = uvb(R2[iy, ix, :],u,v)
    else:
      if ix0 != ix or iy0 != iy:
        # ub, vb = sum(R2[iy, ix, :] * u), sum(R2[iy, ix, :] * v)
        ub, vb = uvb(R2[iy, ix, :],u,v)
        ix0, iy0 = ix, iy
    xp[s], yp[s] = xp[s]+BF*delt * ub, yp[s]+BF*delt * vb
    l_xp.append(xp[s])	
    l_yp.append(yp[s])	
  pdate = pdate + timedelta(hours=BF)
  ymdh = int(pdate.strftime('%Y%m%d%H'))
  itime+=1
  o_ymdh.append('ymd='+str(ymdh))
  o_time.append('hour='+str(itime))
  o_xp.append(xp[s])
  o_yp.append(yp[s])
  if pdate.strftime('%Y%m%d') != ymd0:
    df, ymd0 = opendf(pdate)
df=DataFrame({'ymdh':o_ymdh,'xp':o_xp,'yp':o_yp,'Hour':o_time})
col=['xp','yp','Hour','ymdh']
name='trj_results/'+'trj'+nam[0]+DATE+'.csv'
df[col].set_index('xp').to_csv(name)

# output the line segments for each delta_t
dfL=DataFrame({'TWD97_x':l_xp,'TWD97_y':l_yp})
dfL.set_index('TWD97_x').to_csv(name.replace('.csv','L.csv'))

#make kml file
dir='NC'
if not BACK:dir='RC'
os.system('csv2kml.py -f '+name+' -n '+dir+' -g TWD97')
os.system('csv2bln.cs '+name)
