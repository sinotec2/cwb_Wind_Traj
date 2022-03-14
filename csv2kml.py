#!/opt/anaconda3/envs/py37/bin/python
from pandas import *
import twd97,sys
from pyproj import Proj
import numpy as np


def getarg():
    """ read the setting of plot from argument(std input)"""
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fname", required = True, type=str,help = "csv with x,y,Name,Label")
    ap.add_argument("-n", "--NorH", required = True, type=str,help = "Normal/Highlight/Dot/Reddot/Blackdot/Polygon")
    ap.add_argument("-g", "--GEOG", required = True, type=str,help = "LL or TWD97")
    args = vars(ap.parse_args())
    return args['fname'],args['NorH'],args['GEOG']

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False


fname,NorH,geog=getarg()
NorH,geog=NorH.upper(),geog.upper()
LL,TWD97=False,False
if geog=='LL':
  LL=True
elif geog=='TWD97':
  TWD97=True
else:
  sys.exit('GEOG must be LL or TWD97, should not be:'+geog)
if NorH[0] not in 'HNDRBP' :sys.exit('NorH not right')
NH=NorH.replace('C','').replace('L','')

#B:target with black arona
#R:red target
#D:dot target
stl={'N':'normalPlacemark','H':"highlightPlacemark",'D':"dotPlacemark",'B':"blackdotPlacemark",'R':"reddotPlacemark"}
#csv with 4 columns,
#xp,yp,Hour,ymdh
#266567.0,2549275.0,hour=0,ymd=2019062315
df=read_csv(fname)#,encoding='big5')
if isfloat(df.columns[0]):
  df=read_csv(fname,header=None)#,encoding='big5')
  STR='XYND'
  df.columns=[STR[i] for i in range(len(df.columns))]
if len(df.columns)==2:
  df['Name']=['pt'+str(i) for i in df.index]
  df['Desc']=['Desc'+str(i) for i in df.index]
TITLE=fname
head0='<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document>'+ \
'<name>'+TITLE+'</name><description>'+TITLE+'</description>'
#store the lines
line=[head0]

#icons may be used in KML file:
sty_pngs={'H':'paddle/red-stars.png','N':'paddle/wht-blank.png','B':'pal4/icon25.png','D':'pal4/icon57.png','R':'pal4/icon49.png',}
if 'P' not in NorH:
  #declarations of Placemark styles
  for s in stl:
    stl_url='<Style id="'+stl[s]+'"><IconStyle><Icon><href>http://maps.google.com/mapfiles/kml/'+\
	sty_pngs[s]+'</href> </Icon> </IconStyle></Style>'
    exec('styl'+s+"='"+stl_url+"'")
    exec('line.append(styl'+s+')')
else:
  #colors prepare 10 colors, only one will be adopted, which is determined by random generator
  N=10
  col = '#00FF0A #3FFF0A #7FFF0A #BFFF0A #FFFF0A #FECC0A #FD990A #FC660A #FB330A #FA000A'.replace('#', '').split()
  if len(col) != N: print ('color scale not right, please redo from http://www.zonums.com/online/color_ramp/')
  aa = '28'  # ''28'~ 40%, '4d' about 75%
  rr, gg, bb = ([i[j:j + 2] for i in col] for j in [0, 2, 4])
  col = [aa + b + g + r for b, g, r in zip(bb, gg, rr)]
  col0 = '4d6ecdcf'
  col_line0 = 'cc2d3939'
  st_head = ''
  st_med = '</color><width>1</width></LineStyle><PolyStyle><color>'
  st_tail = '</color></PolyStyle></Style>'
  for i in range(N):
    st_head += '<Style id="level' + str(i) + '"><LineStyle><color>' + col[i] + st_med + col[i] + st_tail
  line.append(st_head)

col=df.columns
#plot the placemarks
if TWD97:
  Latitude_Pole, Longitude_Pole = 23.61000, 120.990
  Xcent, Ycent = twd97.fromwgs84(Latitude_Pole, Longitude_Pole)  
  pnyc = Proj(proj='lcc', datum='NAD83', lat_1=10, lat_2=40, 
        lat_0=Latitude_Pole, lon_0=Longitude_Pole, x_0=0, y_0=0.0)
  x_lpc,y_lpc=np.array(df[col[0]])-Xcent,np.array(df[col[1]])-Ycent
# ll = np.array([twd97.towgs84(i, j) for i, j in zip(df.loc[:,col[0]],df.loc[:,col[1]])])
  long,lati=pnyc(x_lpc, y_lpc, inverse=True)
  lonlat=[str(lon)+','+str(lat)+',0' for lon,lat in zip(long,lati)]
elif LL:
  lonlat=[str(lon)+','+str(lat)+',0' for lon,lat in zip(df.loc[:,col[0]],df.loc[:,col[1]])]

if 'P' not in NorH:
  for i in range(len(df)):
    nam=df.loc[i,col[2]]
    desc=df.loc[i,col[3]]
    sturl= '<styleUrl>#'+stl[NH]+'</styleUrl>'
    line.append('<Placemark><name>'+nam+'</name><description>'+desc+'</description>'+sturl+'<Point><coordinates>')
    line.append(lonlat[i]+'</coordinates></Point></Placemark>')
  #if with CorL in NorH, connection the marks with line
  if 'C' in NorH:
    line.append('<Placemark> <LineString> <coordinates>')
    for i in range(len(df)):
      line.append(lonlat[i])
    line.append('</coordinates></LineString></Placemark>')
  if 'L' in NorH:
    df=read_csv(fname.replace('.csv','L.csv'))#,encoding='big5')
    if len(df)>1000:
      ii=int(len(df)/1000)
      if ii>1:
        df=df.loc[[i for i in range(0,len(df),ii)]].reset_index(drop=True)
    col=df.columns
    line.append('<Placemark> <LineString> <coordinates>')
    if TWD97:
#     ll = np.array([twd97.towgs84(i, j) for i, j in zip(df.loc[:,col[0]],df.loc[:,col[1]])])
#     lonlat=[str(lon)+','+str(lat)+',0' for lon,lat in zip(ll[:,1],ll[:,0])]
      x_lpc,y_lpc=np.array(df[col[0]])-Xcent,np.array(df[col[1]])-Ycent
      long,lati=pnyc(x_lpc, y_lpc, inverse=True)
      lonlat=[str(lon)+','+str(lat)+',0' for lon,lat in zip(long,lati)]
    if LL:
      lonlat=[str(lon)+','+str(lat)+',0' for lon,lat in zip(df.loc[:,col[0]],df.loc[:,col[1]])]
    for i in range(len(df)):
      line.append(lonlat[i])
    line.append('</coordinates></LineString></Placemark>')
else: #Polygon case
  headP = '</styleUrl><Polygon><outerBoundaryIs><LinearRing><tessellate>1</tessellate><coordinates>'
  tailP = '</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>'
  ps=[i[-2:] for i in list(df[col[3]])]
  npoly,i0,i1=1,[0],[len(df)]
  if 'p0' in ps:
    npoly=ps.count('p0')
    if npoly>1:
      i1=[]
      for i in range(npoly-1):
        i0.append(i0[-1]+1+ps[i0[-1]+1:].index('p0'))
        i1.append(i0[-1])
      i1.append(len(df))
  for ip in range(npoly):
#   level=int(np.random.rand()*10)
    level=5
    nam=df.loc[i0[ip],col[2]].replace('p0','')
    desc=df.loc[i0[ip],col[3]].replace('p0','')
    line.append('<Placemark><name>' + nam + '</name><description>'+desc+'</description><styleUrl>#level' + str(level) + headP)
    print (i0[ip],i1[ip])
    for i in range(i0[ip],i1[ip]):
      line.append(lonlat[i])
    line.append(tailP)
line.append('</Document></kml>')
with open(fname+'.kml','w') as f:
  [f.write(l+'\n') for l in line]


