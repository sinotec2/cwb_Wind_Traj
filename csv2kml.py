#!/cluster/miniconda/envs/py37/bin/python
from pandas import *
import twd97,json,sys
def getarg():
    """ read the setting of plot from argument(std input)"""
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fname", required = True, type=str,help = "csv xy data")
    ap.add_argument("-n", "--NorH", required = True, type=str,help = "Normal or Highlight or line")
    ap.add_argument("-g", "--GEOG", required = True, type=str,help = "LL or TWD9")
    args = vars(ap.parse_args())
    return args['fname'],args['NorH'],args['GEOG']

fname,NorH,geog=getarg()
NorH=NorH.upper()
LL,TWD97=False,False
if geog=='LL':
  LL=True
if geog=='TWD97':
  TWD97=True
if 'N' not in NorH and 'H' not in NorH :sys.exit('NorH not right')
NH=NorH.replace('L','')
stl={'N':'normalPlacemark','H':"highlightPlacemark",'L':"normalPlacemark"}
a=read_csv(fname)#,encoding='big5')
TITLE=fname
head0='<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Document>'+ \
'<name>'+TITLE+'</name><description>'+TITLE+'</description>'
line=[head0]
stylH='<Style id="highlightPlacemark"><IconStyle><Icon><href>http://maps.google.com/mapfiles/kml/paddle/red-stars.png</href> </Icon> </IconStyle></Style>'
stylN='<Style id="normalPlacemark"><IconStyle><Icon><href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href> </Icon> </IconStyle></Style>'
line.append(stylH)
line.append(stylN)
col=a.columns
lonlat=[]
for i in range(len(a)):
  nam=a.loc[i,col[2]]
  desc=a.loc[i,col[3]]
  line.append('<Placemark><name>'+nam+'</name><description>'+desc+'</description><styleUrl>#'+stl[NH]+'</styleUrl><Point><coordinates>')
  if TWD97:
    x,y=a.loc[i,col[0]],a.loc[i,col[1]] 
    lat,lon=twd97.towgs84(x,y)
  if LL:
    lon,lat=a.loc[i,col[0]],a.loc[i,col[1]] 
  slonlat=str(lon)+','+str(lat)+',0'
  line.append(slonlat+'</coordinates></Point></Placemark>')
  lonlat.append(slonlat)
if 'L' in NorH:
  line.append('<Placemark> <LineString> <coordinates>')
  for i in range(len(a)):
    line.append(lonlat[i]+'\n')
  line.append('</coordinates></LineString></Placemark>')
line.append('</Document></kml>')
with open(fname+'.kml','w') as f:
  [f.write(l) for l in line]

