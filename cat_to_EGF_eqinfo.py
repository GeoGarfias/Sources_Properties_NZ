# code to write out obspy catalog object as a text file in the format needed for EGF codes
# written out as .mat file to be read into matlab
# date 22/11/19
# author EW-S

# MODIFIED BY IC JUAREZ-GARFIAS
from obspy import read_events
import csv

# can give either a seisan select outfile, or an obspy xml catalog
select_file='/Users/home/juarezilma/Master/NETWORKS/DWARFS_N/Catalogue/DWARN_Mags_uncerts.xml'

catalog=read_events(select_file)  #either read in the QuakeML file or the select.out file from Seisan. 
#The name of the output file
outfile='/Users/home/juarezilma/Master/NETWORKS/DWARFS_N/Catalogue/eqinfo_ND.csv'

qid = []
qlat = []
qlon = []
qyr = []
qmon = []
qday = []
qhr = []
qmn = []
qsc = []
ml = []
qdep = []
name = []
# Add this if you want the focal mechanisms
# strike = []
# dip = []
# rake = []



for i, event in enumerate(catalog):
    qid.append(str(i))
    qlat.append(event.origins[0].latitude)
    qlon.append(event.origins[0].longitude)
    qyr.append(event.origins[0].time.strftime('%Y'))
    qmon.append(event.origins[0].time.strftime('%m'))
    qday.append(event.origins[0].time.strftime('%d'))
    qhr.append(event.origins[0].time.strftime('%H'))
    qmn.append(event.origins[0].time.strftime('%M'))
    qsc.append(event.origins[0].time.strftime('%S'))
    ml.append(event.magnitudes[0].mag)
    qdep.append(event.origins[0].depth/1000)
    # if len(event.magnitudes) > 0:
    #     ml.append(event.magnitudes[-1].mag)
    # else:
    #     ml.append('0')
    # qdep.append(event.origins[0].depth/1000)
    # Add this if you want to add the focal mechanisms to the eqinfo.csv file
    # if len(event.focal_mechanisms) == 1:
    #     strike.append(event.focal_mechanisms[0].nodal_planes.nodal_plane_1.strike)
    #     dip.append(event.focal_mechanisms[0].nodal_planes.nodal_plane_1.dip)
    #     rake.append(event.focal_mechanisms[0].nodal_planes.nodal_plane_1.rake)
    # else:
    #     strike.append('0')
    #     dip.append('0')
    #     rake.append('0')
    name.append(str(event.origins[0].time.strftime('%Y%m%d%H%M%S')))
    
rows=zip(qid, qlat, qlon, qyr, qmon, qday, qhr, qmn, qsc, ml, qdep, name)
#rows=zip(qid, qlat, qlon, qyr, qmon, qday, qhr, qmn, qsc, ml, qdep, name)

with open(outfile, mode='w') as f:
    writer=csv.writer(f)
    writer.writerow(['id', 'qlat', 'qlon', 'qyr', 'qmon', 'qdy', 'qhr', 'qmn', 'qsc', 'ml', 'qdep', 'name'])
with open(outfile, mode='a') as f:  
    writer=csv.writer(f)
    for row in rows:
        writer.writerow(row)
    
