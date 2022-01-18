# Script made by ICJuarez-Garfias
# Script  to read all the waveforms and make a nordic file
# For Southern Dwarfs 06/08/20
#==================================================================================================================
# M O D U L E S 
import numpy as np
import pandas as pd
import datetime, os
from obspy.core import UTCDateTime
from obspy.core.stream import read
from obspy import read_events, Stream
from obspy.io.nordic.core import write_select, _write_nordic
from obspy.clients.fdsn import Client as FDSN_Client
from obspy.core.event import Event, Catalog
#==================================================================================================================
# F U N C T I O N S

def egfs_dataframe(egfs_data,complete_cat):
    # This function create a new data frame with all the information of all the possible EGF of a Target event.
    # Compares a list of IDs with a complete csv catalog
    # The new data frame is saved in the EGF_events variable.
    EGF_events = pd.DataFrame()
    for index_1, row_1 in egfs_data.iterrows():
        for index_2, row_2 in complete_cat.iterrows():
            if row_1[0]==row_2['id']:
                for column_name, column_data in row_2.iteritems():
                    if not column_name in row_1:
                        row_1[column_name] = column_data
                EGF_events = EGF_events.append(row_1, ignore_index=True)
    return EGF_events
# End function   


def find_event(cat, egf_time):
    event = Event()
    for i in range(len(cat)):
        if abs(cat[i].origins[-1].time - egf_time) < 1:
            event = cat[i]
        else:
            continue
    return event

def get_nordic(egf_table, t_event, xml_cat):
    cat = Catalog()
    wave_names = []
    # The easiest way to solve this ----------------------------------------------------------------------------------------------------------
    SAMBA_list = ['COSA', 'COVA', 'EORO', 'FRAN', 'GOVA', 'LABE', 'LARB', 'MTBA', 'MTFO', 'POCR2', 'SOLU', 'WHAT2', 'WHYM']
    ND_list = ['MOSQ', 'BELL', 'JACK', 'HURA', 'SHOT', 'MORG', 'TURI2', 'BRUN', 'GLOR']
    SD_list = ['OLIV', 'CASC', 'DATA', 'GORG', 'BARN', 'RICH', 'DELT', 'THOM', 'MONK', 'POMM']
    COSA_list = ['TEPE', 'HUVA', 'STBA', 'NOBU', 'ASPR', 'MORV', 'KING', 'HAAS']
    Alfa_list = ['BURA2', 'CAMEL', 'CASH', 'FBLA2', 'HATT', 'MTHA2', 'TURI', 'UMAT'] 
    Wizard_list = ['WZ01', 'WZ02', 'WZ03', 'WZ04', 'WZ05', 'WZ06', 'WZ07', 'WZ08', 'WZ09', 'WZ10', 'WZ11', 'WZ12', 'WZ13', 'WZ14', 'WZ15', 'WZ16', 'WZ17', 'WZ18', 'WZ19', 'WZ20']
    GeoN_list = ['JCZ', 'WKZ', 'MSZ', 'MLZ', 'RPZ', 'LBZ', 'EAZ', 'FOZ', 'WVZ', 'INZ', 'LTZ', 'DSZ', 'THZ', 'OXZ', 'NSBS', 'HDWS', 'MECS', 'LPLS', 'WNPS', 'GLNS', 'QTPS', 'TAFS', 'TWAS', 'MCNS', 'FGPS', 'FJDS', 'WHFS', 'WHAS', 'WHFS', 'HAFS', 'KOKS', 'APPS', 'HMCS', 'GMFS', 'ARPS', 'IFPS', 'SJFS', 'GLWS', 'TKAS', 'RDCS', 'WBCS', 'INGS', 'CSHS']
    # ----------------------------------------------------------------------------------------------------------------------------------------
    for index1, egf in egf_table.iterrows():
        st = Stream()
        date = datetime.datetime(int(egf['qyr']),int(egf['qmon']),int(egf['qdy']),int(egf['qhr']),int(egf['qmn']),int(egf['qsc']))
        JDay = date.strftime('%j')
        #Change to UTCDateTime
        eventdt = UTCDateTime(date)
        print('----------------------------------------------------------------------------------')
        print('----------------------------------------------------------------------------------')
        print('Working in the EGF event: ', eventdt)
        print('----------------------------------------------------------------------------------')
        print('----------------------------------------------------------------------------------')
        # Time to get the waveforms
        print('========================================== NOW IM READING THE NETWORKS DATA  =======================================')
        egf_event = find_event(xml_cat, eventdt)
        stat_list = []
        for j in range(len(egf_event.picks)):
            stat_list.append(egf_event.picks[j].waveform_id.station_code)
        stat_list = set(stat_list)
        print(stat_list)
        for stat in stat_list:
            if stat in SD_list or stat in ND_list:
                print('G E T T I N G   S O U T H   D W A R F S   DATA')
                try:
                    st += read('/Volumes/GeoPhysics_11/users-data/chambeca/DWARFS/DWARFS_archive/Y' +  eventdt.strftime('%Y') + '/R' + JDay + '.01/' + stat + '*',starttime=eventdt, endtime=eventdt+50)
                except:
                    print('NO DATA')
                    pass
            elif stat in COSA_list:
                print('G E T T I N G   C O S A   DATA')
                try:
                    st += read('/Volumes/GeoPhysics_11/users-data/chambeca/COSA_archive/Y' +  eventdt.strftime('%Y') + '/R' + JDay + '.01/' + stat + '*',starttime=eventdt, endtime=eventdt+50)
                except:
                    pass
            elif stat in SAMBA_list:
                print('G E T T I N G   S A M B A   DATA')
                st += read('/Volumes/GeoPhysics_09/users-data/chambeca/SAMBA_archive/day_volumes_S/Y' +  eventdt.strftime('%Y') + '/R' + JDay + '.01/' + stat + '*',starttime=eventdt, endtime=eventdt+50)
            elif stat in GeoN_list:
                print('G E T T I N G   G E O N E T   DATA')
                client = FDSN_Client("GEONET")
                client_nrt = FDSN_Client("https://service-nrt.geonet.org.nz")
                try:
                    st += client.get_waveforms(network="NZ", station=stat, location='*', channel='HH*', starttime=eventdt, endtime=eventdt+50) 
                except:
                    pass
        print('----------------------------------------------------------------------------------')
        print('The traces are:')
        print(st)
        print('----------------------------------------------------------------------------------')
        for tr in st:
            tr.data = tr.data.astype(np.int32)
        for tr in st:
            if type(tr.data) == np.ma.core.MaskedArray:
                try:
                    print('Masked data found for ' + tr.stats.station + ' ' + tr.stats.channel +
                        ': padding empty spaces with zeros')
                    tr.data = tr.split().detrend('simple').merge(fill_value=0)[0].data
                except: 
                    st.remove(tr)
        st.write('/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/Seisan_files/Target_Events/' + t_event[0:14] + '/Wav_files/' + eventdt.strftime('%Y') + eventdt.strftime('%m') + eventdt.strftime('%d') + eventdt.strftime('%H') + eventdt.strftime('%M') + eventdt.strftime('%S') + '_SD_multiplexed', format='MSEED')    
        _write_nordic(event=egf_event, outdir='/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/Seisan_files/Target_Events/' + t_event[0:14] + '/S_files/', filename=None, userid="ICJG", evtype="L", wavefiles='/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/Seisan_files/Target_Events/'+ t_event[0:14] + '/Wav_files/' + eventdt.strftime('%Y') + eventdt.strftime('%m') + eventdt.strftime('%d') + eventdt.strftime('%H') + eventdt.strftime('%M') + eventdt.strftime('%S') + '_SD_multiplexed')
        # To creat the select.out file
        wave_list = eventdt.strftime('%Y') + eventdt.strftime('%m') + eventdt.strftime('%d') + eventdt.strftime('%H') + eventdt.strftime('%M') + eventdt.strftime('%S') + '_SD_multiplexed'
        wave_names.append(wave_list)
        cat.append(egf_event)
        #END FOT EGF 
    write_select(catalog=cat, filename='/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/Seisan_files/Target_Events/' + t_event[0:14] + '/select_' + t_event[0:14] + '.out', userid='ICJG', evtype='L', wavefiles=wave_names)
    # End function 

#==================================================================================================================
# M A I N  C O D E
# Inputs --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# The line below makes a list of all the files at the indicated path. The output is a list of the files 'event_EGF_list.csv', those files includes the id of the EGF events of a target event.
if __name__ == '__main__':
    print('The process starts')
    #id_list = np.loadtxt('/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/Seisan_files/lista.txt')
    #id_list = id_list.tolist()
    id_list = ['20190523114632','20190709071545']
    #id_list =  os.listdir('/Users/home/juarezilma/Master/NETWORKS/DWARFS_N/Seisan_files/All_list/') # Read the .txt files with the id information of the EGF
    print(id_list)
    sdwarfs_cat = pd.read_csv('/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/Catalogue/eqinfo_SD.csv') # The complete catalogue in csv version
    catalog = read_events('/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/Catalogue/DWARS_cat575cc_hypoDD_mags.xml') # The complete catalogue in xml version
    print('I just read the xml file')

    for target_event in id_list:
        target_event = str(target_event)
        print('==================================================================================================')
        print('==================================================================================================')
        print('==================================================================================================')
        print('Working with the main event:', target_event[0:14])
        print('==================================================================================================')
        print('==================================================================================================')
        print('==================================================================================================')
        egf_id = pd.read_csv('/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/Seisan_files/All_list/' + target_event[0:14] + '_EGF_list.csv', header=None) #List of id
        # This function create a new data frame with all the information of all the possible EGF of a Target event
        egf_events_list = egfs_dataframe(egf_id,sdwarfs_cat)
        print(egf_events_list)
        print('==================================================================================================')
        print('==================================================================================================')
        # Now its time to generate the nordic files
        # Calling the function
        get_nordic(egf_events_list,target_event,catalog)
    # End code
#==================================================================================================================