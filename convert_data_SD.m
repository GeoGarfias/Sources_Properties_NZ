% script to convert seismograms into event files like Xiaowei:


load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/cluster_SD.mat
load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/eqinfo_SD.mat

%for i=1;
for i=48:length(main_events);
    load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/cluster_SD.mat
    load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/eqinfo_SD.mat
    mevent = ['ev' num2str(main_events(i))]
    if ismember(mevent,mname_cat,'rows') == 1;
        clusterj_folder = ['/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/DATA/' mevent(3:16)];
        cd (clusterj_folder)
        cluster_file = [clusterj_folder '/SD_j' mevent(3:16) '_data.mat'];
        load (cluster_file)
        jlist=eval([ 'cat_' mevent '.index']);
        for j=jlist
            clear data
            event_names = num2str(eqinfo.name(j))
            xx = ['''*ev' event_names(1:12) '*SD''']
            event_i = eval([' char(who(' xx '))'])
            if isempty(event_i)  == 0
                for k=1:size(event_i,1)
                    clear smg_event
                    smg_event = eval(event_i(k,:));
                    clear station_name channel_name comp_number
                    
                    station_name = event_i(k,4:8);
                    channel_name = event_i(k,10);
                    
                    if strcmp(event_i(k,20:22),'AOI') == 1
                        station_name = event_i(k,20:end-2);
                        channel_name = event_i(k,end-1);
                    end
                    
                    if strcmp(channel_name,'Z') == 1
                        comp_number = 1
                    elseif strcmp(channel_name,'N') == 1
                        comp_number = 2
                    elseif strcmp(channel_name,'E') == 1
                        comp_number = 3
                    end
                    if strcmp(channel_name,'1') == 1
                        comp_number = 1
                    elseif strcmp(channel_name,'2') == 1
                        comp_number = 2
                    elseif strcmp(channel_name,'3') == 1
                        comp_number = 3
                    end
                    
                    
                    if k==1
                        eventj.id    =    eqinfo.id(j) ;
                        eventj.name    =    eqinfo.name(j,:) ;
                        eventj.qlat  =    eqinfo.qlat(j) ;
                        eventj.qlon  =    eqinfo.qlon(j) ;
                        eventj.qdep  =    eqinfo.qdep(j) ;
                        eventj.qyr   =    eqinfo.qyr(j) ;
                        eventj.qmon  =    eqinfo.qmon(j) ;
                        eventj.qdy   =    eqinfo.qdy(j) ;
                        eventj.qhr   =    eqinfo.qhr(j) ;
                        eventj.qmn   =    eqinfo.qmn(j) ;
                        eventj.qsc   =    eqinfo.qsc(j) ;
                        eventj.qmb   =    eqinfo.mb(j) ;
			eventj.mw    =    eqinfo.mw(j) ;
                    end % if k==1
                    
                    eventj.chan{k}= channel_name;
                    eventj.comp(k)= comp_number;
                    eventj.seis{k} = smg_event.data;
                    eventj.sta{k} = station_name;
                    eventj.mean(k) = smg_event.mean;
                    eventj.range(k) = smg_event.range;
                    eventj.dt(k) = smg_event.sampint;
                    eventj.samprate(k) = 1/smg_event.sampint;
                    eventj.B(k) = smg_event.B;
                    try
                        if smg_event.a > 0
                            eventj.a(k) = smg_event.a;
                        end
                    catch
                        % do nothing
                    end
                    try
                        if smg_event.T3 > 0
                            eventj.T3(k) = smg_event.T3;
                        end
                    catch
                        % do nothing
                    end
                    %eventj.KT0{k} = smg_event.KT0;
                    %eventj.T1(k) = smg_event.T1;
                    %eventj.KT1{k} = smg_event.KT1;
                    eventj.stlat(k) = smg_event.slat;
                    eventj.stlon(k) = smg_event.slon;
                    eventj.stelev(k) = smg_event.selev;
                    eventj.pazi(k) = smg_event.azim;
                    eventj.syr(k) = smg_event.year;
                    eventj.sjday(k) = smg_event.jday;
                    eventj.shr(k) = smg_event.hour;
                    eventj.smn(k) = smg_event.min;
                    eventj.ssc(k) = smg_event.sec;
                    
                    
                end % for k
            else
                continue
            end % if
            clear data
            
            data = eventj;
            
            eval([' save ev' num2str(eqinfo.name(j)) '.mat data'])
            clear eventj
        end % for j
    end % if ismember==1
    clearvars -except main_events
end % % length(main_events)













