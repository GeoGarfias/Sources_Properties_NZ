    % Script to read in sac files
% add in path for rsac.m

% must be more efficient way to do this e.g.:
%
clear all
close all

load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/cluster_SD.mat
load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/eqinfo_SD.mat

%for m=15;
for m=28:length(main_events);
    mevent = ['ev' num2str(main_events(m))]
    %if ismember(mevent,mname_cat,'rows') == 1;
    load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/cluster_SD.mat
    load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/eqinfo_SD.mat
        hlist=eval(['cat_' mevent '.index'])
        clusterj_folder = ['/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/DATA/' mevent(3:16)];
        cd (clusterj_folder)
        cd SAC_files 
        for q=hlist
        %cd SAC_files
            myFiles = dir(fullfile(strcat('*.sac_cut'))); %gets all wav files in
            fid=fopen('a_file_list.m', 'wt')'
            %writes results out to text file b/c matlab sucks and cant assign variables properly
            for k = 1:length(myFiles)
                baseFileName = myFiles(k);
                fprintf(1, 'Now reading %s\n', baseFileName.name);
                FileName=strsplit(baseFileName.name, '.');
                test=join([string(FileName(1,1)),'_SD=rsac(' string(baseFileName.name) ');' '\n']);
                test1=strrep(test, ' ', '');
                test2=strrep(test1, '(', '("');
                test3=strrep(test2, ')', '")');
                fprintf(fid, test3);
            end
            fclose(fid)
            a_file_list %then run the file
            % then need to evaluate contents of a_file_list.txt 
        end
        % make a list of ALL the seismograms read in with rsac.
        % so find some part of name they all have in common eg. 2016
        xx = ['''*_SD*''']
        record_list =eval([' char(who(' xx '))'])

        % For each of these files, get into useable format.
        % Includes Magnitude as mag.
        cd ..
        stnme = record_list;
        [len1 len2] = size(stnme);
        for ii=1:len1  % records
            strnme=stnme(ii,:);
            tempvar=eval([stnme(ii,:)]);
            eval(['clear ' stnme(ii,:)])
            
            % picking record to start at event origin time
            B=tempvar(6,3);
            E=tempvar(7,3);
            NPTS=tempvar(80,3);
            DELTA=tempvar(1,3);
            
            %% For PiLAB just want whole record.
            I=1;

            SMG=tempvar(I:NPTS,2);
            eval([strnme '.data=SMG;']);
            eval([strnme '.mean=mean(SMG);']);
            eval([strnme '.range=max(SMG)-min(SMG);']);
            eval([strnme '.numdata=length(SMG);']);
            eval([strnme '.sampint=DELTA;']);
            eval([strnme '.B=B;']);
        
            % picking P and S times
            % Adjust this depending on what headers are filled in as what in
            % your sac files.
            a=tempvar(9,3);
            if a > 0
                eval([strnme '.a=a;']);
                KA = char(tempvar(151,3));
                eval([strnme '.KA=KA;']);
            end
            T0=tempvar(11,3);
            if T0 > 0
                eval([strnme '.T0=T0;']);
                KT0 = char(tempvar(159,3));
                eval([strnme '.KT0=KT0;']);
            end
            T1=tempvar(12,3);
            if T1 > 0
                eval([strnme '.T1=T1;']);
                KT1 = char(tempvar(167,3));
                eval([strnme '.KT1=KT1;']);
            end
            T2=tempvar(13,3);
            if T2 > 0
                eval([strnme '.T2=T2;']);
                KT2 = char(tempvar(175,3));
                eval([strnme '.KT2=KT2;']);
            end
            T3=tempvar(14,3);
            if T3 > 0
                eval([strnme '.T3=T3;']);
            end
            %%
            % coord
            if tempvar(32,3)~= -12345
                eval([strnme '.slat=tempvar(32,3);'])
                eval([strnme '.slon=tempvar(33,3);']);
                eval([strnme '.selev=tempvar(34,3);']);
            end
            eval([strnme '.elat=tempvar(36,3);']);
            eval([strnme '.elon=tempvar(37,3);']);
            eval([strnme '.edep=tempvar(39,3);']);
            eval([strnme '.mag=tempvar(40,3);']);
            if tempvar(51,3)~=-12345
                dist1=tempvar(51,3);
                azi1=tempvar(52,3);
                bazi1=tempvar(53,3);
                gcarc1=tempvar(54,3);
            else
                if tempvar(32,3)~= -12345
                    e_strnme = eval(strnme);
                    [gcarc1 azi1] = distance(e_strnme.elat,e_strnme.elon,e_strnme.slat,e_strnme.slon);
                    bazi1 = azimuth(e_strnme.slat,e_strnme.slon,e_strnme.elat,e_strnme.elon);
                    dist1 = gcarc1*111.1;
                else
                    dist1  = NaN;
                    azi1   = NaN;
                    bazi1  = NaN;
                    gcarc1 = NaN;
                end
            end
            eval([strnme '.dist=dist1;']);
            eval([strnme '.azim=azi1;']);
            eval([strnme '.bazim=bazi1;']);
            eval([strnme '.gcarc=gcarc1;']);
            % picking dates and times
            eval([strnme '.year=tempvar(71,3);']);
            eval([strnme '.jday=tempvar(72,3);']);
            eval([strnme '.hour=tempvar(73,3);']);
            eval([strnme '.min=tempvar(74,3);']);
            eval([strnme '.sec=tempvar(75,3);']);
            eval([strnme '.msec=tempvar(76,3);']);
        
            % Now save some variables for using later.
            azi(ii) = azi1;
            dist(ii) = dist1;
            gcarc(ii) = gcarc1;
        
            clear e_strnme gcarc1 dist1 bazi1 azi1
        
        end;
        eval(['save SD_j' num2str(mevent) '_data.mat *SD'])
    %end % if ismember==1
    clearvars -except main_events
end % length(main_events)


