% Complete new rewrite of stacking scripts.
% Started for OK Guthrie Nov 2016
% Being by compiling ALL ratios into one matrix with correct names
% for each row. 

% Dec 2017: added fix to remove P that got into S by mistake....
%
%   Agu 2020: ICJ-G for SAMBA DATA

close all;
clear all;
load /Users/home/juarezilma/Master/NETWORKS/SAMBA/EGF_matlab/WORK/cluster_SAMBA.mat;
%============================================================================================================
% set a default value of l in case no events..... (k per ratio)
k=0;PS=0;
xclimit = [0.5 0.6 0.7 0.8]; %Use with suffix 1a,
%list_event = [318 466]
for j=466;    
    % initiate some variables for counting number of ratios of
    clear nnsamp SAMBA* go_ahead  % related to sample decimation for STF
    eval(['load /Users/home/juarezilma/Master/NETWORKS/SAMBA/EGF_matlab/WORK/cluster_SAMBA.mat']);
    n100=0;  n200=0; n50=1; n80=1; n250=0;
    all_mname_list = eval(['mname_cat']);
    cluster = eval(['cat_' all_mname_list(j,:)]);
    clusterj_folder = ['/Users/home/juarezilma/Master/NETWORKS/SAMBA/EGF_matlab/DATA/Main_events/' cluster.name(1,3:16)];
    cd (clusterj_folder)
    nsec = nsec_value(cluster.mag(1))
    nsecX = nsecX_value(cluster.mag(1)); %Avoiding nsecX for the moment, making nsecX = nsec - icjg
    mag=cluster.mag(1) % Doing this to save it later - icjg
    mo = 10.^(((mag)*1.5)+9.1); % Doing this to save it later - icjg
    % Load in data:
    ['mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_2_egf_SN9.mat']
    if exist(['mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_2_egf_SN9.mat'])==2
        % make sure do not reset j when load in data file.....
        clear j_current mo_current nsec_current mag_current nsecX_current;
        j_current = j; mag_current = mag; mo_current = mo; nsec_current = nsec; nsecX_current = nsecX;
        eval(['load mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_2_egf_SN9.mat'])        
        eval(['load /Users/home/juarezilma/Master/NETWORKS/SAMBA/EGF_matlab/WORK/cluster_SAMBA.mat']);
        cd (clusterj_folder)
        j = j_current; mag=mag_current; mo=mo_current; nsec = nsec_current;nsecX = nsecX_current;
        all_mname_list(j,:)  
        % set limit for calculating LP level to normalise spectral ratios.
        flimit = 10/nsec;
        cluster = eval(['cat_' all_mname_list(j,:)]);
        %  For P and S: PS=1 for P, PS=2 for S
        PS_letter = ['P';'S'];
        for PS = 1:2
            % first get list of ratios in file ending in PS
            ratio_list_name = ['SAMBA_*' PS_letter(PS)];
            ratio_list_cell = who(ratio_list_name);           
            % Now can continue as for SCSN, NZ etc.
            % convert to a character array. -
            ratio_list = cell2mat(ratio_list_cell);
            [rnlen1 rnlen2] = size(ratio_list);
            % Check anything in list....
            if rnlen1 > 0
                % example nsgm: SAMBA_COSAX_Z_ev20120510115874_ev20151006041526_P
                mainjPS.hypo.mag = cluster.mag(1);
                mainjPS.hypo.mo = mo;
                mainjPS.hypo.nsec = nsec;
                mainjPS.hypo.nsecX = nsecX;
                mainjPS.hypo.lat = cluster.lat(1);
                mainjPS.hypo.lon = cluster.lon(1);
                mainjPS.hypo.depth = cluster.depth(1);
                mainjPS.hypo.datenumber = cluster.datenumber(1);               
                % k = all ratios, kk = all with f_resamp.
                kk=1;
                for k = 1:rnlen1
                    clear nsgm sgm
                    nsgm = ratio_list_cell{k};
                    sgm = eval(nsgm);
                    % only proceed if resampled....
                    if isfield(sgm,'f_resamp') == 1
                        %%% Note IF S actually has picks in data file
                        % this piece of code is written to find any S phases that are really P hangovers...
                        % December 2017
                        %%%%%%%%%%%%%%%%%%%
                        if PS == 1
                            % use data
                            mainjPS.all.stn(kk,:) = nsgm(7:11) ;  % Station
                            mainjPS.all.comp(kk,:) = nsgm(13) ;  % comp
                            mainjPS.all.egf(kk,:) = nsgm(end-15:end-2)  ; % EGF
                            mainjPS.all.f_resamp(kk,:) = sgm.f_resamp;
                            mainjPS.all.fratio_resamp(kk,:) = sgm.fratio_resamp;
                            % and normalise...
                            nflimit = 1;
                            for nfl=1:length(sgm.f_resamp)
                                % So USE nanmean of everything UNDER flimit
                                if sgm.f_resamp(nfl) < flimit
                                    nflimit = nflimit+1;
                                end
                            end   % end of nfl calculate loop
                            nflimit = nflimit-1;   % want to be one under not one over
                            norm_factor = nanmean(mainjPS.all.fratio_resamp(kk,1:nflimit));
                            mainjPS.all.norm_fratio_resamp(kk,:) = mainjPS.all.fratio_resamp(kk,:) ./ norm_factor;
                            clear norm_factor
                            counter(kk) = sgm.samp;
                            % deal with mixed sampling rates for stf
                            if isfield(mainjPS.all,'stf')==1
                                % ie if not first station....
                                mainjPS.all.stf(kk,:)=NaN*mainjPS.all.stf(kk-1,:);
                                mainjPS.all.norm_stf(kk,:)=NaN*mainjPS.all.stf(kk-1,:);
                            end  % need to put in NaNs for missing STFs.
                            %if isfield(sgm,'stf') == 1 && sgm.samp(1) == 100 % ignore low samp rate stations....
                            if isfield(sgm,'stf') == 1 && sgm.samp(1) > 50 % ignore low samp rate stations....
                                % Using the nex function to change the sample rate to 100 Hz for everything - icjg
                                % Change it from 50, 80, 200 to 100 - icjg
                                mainjPS.all.stf(kk,:) = sample_rate(sgm.samp,sgm.stf,n100,n200,n50,n80,n250)             
                                mainjPS.all.norm_stf(kk,:) = mainjPS.all.stf(kk,:) ./ range(mainjPS.all.stf(kk,:));
                            end   % if stf exists
                            %----------------------------------------------------------------------------------------------------
                            % Get xcorrelation.
                            [mxcoraf3 mxcorabf3] = max(sgm.xcorr2f3);
                            mainjPS.all.xcorrf3(kk) = mxcoraf3(1);
                            clear mx*
                            % get station info.
                            mainjPS.all.pazi(kk) = sgm.pazi;
                            mainjPS.all.sdist(kk) = sgm.del;
                            mainjPS.all.epi_d(kk) = sgm.epi_d;
                            mainjPS.all.slat(kk) = sgm.slat;
                            mainjPS.all.slon(kk) = sgm.slon;
                            % Not sure what pazi is so recalculate here.
                            [d1 a1]=distance(sgm.elat(1),sgm.elon(1),sgm.slat,sgm.slon);
                            mainjPS.all.sazi(kk) = a1;
                            clear d1 a1
                            mainjPS.all.egf_ml(kk) = sgm.ml(2);
                            mainjPS.all.egf_lat(kk) = sgm.elat(2);
                            mainjPS.all.egf_lon(kk) = sgm.elon(2);
                            mainjPS.all.samp(kk) = sgm.samp;
                            kk=kk+1;
                        end %% if PS =1
                        if PS == 2
                            % use data
                            mainjPS.all.stn(kk,:) = nsgm(7:11) ;  % Station
                            mainjPS.all.comp(kk,:) = nsgm(13) ;  % comp
                            mainjPS.all.egf(kk,:) = nsgm(end-15:end-2)  ; % EGF                          
                            mainjPS.all.f_resamp(kk,:) = sgm.f_resamp;
                            mainjPS.all.fratio_resamp(kk,:) = sgm.fratio_resamp;
                            % and normalise...
                            nflimit = 1;
                            for nfl=1:length(sgm.f_resamp)
                                % So USE nanmean of everything UNDER flimit
                                if sgm.f_resamp(nfl) < flimit
                                    nflimit = nflimit+1;
                                end
                            end   % end of nfl calculate loop
                            nflimit = nflimit-1;   % want to be one under not one over
                            norm_factor = nanmean(mainjPS.all.fratio_resamp(kk,1:nflimit));
                            mainjPS.all.norm_fratio_resamp(kk,:) = mainjPS.all.fratio_resamp(kk,:) ./ norm_factor;
                            clear norm_factor
                            counter(kk) = sgm.samp;
                            % deal with mixed sampling rates for stf
                            if isfield(mainjPS.all,'stf')==1
                                % ie if not first station....
                                mainjPS.all.stf(kk,:)=NaN*mainjPS.all.stf(kk-1,:);
                                mainjPS.all.norm_stf(kk,:)=NaN*mainjPS.all.stf(kk-1,:);
                            end  % need to put in NaNs for missing STFs.  
                            %if isfield(sgm,'stf') == 1 && sgm.samp(1) == 100 % ignore low samp rate stations....
                            if isfield(sgm,'stf') == 1 && sgm.samp(1) > 50 % ignore low samp rate stations....
                                 % Using the nex function to change the sample rate to 100 Hz for everything - icjg
                                % Change it from 50, 80, 200 to 100 - icjg
                                mainjPS.all.stf(kk,:) = sample_rate(sgm.samp,sgm.stf,n100,n200,n50,n80,n250)             
                                mainjPS.all.norm_stf(kk,:) = mainjPS.all.stf(kk,:) ./ range(mainjPS.all.stf(kk,:));
                            end   % if stf exists
                            % get xcorrelation.
                            [mxcoraf3 mxcorabf3] = max(sgm.xcorr2f3);
                            mainjPS.all.xcorrf3(kk) = mxcoraf3(1);
                            clear mx*
                            % get station info.
                            mainjPS.all.pazi(kk) = sgm.pazi;
                            mainjPS.all.sdist(kk) = sgm.del;
                            mainjPS.all.epi_d(kk) = sgm.epi_d;
                            mainjPS.all.slat(kk) = sgm.slat;
                            mainjPS.all.slon(kk) = sgm.slon;
                            % Not sure what pazi is so recalculate here.
                            [d1 a1]=distance(sgm.elat(1),sgm.elon(1),sgm.slat,sgm.slon);
                            mainjPS.all.sazi(kk) = a1;
                            clear d1 a1
                            mainjPS.all.egf_ml(kk) = sgm.ml(2);
                            mainjPS.all.egf_lat(kk) = sgm.elat(2);
                            mainjPS.all.egf_lon(kk) = sgm.elon(2);
                            mainjPS.all.samp(kk) = sgm.samp;
                            kk=kk+1;
                            clear data
                        end % if PS == 2
                    end % if resampled.
                end % for k over ratios
            end % if any ratios exist
            % and save resulting structure to proper name
            % save number of different samples rates:
            mainjPS.all.n100=n100;
            % and get number of different stations and egfs and their names:
            % And check whether there really are any stations to continue with...
            if isfield(mainjPS.all,'stn')==1
                [xs ys zs]=unique(mainjPS.all.stn,'rows');
                mainjPS.all.stnN = length(ys); %clear xs ys zs
                mainjPS.all.stn_list = xs;   % list of unique stations
                mainjPS.all.stn_azi = mainjPS.all.sazi(ys)';  % corresponding
                % list of azimuths
                mainjPS.all.stn_sdist = mainjPS.all.sdist(ys)'; % corresponding
                % list of distances
                [xe ye ze]=unique(mainjPS.all.egf,'rows');
                mainjPS.all.egfN = length(ye); %clear xe ye ze
                mainjPS.all.egf_list = xe;  % list of unique egfs
                clear xs ys zs xe ye ze
                % Now get ratios per station
                for l = 1:mainjPS.all.stnN
                    kk=1; % counter for ratios.
                    for k=1:length(mainjPS.all.samp)
                        if strcmp(mainjPS.all.stn(k,:),mainjPS.all.stn_list(l,:))==1
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.sazi = mainjPS.all.sazi(k);']);
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.slat = mainjPS.all.slat(k);']);
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.slon = mainjPS.all.slon(k);']);
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.comp(kk,:) = mainjPS.all.comp(k,:);']);
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.egf(kk,:) = mainjPS.all.egf(k,:);']);
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.f_resamp(kk,:) = mainjPS.all.f_resamp(k,:);']);
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.fratio_resamp(kk,:) = mainjPS.all.fratio_resamp(k,:);']);
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.norm_fratio_resamp(kk,:) = mainjPS.all.norm_fratio_resamp(k,:);']);
                            eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.xcorrf3(kk) = mainjPS.all.xcorrf3(k);']);
                            if isfield(mainjPS.all,'stf')==1
                                eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.stf(kk,:) = mainjPS.all.stf(k,:);']);
                                eval(['mainjPS.stn.' mainjPS.all.stn_list(l,:) '.norm_stf(kk,:) = mainjPS.all.norm_stf(k,:);']);
                            end
                            kk=kk+1;
                        end  % if stn matches
                    end % for k ratios
                end % for l stations.      
                % Now get ratios per egf
                for l = 1:mainjPS.all.egfN
                    kk=1; % counter for ratios.
                    for k=1:length(mainjPS.all.samp)
                        if strcmp(mainjPS.all.egf(k,:),mainjPS.all.egf_list(l,:))==1
                            eval(['mainjPS.egf.ev' mainjPS.all.egf_list(l,:) '.comp(kk,:) = mainjPS.all.comp(k,:);']);
                            eval(['mainjPS.egf.ev' mainjPS.all.egf_list(l,:) '.stn(kk,:) = mainjPS.all.stn(k,:);']);
                            eval(['mainjPS.egf.ev' mainjPS.all.egf_list(l,:) '.f_resamp(kk,:) = mainjPS.all.f_resamp(k,:);']);
                            eval(['mainjPS.egf.ev' mainjPS.all.egf_list(l,:) '.fratio_resamp(kk,:) = mainjPS.all.fratio_resamp(k,:);']);
                            eval(['mainjPS.egf.ev' mainjPS.all.egf_list(l,:) '.norm_fratio_resamp(kk,:) = mainjPS.all.norm_fratio_resamp(k,:);']);
                            eval(['mainjPS.egf.ev' mainjPS.all.egf_list(l,:) '.xcorrf3(kk) = mainjPS.all.xcorrf3(k);']);
                            if isfield(mainjPS.all,'stf')==1
                                eval(['mainjPS.egf.ev' mainjPS.all.egf_list(l,:) '.stf(kk,:) = mainjPS.all.stf(k,:);']);
                                eval(['mainjPS.egf.ev' mainjPS.all.egf_list(l,:) '.norm_stf(kk,:) = mainjPS.all.norm_stf(k,:);']);
                            end
                            kk=kk+1;
                        end  % if egf matches
                    end % for k ratios
                end % for l egfs.
            end % end of if actually any ratios in all....
            eval([all_mname_list(j,:) '_' PS_letter(PS) ' = mainjPS;']);
            clear mainjPS
        end % for PS
        
        eval(['save stack_' cluster.name(1,:) '_1_SN9.mat  ev*'])
        clear e* n* r* x* s* m* l k* i* h* f* d* c* b* a* PS* OK* Ny* FRANGE* j_input j_current
    end % if input file exists for this event
end % for j over main
%_________________________________________________


