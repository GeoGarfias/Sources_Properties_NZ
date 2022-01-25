% First time trying to make this to work better - ICJ-Garfias 20/08/20
% ========================================================================================================
% Rachel Abercrombie.
% New version Nov 2016 to not read in ALL EGFS at once, and only
% Version using Germans matlab code deconvolution.
% ========================================================================================================
% P R E L I M I N A R S
% Before using need to have run make_egf_clusters which makes lists main and EGFs within
% Needs these matlab routines to run:
%  adaptspec.m
%  adaptspec2.m
%  avgspec.m
%  avgsspec2.m
%  bpassc.m 
%  CLIP_DETECT.m
%  create_fvector.m
%  eigenspec.m
%  mt_deconv.m
%  mtspec.m
%  qiinv.m
%  nsec_value.m
%  save_main_data.m
%  save_egf_data.m
%  start_signal.m
%  check_samplerate.m
%
% ========================================================================================================
% save seismograms in mname_data that have spectra and STFs
% Now loops over the list of seismograms in the file, NOT over the list of stations and components.....
% Perform MT deconvolution of pulses for EGF
%  uses same index for egf as main so mname1 is a list of the main events
%. Also deals with sampling rate differences. First interpolates to highest
%  rate, and then keeps minimum sampling rate for max real f.
% ========================================================================================================
% ========================================================================================================
% C O D E   B E G I N
% Clear out any old variables
clear x* y* s* t* n*
clear out old ratios:
clear *MT

%figure(3);clf;
%---------------------------------------------------------------------------------------------------------- 
%----------------------------------------------------------------------------------------------------------
% If you are running the code for one event only, 'j' is set for that event. 
% You need to write the event matlab id that is set in the code ' make_clusters.m' and its magnitude dependant.

% Looping variables:
% j = loop over each main-egf cluster
% l = loop over events within EGF cluster
%        EIDe = EGF cluster
%        suffix j = main event, suffix l = EGF event.
% k = loop over all seismograms - Bill first then AA

% Delete previous run log.
!rm junk_egf.out
%----------------------------------------------------------------------------------------------------------
% Currently working for one event at a time
% List of events by id-magnitude 
load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/cluster_SD.mat

for j= 51
%for j=62:length(main_events) 
    load /Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/WORK/cluster_SD.mat
    mevent = ['ev' num2str(main_events(j))]
    % Im adding the next line because I did not save that variable when I run the 'make_clusters.m' code - ICJG
    all_mname_list = eval(['mname_cat']);
    cluster = eval(['cat_' mevent ';']); 
    clusterj_folder = ['/Users/home/juarezilma/Master/NETWORKS/DWARFS_S/EGF_matlab/DATA/' mevent(3:16)];
    save temp_cluster cluster
    clear cat* max* ml* counter msort* YHSbox  mname_j mname_ml mname_nEGF
    load temp_cluster
    eval([mevent '= cluster;']);

    cd (clusterj_folder) % set to work directory
    %_______________________________________________________________________________________________________________
    % Variables for trying to vary P and S picks:
    x=0;      % samples added to EGF P pick
    y=0;      % samples added to EGF S pick
    % Reset counters etc. for each new main event fname..
    % set counters for xcorrelation output
    count_xcorr0 = 1;
    count_xcorrx = 1;
    % Set pcnt for determining start time of window before P or S..
    pcnt = 10; 
    dummy = 1;   % For checking on sample rates....
    % You can add a correlation limit to make it faster - ICJG
    xcor_limit = 0.0  
    % Set vs and vp for approx. S-P time calculation...
    vp = 5;
    vs = vp/sqrt(3);
    %_______________________________________________________________________________________________________________
    [mlen2 mlen1]=size(cluster.id)
    
    % If want to know how to calculate the nsec value, check the function
    % The cluster.mag is the mw magnitude calculate in the make_cluster
    % code
    nsec = nsec_value(cluster.mw(1))
    % Checking if the file already exist
    if exist(['mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_1_egf.mat'])==0    
        % make sure do not reset j when load in data file.....
        clear j_current;
        j_current = j;
        % Now load in data for main and EGFs and rename.
        cluster.name(1,:)
        clear input_fname
        input_fname = [num2str(cluster.name(1,:))],%keyboard
        % Check if the file 'evYYYMMDDHHmmss.mat' for the main event exist
        if exist([input_fname '.mat'])==2
            eval(['load ' input_fname '.mat'])
            eval([cluster.name(1,:) '= data; clear data;'])
            [mlen2 mlen1]=size(cluster.id)           
            j = j_current;         
            % So now have data in Xiaowei format. 
            % Need to get it to mine (RACHEL) so as to work with
            % Recalculate nsec, just in case. 
            nsec = nsec_value(cluster.mw(1))
            % Set the low pass filters... for cross correlation
            highpassf1 = 1.0
            lowpassf1 = 40    % set a high value
            lowpassf2 = 10    % set a low value  
            lowpassf3 = (10/nsec + highpassf1)
            mname = cluster.name
            [mlen1 mlen2]=size(cluster.name)
            
            %Now make lists of available seismograms for main event....
            cd (clusterj_folder)
            clear main_event
            main_event = eval([cluster.name(1,:)]);
            % Now start iterating over egfs.
            for l=1:mlen1 %original slow way.
                cd (clusterj_folder)
                clear input_fname
                input_fname = [num2str(cluster.name(l,:))],%keyboard
                % Check if the file 'evYYYMMDDHHmmss.mat' exist - THOSE ARE THE EGF OF THE MAIN EVENT
                if exist([input_fname '.mat'])==2
                    eval(['load ' input_fname '.mat'])
                    eval([cluster.name(l,:) '= data; clear data;'])
                    clear egf_event
                    egf_event = eval([cluster.name(l,:)]);

                    % M A K I N G   T H E   L O O P  over the traces of the M A I N   E V E N T
                    for k=1:length(main_event.syr)
                        clear tstp sgmlp sgmls sgmjp sgmjs sig_startj sig_startl
                        clear nsgmj nsgml

                        msta = main_event.sta{k};
                        if isempty(msta) ==0
                            mchan = main_event.chan{k};
                            nsgmj = [cluster.name(1,:) '_' msta '_' mchan]
                            clear sgmj
                            sgmj = save_main_data(main_event,k);
                        
                            main_index = k;                          
                            clear nsgml
                            % M A K I N G   T H E   L O O P   over the E G F
                            for l_e = 1:length(egf_event.syr)    % over stations/chans (ie seismograms)
                                if strcmp(egf_event.sta{l_e},msta)==1 && strcmp(egf_event.chan{l_e},mchan)==1
                                    clear sgml
                                    nsgml = [cluster.name(l,:) '_' cell2mat(egf_event.sta(l_e)) '_' cell2mat(egf_event.chan(l_e))]; 
                                    egf_index = l_e; 
                                    % This line called the function 'save_egf_data.m' - That function save all the needed info of the EGF event in a MATLAB structure - ICJG
                                    sgml = save_egf_data(egf_event,l_e);
                                    % Now make seismograms: sig_sgmj and sig_sgml and noi_sgmj,noi_sgml
                                    % initialise these to zero to make if loops easier later...
                                    %  Sort out sampling rate variables
                                    sampj = round(1./sgmj.sampint);
                                    % No point bothering with the really low sample rate data......
                                    if sampj > 20
                                        sampintj = 1/sampj;
                                        %______________________________________________________________________________________________
                                        % I will comment all the lines with filter f1 and f2 - dont need them at the moment ICJG
                                        %______________________________________________________________________________________________                                  
                                        % The signal for the main event is been filtered in the next 3 lines with the 3 different lowpass -ICJG
                                        % if lowpassf1 > sampj*0.4
                                        %     lowpassf1_use = sampj*0.4
                                        % else
                                        %     lowpassf1_use = lowpassf1;
                                        % end
                                        clear f1sgmj f2sgmj f3sgmj f1sig* f2sig* f3sig*
                                        % Filtering the main event signal
                                        % f1sgmj = bpassc(sgmj.data,1,length(sgmj.data),2,highpassf1,lowpassf1_use,sampj);
                                        % f2sgmj = bpassc(sgmj.data,1,length(sgmj.data),2,highpassf1,lowpassf2,sampj);
                                        f3sgmj = bpassc(sgmj.data,1,length(sgmj.data),2,highpassf1,lowpassf3,sampj);                
                                        samp_count(dummy) = sampj;
                                        dummy = dummy+1;                                      
                                        sampl = round(1./sgml.sampint);
                                        sampintl = 1/sampl;
                                        % Filtering the egf signal
                                        clear f1sgml f2sgml f3sgml
                                        % f1sgml = bpassc(sgml.data,1,length(sgml.data),2,highpassf1,lowpassf1_use,sampl);
                                        % f2sgml = bpassc(sgml.data,1,length(sgml.data),2,highpassf1,lowpassf2,sampl);
                                        f3sgml = bpassc(sgml.data,1,length(sgml.data),2,highpassf1,lowpassf3,sampl);
                                                                               
                                        % Variables for trying to vary P and S picks:
                                        x=0;      % samples added to EGF P pick
                                        y=0;      % samples added to EGF S pick
                                        
                                        % Determine whether P or S based on component....
                                        for PS=1:2
                                            PS_letter = ['P';'S'];                                    
                                            % correction line added Dec 2017. Anything pre this may have P in S waves..
                                            clear sig_startj sig_startl
                                            %% Added loop to redo with egf offset by xcorr offset.
                                            for kk = 1:2
                                                epi_d = cluster.epi_d(l);
                                                hyp_sep=NaN;
                                                %--------------------------------------------------------------------------------------
                                                % Only do P windows if tstp > nsec
                                                % No way of dealing with this at the moment......! PROBLEM                                        
                                                if PS == 1
                                                    if isfield(sgmj,'p') == 1    && ~isnan(sgmj.p)==1
                                                        if isfield(sgml,'p') == 1     && ~isnan(sgml.p)==1
                                                            tstp = sgmj.s - sgmj.p
                                                            %_________________________________________________
                                                            % Add this to use a smaller window for P - ICJG
                                                            if nsec == 5
                                                                nsec_p = nsec-2
                                                            else 
                                                                nsec_p = nsec
                                                            end
                                                            %_________________________________________________
                                                            if tstp >= nsec_p 
                                                                % I made the function start_signal.m to replace the next lines
                                                                %[sig_startj, sig_startl] = start_signal(sgmj, sgmj.p, sgml.p, nsec, sampj, sampl)
                                                                sgmjp=fix(sgmj.p / sgmj.sampint);
                                                                sgmlp=fix(sgml.p / sgmj.sampint) +x;
                                                                ic = fix(nsec*sampj * pcnt / 100);
                                                                sig_startj = round(sgmjp-ic);
                                                                ic = fix(nsec*sampl * pcnt / 100);
                                                                sig_startl = round(sgmlp-ic);                                                          
                                                            end   % S-P big enough.
                                                            %end % if S pick exists..
                                                        end   % isfield sgml.p
                                                    end   % isfield sggmj.p
                                                end   % if PS ==1                                           
                                                if PS == 2
                                                    if isfield(sgmj,'s') == 1     && ~isnan(sgmj.s)==1
                                                        if isfield(sgml,'s') == 1     && ~isnan(sgml.s)==1
                                                            %I made the function start_signal.m to replace the next lines
                                                            %[sig_startj, sig_startl] = start_signal(sgmj, sgmj.s, sgml.s, nsec, sampj, sampl)
                                                            sgmjs=fix(sgmj.s / sgmj.sampint);
                                                            sgmls=fix(sgml.s / sgmj.sampint) +y;
                                                            ic = fix(nsec*sampj * pcnt / 100);
                                                            sig_startj = round(sgmjs-ic);
                                                            ic = fix(nsec*sampl * pcnt / 100);
                                                            sig_startl = round(sgmls-ic);                                                      
                                                        end
                                                    end
                                                end  % if PS == 2
                                                %--------------------------------------------------------------------------------------
                                                % EWS added this line in to deal with too short seismograms
                                                clear check_l check_j
                                                nl  = fix(nsec*sampl);
                                                if (exist('sig_startj') == 1) && (exist('sig_startl') == 1)
                                                    check_l = sig_startl + nl - length(sgml.data)
                                                    check_j = sig_startj + nl - length(sgmj.data)
                                                end
                                                %--------------------------------------------------------------------------------------
                                                % Continue with EGF deconvolution if have a sig_startj and a sig_startl
                                                if (exist('sig_startj') == 1) && (exist('sig_startl') == 1 && check_l <= 0 && check_j <= 0)
                                                    noi_startj = round(sig_startj - sampj*nsec*1.05);                                               
                                                    nj  = fix(nsec*sampj);             % number of samples  
                                                    % Put in to deal with picks beyond length of seismogram when P missed        
                                                    if noi_startj < length(sgmj.data) && noi_startj+ nj > 5                                                       
                                                        if noi_startj < 1
                                                            nzero = zeros(abs(noi_startj)+1,1);
                                                            noi_sgmj = [nzero; sgmj.data(1:noi_startj+nj)-mean(sgmj.data(1:noi_startj+nj))];
                                                            noi_sgmj = noi_sgmj'
                                                        else
                                                            noij = noi_startj: noi_startj+nj;  % length for seismogram-noise
                                                            noi_sgmj = sgmj.data(noij)-mean(sgmj.data(noij));        % main seismogram - noise
                                                        end   % if noi_startj < 1                                                      
                                                        if sig_startj < 1
                                                            nzero = zeros(abs(sig_startj)+1,1)';                                                          
                                                            sig_sgmj = [nzero sgmj.data(1:sig_startj+nj)-mean(sgmj.data(1:sig_startj+nj))];
                                                        else
                                                            if sig_startj+nj < length(sgmj.data)
                                                                sigj = sig_startj: sig_startj+nj;  % length for seismogram
                                                            else
                                                                sigj=sig_startj:length(sgmj.data);
                                                            end
                                                        end
                                                        %--------------------------------------------------------------------------------------
                                                        % remove means
                                                        sig_sgmj = sgmj.data(sigj)-mean(sgmj.data(sigj));        % main seismogram
                                                        
                                                        % And now check for clipping:
                                                        % use routine from Daniel Trugman, March 2016
                                                        % Use default setting as best he has found to date.
                                                        %--------------------------------------------------------------------------------------
                                                        clear clip;
                                                        clip = CLIP_DETECT(sig_sgmj);
                                                        if clip == 0    
                                                            % Removing the mean from the main signal filtered with 3 different - ICJG                                                    
                                                            % f1sig_sgmj = f1sgmj(sigj)-mean(f1sgmj(sigj));        % main seismogram
                                                            % f2sig_sgmj = f2sgmj(sigj)-mean(f2sgmj(sigj));        % main seismogram
                                                            f3sig_sgmj = f3sgmj(sigj)-mean(f3sgmj(sigj));        % main seismogram
                                                            %   Then EGF
                                                            noi_startl = round(sig_startl - sampl*nsec*1.05);
                                                            nl  = fix(nsec*sampl);             % number of samples 
                                                            % Put in to deal with picks beyond length of seismogram when P missed
                                                            if noi_startl < length(sgml.data) && noi_startl+ nl > 5
                                                                if noi_startl < 1
                                                                    nzero = zeros(abs(noi_startl)+1,1)';
                                                                    noi_sgml = [nzero sgml.data(1:noi_startl+nl)'-mean(sgml.data(1:noi_startl+nl))];
                                                                else
                                                                    noil = noi_startl: noi_startl+nl;  % length for seismogram-noise
                                                                    noi_sgml = sgml.data(noil)-mean(sgml.data(noil));        % main seismogram - noise
                                                                end                                                               
                                                                if sig_startl+nl > length(sgml.data)  % Deal with short seismograms
                                                                    szero = zeros(sig_startl+nl - length(sgml.data),1); %columns
                                                                    sigl = sig_startl:length(sgml.data); %columns
                                                                    sig_sgml = [sgml.data(sigl)-mean(sgml.data(sigl));  szero];
                                                                    % f1sig_sgml = [f1sgml(sigl)-mean(f1sgml(sigl)); szero];        % EGF seismogram
                                                                    % f2sig_sgml = [f2sgml(sigl)-mean(f2sgml(sigl)); szero];        % EGF seismogram
                                                                    f3sig_sgml = [f3sgml(sigl)-mean(f3sgml(sigl)); szero];        % EGF seismogram
                                                                elseif sig_startl < 1
                                                                    szero = zeros(abs(sig_startl)+1,1)';
                                                                    sig_sgml = [szero sgml.data(1:sig_startl+nl)-mean(sgml.data(1:sig_startl+nl))];
                                                                else
                                                                    sigl = sig_startl: sig_startl+nl;  % length for seismogram
                                                                    sig_sgml = sgml.data(sigl)-mean(sgml.data(sigl));        % EGF seismogram
                                                                    % f1sig_sgml = f1sgml(sigl)-mean(f1sgml(sigl));        % EGF seismogram
                                                                    % f2sig_sgml = f2sgml(sigl)-mean(f2sgml(sigl));        % EGF seismogram
                                                                    f3sig_sgml = f3sgml(sigl)-mean(f3sgml(sigl));        % EGF seismogram
                                                                end      % End of if seismogram too short.
                                                                %--------------------------------------------------------------------------------------                                                             
                                                                % Now need to interpolate to sort out sample rate if different.
                                                                % CHECK SAMPLING RATES ARE EQUAL
                                                                min_samp = min([sampj sampl]);
                                                                max_samp = max([sampj sampl]);
                                                                samprat = max_samp/min_samp;                                                              
                                                                if sampl ~= sampj
                                                                    ['sample rates not equal - Any Cross correlation stuff will not work....']
                                                                    pause
                                                                    if sampj == min_samp
                                                                        sig_sgmj = interp(sig_sgmj,samprat);
                                                                        noi_sgmj = interp(noi_sgmj,samprat);
                                                                    elseif sampl == min_samp
                                                                        sig_sgml = interp(sig_sgml,samprat);
                                                                        noi_sgml = interp(noi_sgml,samprat);
                                                                    end     % end of interpolation
                                                                end   % end of if sample rates not equal                                                               
                                                                if length(sig_sgmj) > length(sig_sgml)
                                                                    sig_sgmj=sig_sgmj(1:length(sig_sgml))
                                                                else if length(sig_sgml) > length(sig_sgmj)
                                                                        sig_sgml=sig_sgmj(1:length(sig_sgmj))
                                                                    end
                                                                end                                                               
                                                                if length(sig_sgmj) ~=length(sig_sgml)
                                                                    ['seismogram lengths not equal']
                                                                    length(sig_sgmj),length(sig_sgml)
                                                                    keyboard
                                                                end                                                           
                                                                %--------------------------------------------------------------------------------------
                                                                % FUN BEGINS
                                                                % CORRELATION IS ABOUT TO START
                                                                % I can do this in a function (MAYBE) - ICJG
                                                                % Perform cross correlation to raw signal
                                                                [xcora xcorb] = xcorr(sig_sgmj,sig_sgml,'coeff');
                                                                [mxcora mxcorab] = max(xcora);nxcor = length(xcora);
                                                                xcor_offset = mxcorab - ((nxcor -1)/2 +1);                                                             
                                                                % [xcoraf1 xcorbf1] = xcorr(f1sig_sgmj,f1sig_sgml,'coeff');
                                                                % [mxcoraf1 mxcorabf1] = max(xcoraf1); nxcorf1 = length(xcoraf1);
                                                                % xcorf1_offset = mxcorabf1 - ((nxcorf1 -1)/2 +1);                                                               
                                                                % [xcoraf2 xcorbf2] = xcorr(f2sig_sgmj,f2sig_sgml,'coeff');
                                                                % [mxcoraf2 mxcorabf2] = max(xcoraf2); nxcorf2 = length(xcoraf2);
                                                                % xcorf2_offset = mxcorabf2 - ((nxcorf2 -1)/2 +1);  
                                                                % Correlation for main signal filter3 - THE ONE IM USING - ICJG                                                             
                                                                [xcoraf3 xcorbf3] = xcorr(f3sig_sgmj,f3sig_sgml,'coeff');
                                                                [mxcoraf3 mxcorabf3] = max(xcoraf3); nxcorf3 = length(xcoraf3);
                                                                xcorf3_offset = mxcorabf3 - ((nxcorf3 -1)/2 +1);
                                                                %--------------------------------------------------------------------------------------                    
                                                                % get name for output variables...
                                                                clear output_name
                                                                % Need to make from end as station variable length..
                                                                output_name = ['SD_' nsgmj(18:end) '_' nsgmj(3:16)  ];
                                                                output_namel= nsgml(3:16);
                                                                output_namell = ['SD_'  nsgml(18:end) '_' nsgml(3:16)];
                                                                
                                                                eval([output_name '_' output_namel '_' PS_letter(PS) '.xcorr' int2str(kk) '(:,1) = xcora;']);
                                                                eval([output_name '_' output_namel '_' PS_letter(PS) '.xcorr' int2str(kk) '(:,2) = xcorb;']);                                                                
                                                                % eval([output_name '_' output_namel '_' PS_letter(PS) '.xcorr' int2str(kk) 'f1(:,1) = xcoraf1;']);
                                                                % eval([output_name '_' output_namel '_' PS_letter(PS) '.xcorr' int2str(kk) 'f1(:,2) = xcorbf1;']);                                                              
                                                                % eval([output_name '_' output_namel '_' PS_letter(PS) '.xcorr' int2str(kk) 'f2(:,1) = xcoraf2;']);
                                                                % eval([output_name '_' output_namel '_' PS_letter(PS) '.xcorr' int2str(kk) 'f2(:,2) = xcorbf2;']);                                                               
                                                                eval([output_name '_' output_namel '_' PS_letter(PS) '.xcorr' int2str(kk) 'f3(:,1) = xcoraf3;']);
                                                                eval([output_name '_' output_namel '_' PS_letter(PS) '.xcorr' int2str(kk) 'f3(:,2) = xcorbf3;']); 

                                                                eval([output_name '_' output_namel '_' PS_letter(PS) '.xcor_offset' int2str(kk) '(:,1) = xcor_offset;']);
                                                                % eval([output_name '_' output_namel '_' PS_letter(PS) '.xcor_offset' int2str(kk) 'f1(:,1) = xcorf1_offset;']);
                                                                % eval([output_name '_' output_namel '_' PS_letter(PS) '.xcor_offset' int2str(kk) 'f2(:,1) = xcorf2_offset;']);
                                                                eval([output_name '_' output_namel '_' PS_letter(PS) '.xcor_offset' int2str(kk) 'f3(:,1) = xcorf3_offset;']);
                                                                %--------------------------------------------------------------------------------------
                                                                %% Add limit to cross correlation for which a deconvolution is
                                                                %performed - to save time!!!
                                                                % ie. use the crazy one (0.99) on the raw data and the xcor_limit
                                                                % on the lowest f option.
                                                                % comment out this loop so actually see what is happening.......
                                                                if (max(xcoraf3) >= xcor_limit) && (max(xcoraf3) < 0.99) 
                                                                    % Set length etc. for deconvolution. ** Done after interpolation.
                                                                    %     Problem - if different samp, then signals are different lengths.
                                                                    %     so keep the shortest length for everything.
                                                                    n = min([length(sig_sgmj) length(sig_sgml)]);
                                                                    %     Alert if there is a big difference:
                                                                    if abs(length(sig_sgmj) - length(sig_sgml)) > 10
                                                                        length(sig_sgmj) - length(sig_sgml)
                                                                    end
                                                                    p  = nextpow2(n);
                                                                    n2 = 2^p;
                                                                    % Set lowest original Nyquist frequency for low passing STF.
                                                                    Nyf_samp = min([sampj sampl]);   % Find minimum samples of data
                                                                    lp = (80/100)*Nyf_samp/2;
                                                                    %--------------------------------------------------------------------------------------                                                      
                                                                    % Perform deconvolution:
                                                                    % Call to new matlab scripts....
                                                                    % uses interpolation sample rate so does not filter STF...
                                                                    clear tfun spec_ratio spec1 spec2 freq                                                                 
                                                                    [tfun,spec_ratio,spec1,spec2,freq] = mt_deconv(sig_sgmj,sig_sgml,noi_sgmj,noi_sgml,1./Nyf_samp);                                                                   
                                                                    % remove nonsense from missing data....
                                                                    if sum(sum(isnan(tfun))) == 0                                                                  
                                                                        % rename results for saving.
                                                                        eval([output_name ' = sgmj;'])
                                                                        eval([output_namell ' = sgml;'])
                                                                        if exist(['mname_' cluster.name(1,:) '_data.mat']) == 2
                                                                            % save data to data file
                                                                            eval(['save  mname_' cluster.name(1,:) '_data ' output_name ' ' output_namell ' -append'])
                                                                        else
                                                                            eval(['save  mname_' cluster.name(1,:) '_data ' output_name ' ' output_namell ])%' -append'])                                                                     
                                                                        end
                                                                        %%diary off
                                                                        eval(['clear ' output_name ' ' output_namell])%% clear out
                                                                        %--------------------------------------------------------------------------------------
                                                                        % Set the first spec in this to be spec_ratio1. Then save spec_ratio2 but do not use or do any fitting with. So should be able to work with all subsequent scripts with no change....
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.fratio = spec_ratio(:,1);'])                                                         
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.samp = Nyf_samp;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.f = freq;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.stf = tfun;']);
                                                                        eval([output_name  '_MT.spec' PS_letter(PS) '(:,1) = freq;']);
                                                                        eval([output_name  '_MT.spec' PS_letter(PS) '(:,2) = spec1(:,1);']);
                                                                        eval([output_name  '_MT.samp = sampj;']);
                                                                        eval([output_name  '_MT.spec' PS_letter(PS) '_noise(:,1) = freq;']);
                                                                        eval([output_name  '_MT.spec' PS_letter(PS) '_noise(:,2) = spec1(:,2);']);
                                                                        eval([output_name(1:end-14) output_namel '_MT.spec' PS_letter(PS) '(:,1) = freq;']);
                                                                        eval([output_name(1:end-14) output_namel '_MT.spec' PS_letter(PS) '(:,2) = spec2(:,1);']);
                                                                        eval([output_name(1:end-14) output_namel '_MT.samp = sampl;']);
                                                                        eval([output_name(1:end-14) output_namel  '_MT.spec' PS_letter(PS) '_noise(:,1) = freq;']);
                                                                        eval([output_name(1:end-14) output_namel  '_MT.spec' PS_letter(PS) '_noise(:,2) = spec2(:,2);']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.epi_d = epi_d;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.del = sgmj.del;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.hyp_sep = hyp_sep;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.ml = [sgmj.ml sgml.ml];']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.mw = [sgmj.mw sgml.mw];']);
                                                                        % elat,elon= event
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.elat = [sgmj.lat sgml.lat];']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.elon = [sgmj.lon sgml.lon];']);
                                                                        % slat,slon= station
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.slat = sgmj.slat;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.slon = sgmj.slon;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.pazi = sgmj.pazi;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.main_index = k;']);
                                                                        eval([output_name '_' output_namel '_' PS_letter(PS) '.egf_index = l_e;']);
                                                                        %--------------------------------------------------------------------------------------
                                                                        % MAKING FIGURES - ICJG
                                                                        % Need to add this into a function 
                                                                        % Make a figure to keep...
                                                                        % figure(3);
                                                                        % if kk==1
                                                                        %     clf;
                                                                        % end
                                                                        % subplot(6,2,kk); plot(sig_sgmj);                                                      
                                                                        % xlim([0 length(sig_sgmj)]);
                                                                        % title1a=['Southern DWARFS - ' output_name(4:8) '-' output_name(10) '-' output_name(12:end) '-' output_namel ' ' PS_letter(PS)];
                                                                        % if PS == 1
                                                                        %     title1b=['Start A  = ' int2str(sig_startj) ', Start B = ' int2str(sig_startl) ', window = ' num2str(nsec) 's'];
                                                                        % else
                                                                        %     title1b=['Start A  = ' int2str(sig_startj) ', Start B = ' int2str(sig_startl) ', window = ' num2str(nsec) 's'];
                                                                        % end
                                                                        % title({title1a;'Main event - raw signal'},'FontSize', 7);
                                                                        % ylabel(['Stn Dist: ' num2str(sgmj.del,'%5.1f')],'FontSize', 6)
                                                                        
                                                                        % subplot(6,2,kk+2); plot(sig_sgml,'r');
                                                                        % xlim([0 length(sig_sgml)]);
                                                                        % title({title1b;'EGF - raw signal'},'FontSize', 7);
                                                                        
                                                                        % subplot(6,2,kk+4);plot(f3sig_sgmj./range(f3sig_sgmj),'b');hold on; plot(f3sig_sgml./range(f3sig_sgml),'r');
                                                                        % xlim([0 length(sig_sgmj)]);
                                                                        % xlabel(['Samps, Filter F3 = ' num2str(lowpassf3,3) ' Hz, Offset f3 ' int2str(xcorf3_offset) ' samples'],'FontSize', 7);
                                                                        % title(['Filtered signals'],'FontSize', 7)
                                                                        
                                                                        % % And cross correlation figure.
                                                                        % subplot(6,2,kk+6);plot(xcorbf3,xcoraf3);
                                                                        % title('Correlation function','FontSize', 7);
                                                                        % xlabel(['Max correlation value : raw signal = ' num2str(mxcora,'%5.3f') ' Filtered = ' num2str(mxcoraf3,'%5.3f')], 'FontSize', 7)
                                                                        
                                                                        % subplot(3,4,7+2*kk); loglog(freq,spec_ratio(:,1),'b',freq,spec_ratio(:,2),'r');
                                                                        % axis([1/nsec min(Nyf_samp)/2 max(spec_ratio(:,1))*1.1/1e3 max(spec_ratio(:,1))*1.1]);
                                                                        % title1d=['Spectral Ratio: Mw' num2str(cluster.mag(1)) '/' num2str(cluster.mag(l)) ' ML' num2str(sgmj.ml) '/' num2str(sgml.ml)];
                                                                        % title(title1d,'FontSize', 6);ylabel('Amplitude Ratio'); %: German .m code');
                                                                        % xlabel(['Frequency: ' num2str(1/nsec) ' - '  num2str(min(Nyf_samp)/2) ' Hz'],'FontSize', 7);
                                                                        
                                                                        % subplot(3,4,8+2*kk); plot(tfun);
                                                                        % axis([1 length(tfun) min(tfun)*1.1 max(tfun)*1.1])
                                                                        % title1e=['Epi sep ' num2str(epi_d,'%4.2f') ' km '];
                                                                        % title(title1e,'FontSize', 7);xlabel('STF: Samples','FontSize', 7);
                                                                        %--------------------------------------------------------------------------------------
                                                                        % Now make a matrix of the cross correlation results
                                                                        if kk==1
                                                                            xcor_out0(count_xcorr0,:) = [j l k PS mxcora xcor_offset];
                                                                            count_xcorr0 = count_xcorr0+1;
                                                                        elseif kk==2
                                                                            xcor_outx(count_xcorrx,:) = [j l k PS mxcora xcor_offset];
                                                                            count_xcorrx = count_xcorrx+1;
                                                                        end                                                                     
                                                                        clear stf sig_startj sig_startl stf_MT stf_MT_noij stf_MT_noil sig_sgmj sig_sgml noi_sgmj noi_sgml sampi lp 
                                                                        %figure(3); orient landscape
                                                                        if kk == 1
                                                                            x=-xcorf3_offset;y=-xcorf3_offset;
                                                                        else
                                                                            %eval(['print -dpdf -bestfit ' output_name  '_' output_namel '_' PS_letter(PS) '_nsec' num2str(nsec) '_1_egfMT_mG.pdf'])
                                                                            x=0;y=0;
                                                                        end                                                                 
                                                                    end % removing nonsense when get no real data.... (sum(isnan))
                                                                end % end of if xcor cross correlation limit                                                               
                                                            end   % end of noise starts before end of record  l
                                                        end  % end if clip = 0
                                                    end   % end of noise starts before end of record  j
                                                end  % End of if loop for having sig_start
                                            end  % end for kk repeat loop to do with and without xcorrelation
                                        end  % end of PS =1,2
                                        %  end  % end of both sgmj.p & PS=1 OR sgmj.s and PS=2
                                    end  % end of if sampj > 20
                                    %%  end  % end for if loop having a structure l exist_l
                                    % delete ones for which nothing worked...
                                    delete_pair = 0;
                                    for PS = 1:2
                                        if exist('output_name') == 1
                                            save_ratio = [output_name '_' output_namel '_' PS_letter(PS)];
                                            if exist(save_ratio) == 1
                                                if isfield(eval(save_ratio),'fratio') == 1
                                                    delete_pair = 1;
                                                else
                                                    eval(['clear ' save_ratio])
                                                end  % end of isfield fratio
                                            end  % end of exist save_ratio
                                        end  % end of exist output_name
                                    end  % end of PS
                                    if delete_pair == 0
                                        eval(['clear ' nsgml ' ' nsgml '_MT ' ])
                                    end
                                    break   % end of if get a l to j match
                                end     % end of if get a l to j match
                            end  % end for l_e loop over stn and chan in EGF match k main smgm
                            % keyboard
                        end  % end of if msta is NOT empty
                        
                    end % end of k loop for all seismograms in j
                end    % if a datafile for EGF event
                %% need to save data here for seismograms..
            end   % end of for loop over l for all EGF
            
            ['save egfs now'],
            save egf_temp_variables j k l  all_mname_list pcnt
            clear no* ml* max* mi* m  ic el* PS* y x ti* s* p* no* co*  mx* wl MT* n n2 nj nl nsg* nx*
            clear xcor_offset xcora xcorb recor* nz* mgsort*  ans azi clen* dist dummy e* gc*
            clear mname_*
            eval(['mname_' cluster.name(1,:) ' = cluster'])
            % Trying to sort out save so as to not crash when only P or only S
            % or neither
            cd (clusterj_folder)
            try
                eval(['save mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_1_egf.mat -regexp (P|S|MT)$'])
                eval(['save mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_1_egf.mat -append all* Ny* ca* h* m* n* v* x* pcnt'])
            catch
            end  % end of try saving....
            %eval(['mkdir pdf_plots/' cluster.name(1,:)])
            %eval(['!mv *' cluster.name(1,:) '*.pdf pdf_plots/' cluster.name(1,:)]);
            % diary off
            
            %   clear NZ* a* b* c* d* e* f* g* h* i* k* l* m* n* o* p* q* r* s* t* u* v* w* x* y* z*
            
            % Global clear to speed things up....
            clear
            % end  % end of ???
            if exist('egf_temp_variables.mat') == 2
                load egf_temp_variables
            end
        end   % end of if data file exists....
        
    end   % end of check whether have already done this event.....
    clearvars -except main_events
end  % end of j loop over main events
