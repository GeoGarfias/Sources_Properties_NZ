% Modified by ICJ-G 25/08/20 for SAMBA data
%============================================================================================================
%% Only uses existing ratios.....
%% calculates low and high f from spectral noise ratios
%% if needed using function: get_frange_function.m


%% July 2012: Adapted this so as to have SAME frange for ALL spectra.
%% Use NaN where outside signal to noise bandwidth
%% Need to define high1f and low1f as this range.

iplot = 1;   % subplot index
ifig = 1;    % figure index
   
% Inputs:
% Need lowf, and high f for log resampling - set as
% 3rd point, and 80% Nyquist of lowest sample rate event
%============================================================================================================
% Evaluate the SAMBA cluster
eval(['load /Users/home/juarezilma/Master/NETWORKS/SAMBA/EGF_matlab/WORK/cluster_SAMBA.mat']);

%for j=1:length(all_mname_list)
for j=9
  % Evaluate the SAMBA cluster again (?)
  eval(['load /Users/home/juarezilma/Master/NETWORKS/SAMBA/EGF_matlab/WORK/cluster_SAMBA.mat']);
  all_mname_list = eval(['mname_cat']);
  cluster = eval(['cat_' all_mname_list(j,:)]);

  figure(1);clf

  [mlen1 mlen2]=size(cluster.name)

  clear  nsec;
  nsec = nsec_value(cluster.mag(1))
  % Check not already done...
  go_ahead = 0;  % whether to do this event...

    if exist(['mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_2_egf_250Hz_SN9.mat']) == 0 % 2 because output from second code
        % if mat file exists - load it in:
        ['mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_1_egf.mat']
        if exist(['mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_1_egf.mat']) == 2
        % make sure do not reset j when load in data file.....
        clear j_current;
        j_current = j;
        eval(['load mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_1_egf.mat'])
        go_ahead = 1;
        end
        if go_ahead == 1
            j = j_current;
            all_mname_list(j,:)
            % reload cluster to be sure have right one.
            load /Users/home/juarezilma/Master/NETWORKS/SAMBA/EGF_matlab/WORK/cluster_SAMBA.mat
            % make sure nsec is the same value here
            % Define parameters that govern range of all spec so as to enable stacking.....
            max_f = 80;   % Define this here assumes not > 200 samp/s data...
            low1f = (1/nsec)
            high1f = max_f
            PS_letter = ['P';'S'];  
            for PS = 1:2
                % first get list of ratios in file ending in PS
                ratio_list_name = ['I_*' PS_letter(PS)]; % I for italy - written out at end of previous code. 
                ratio_list_cell = who(ratio_list_name);
                rnlen1 = length(ratio_list_cell)
                for l = 1:rnlen1
                    [num2str(l) ' /' num2str(rnlen1)]
                    %nratio = ratio_list(l,:)
                    nratio = ratio_list_cell{l}
                    ratio = eval(nratio);
                    xcor_limit = 0.5 % this the point to then cut out the rubbish, set as same as in previous code to keep using all the same stuff.
                    if (max(ratio.xcorr2f3(:,1)) >= xcor_limit) &&  (max(ratio.xcorr2f3(:,1))  < 0.99)
                        % these variables are needed for get_frange_function  
                        % as variable length station need to define from end...
                        nsgmj = [nratio(1:end-16) 'MT'] % see email for how to define these for your output.
                        nsgml = [nratio(1:end-31) nratio(end-15:end-1) 'MT']
                        if length(nsgmj) > 0 && length(nsgml) > 0
                            sgml=eval(nsgml);
                            sgmj=eval(nsgmj);
                            % Now need to set f range and log resample
                            % use 2 following lines if no frange set (see get_frange.m)
                            %% Set minimum and maximum value for fresamp so that have same values. Then pick which to use based on actual frange....
                            % check if have spectrum!
                            if isfield(sgmj,['spec' PS_letter(PS)]) == 1 %does it have a .specS/P ?
                                % First check if have frange...
                                if isfield(ratio,'frange') == 0
                                    FRANGE = get_frange_function2I1(nratio,ratio,nsgmj,sgmj,nsgml,sgml,PS_letter,PS)
                                    eval([nratio '.frange = [FRANGE(5) FRANGE(6)];']);
                                    eval([nsgml '.' PS_letter(PS) 'frange = [FRANGE(3) FRANGE(4)];']);
                                    eval([nsgmj '.' PS_letter(PS) 'frange = [FRANGE(1) FRANGE(2)];']);
                                    % remake these so have updated ones in work place. 
                                    ratio=eval(nratio)
                                    sgml=eval(nsgml)
                                    sgmj=eval(nsgmj)
                                else
                                end
                                highf = log10(ratio.f(ratio.frange(2)));
                                % deal with Nyquist frequency
                                Nyf_samp = min(ratio.samp)
                                if highf >= log10(0.8*Nyf_samp/2)
                                    highf = log10(0.8*Nyf_samp/2)
                                end
                                %% Deal with fact that first f point is 0....
                                if ratio.frange(1) < 3
                                    lowf = log10(ratio.f(3));
                                    % fix from Christine...
                                    eval([nratio '.frange(1) = 3;']);
                                else
                                    lowf = log10(ratio.f(ratio.frange(1)));
                                end
                                deltaf = 0.05;      % Not sure about deltaf?!
                                % Decide if worth even calculating resampled ratio.....     
                                if highf > lowf
                                    if highf > log10(high1f)
                                        ['High1f NOT high enough']
                                        keyboard
                                    end
                                    fratio_resamp = log_resample(ratio.fratio,ratio.f,log10(low1f),log10(high1f),deltaf);
                                    % Now need to set part that is outside lowf and highf to NaN.
                                    for fi = 1:length(fratio_resamp)
                                        if (fratio_resamp(fi,1)) < 10.^lowf || (fratio_resamp(fi,1) > 10.^highf)
                                            fratio_resamp(fi,2) = NaN;
                                        end
                                    end
                                    % Save it:
                                    eval([nratio '.f_resamp = fratio_resamp(:,1);']);
                                    eval([nratio '.fratio_resamp = fratio_resamp(:,2);']);
                                    % Re evaluate so as to get new fields
                                    ratio = eval(nratio);
                                end   % end of loop over if highf>lowf
                            else  % If do not meet xcor_limit, so no ratio saved, clear out.....
                                eval(['clear ' nratio])
                            end % if spec PS exists within sgmj or sgml
                        end  % end of if nsgml and nsgmj exist
                    end   % end of if loop for xcora >= xcor_limit
                end  % end for l loop over existing ratios
            end        % End of loop over P and S
            ['save egfs now']
            %save egf_temp_variables j l PS xarea all_mname_list narea
            save egf_temp_variables j l PS all_mname_list
            clear no* ml* ma* mi* m l kk k j ic el* PS* y x ti* s* p* no* co* mx* wl MT* n n2 nj nl nsg* nx*
            clear xcor_offset xcora xcorb recor* nz* mgsort* l* ans azi clen* dist dummy e* gc* jj
            clear mname_*
            eval(['mname_' cluster.name(1,:) ' = cluster'])
            % Check if actually have anything to save...... 
            ratio_list_name_out = ['I_*_*_*_*_*'];
            ratio_list_out_cell = who(ratio_list_name_out);
            if length(ratio_list_out_cell) > 0
                eval(['save mname_' cluster.name(1,:) '_nsec' num2str(nsec) '_2_egf_SN9.mat'])
            end 
            % Global clear to speed things up
            clear
        end   % end of if input file exists - now go_ahead =1
        %end  % end of if cluster exists  
        load egf_temp_variables
    end  % end of if done before
end  % end of j loop over main events




