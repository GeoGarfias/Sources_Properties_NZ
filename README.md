# Sources_Properties_NZ
==================================================
For questions, please email at ilmadelcarmen.juarezgarfias@vuw.ac.nz (academic email) or at icjgarfias@gmail.com (personal email)
=================================================================================================================================================================

In this repository you will find codes to calculate corner frequency, stress drop and directivity following Abercrombie, R. E., Bannister, S., Ristau, J., and Doser, D. 2017 (https://doi.org/10.1093/gji/ggw393) and Abercrombie, R. E., Poli, P., and Bannister, S. 2017 (https://doi.org/10.1002/2017JB014935).

This repository has both MATLAB and Python codes. The first stage of the process is mainly Python based and all the methodology is made in MATLAB.

1. The first code you need to run is called 'make_egf_clusters.m'
   This code groups together earthquakes based on epicental distance and magnitude difference. You can modify the epicentral and magnitude criteria as preference.    This code provides a file called 'cluster_NETWORK_NAME.mat'. This file includes a list of earthquakes (with all the info about them, i.e. location, time,          magnitude, depth) and within each earthquake a list of possible associated (co-located) earthquakes.

   To run this code, you will need a catalogue of earthquakes in a .MAT (MATLAB) file. Nowadays, catalogues are in .XML files, so if you have a .XML catalogue        file, the code named 'cat_to_EGF_eqinfo.py' will provide you with a .CSV catalogue file and then you will have to import the .CSV catalogue file in MATLAB to        create a .MAT file. 
   
   If you have your catalogue in a .XML file follow the next steps before running the 'make_egf_clusters.m' MATLAB code. If you do not have your catalogue in a        .XML file, there is a file example in this repository ('eqinfo_ND.mat') that will give you an idea of the file you need to have to run the         'make_egf_clusters.m' MATLAB code.
   
   1.1. Before running the 'make_egf_clusters.m' MATLAB code, run the 'cat_to_EGF_eqinfo.py'
        This code will generate a .CSV file called 'eqinfo.csv'
   1.2. Import the 'eqinfo.csv' file into MATLAB and make a new 'eqinfo.mat' file. 
   1.3. Finally run the 'make_egf_clusters.m' MATLAB code - This code may take a while to run, everything depends of how many earthquakes are in your catalogue.
   1.4. After running the 'make_egf_clusters.m' MATLAB code do NOT clear the workspace, you need to manually save the next variables in a .MAT file called            'cluster_NETWORKNAME.mat':
               - All cat_ev* structures &
               - all_mname_list

Once you have the cluster file, you need to get individual station component SAC files for all the earthquakes. 
There are three examples files in this repository:
	1. XX_THOMX_Z_ev202001260401.sac_cut     2. XX_THOMX_N_ev202001260401.sac_cut     3. XX_THOMX_E_ev202001260401.sac_cut

2. Run 'load_sac.m' in MATLAB
   Once you had all the sac files, it's time for the data procesing in MATLAB.
   The first code will read the SAC data and transformed to the specific format the MATLAB procesing need.
   
   2.1. Before running, need to edit some lines in the code.
   	Line 9-10 & 16-17: Edit the path for the cluster and eqinfo files.
	Line 13: Variable m represents the event the code is reading. You can choose between a list of events or and individual event. 
	Line 19: Change the path of folder you will save the data.
	Line 43: Write a general 'key' that all your file names share.
  2.2. The output of this codes is a file called 'SD_j20190501033509_data.mat' (where SD = Network code). This file contains matrices that corresponds to every    	  individual SAC file. Those files are named as XX_SSSS_C_evYYYYMMDDHHMMSS_PROJECT (where XX = Network code SSSS = station and C = channel, 		            ND_GLORX_E_ev20180512145612_ND). 
  
3. Run the ‘convert_data.m’ code
   This converts the data structures into the event files with al the data of EGFs combined for the main-event. 
		
  3.1. The edits you need to make on the code are:
  	Line 4-5 & 9-10: Edit the path for the cluster and eqinfo files.
	Line 13: Change the path of folder you will save the data.
	Line 29-30: Make sure that the station and channel names  are being read correctly. 
	Line 8: Variable i represents the event the code is reading. You can choose between a list of events or and individual event. 
	
  3.2.The outputs of this codes are MAT file for each EGF event. This MAT files includes all the SAC information per event (picks, time, magnitude, etc). 

4. Run the 'EGF_decon.m' code
   This reads in the cluster file and then the data file written out in the previous step  and performs a deconvolution between main-events and EGFs-events. It’s    quite slow! Tries all pairings of data and only accepts the ones with cross-correlation value bigger than 0.5.
   
   4.1. The edits you need to make on the code are:
   	Lines 60,64: Edit the path for the cluster.
	Line 62-63: Variable j represents the event the code is reading. You can choose between a list of events or and           individual event. 
	Line 69: Change the path of folder you will save the data.
	Line 88: Choose a cross-correlation limit - 0.5 is suggested.
	Line 90-91: Set an approximate value for vp and vs.
	
	



   6.2
   7.3
   8.4
   9.
        a. Edits: Need to change project name throughout (find and replace FIORD=YOURPROJECT NAME). Also rename and save script in codes/project. 
        b. Nsec (see lines 149-157 and lines 198-204) is a key variable that defines the window length as a function of magnitude. You might have to play around with this to check it’s right, as it uses the given local magnitude (ML) to then convert to Mw and Mo to estimate the event duration based on empirical scaling. It is currently set to use the NZ scaling relationship (From Ristau 2014?). It should not be longer than the S-P time at the closest station. Try it with the current settings (which I set for Fiordland) and check the output. It may need shortening, because local stations will be closer for your study probably. 
        c. f1 and f3 define the bandpass filters for the deconvolution. This currently set to be magnitude dependent and calculated using nsec (around lines 207-229). You can check if this is appropriate in the plots made later in stack_spec_XXX below. If you can see long period signals are present on sites with low signal to noise, then you might need to play around with this.  Having a well defined magnitude scale is likely to help this. 
        d. Variable naming: main event is always a loop over j and is sgmj etc. Currently set to run for 1 event only (e.g. j-=4). Can loop over j to run for multiple mainshocks (e.g. line 77) .  EGF is always loop over l and is sgml etc. k usually loops over the stations and channels. The kk option is to do the deconvolution with the input picks. 
        e. If it gives output up to nsec=???? and then looks like it has finished quickly, it’s because the files have already been generated. The code includes a time-saving check to see if the expected output files (/work/PROJECT/mname_evYYYYMMDDHHMMSS_nsec_1_egf.mat) already exist. If it does, it will break. If you have changed something in the code, and want to re-run, the you’ll have to delete/rename the original output files to rerun. Or remove that checking flag. 
        f. Outputs: (/work/PROJECT/mname_evYYYYMMDDHHMMSS_nsec_1_egf.mat




6. dfdf
7. Run EGF_decon.m:  If the EGF_decon.m presents a error saying that there are no picks, run the del_traces_nopicks.m code. this will delete the traces that not contains picks. We need P and S picks for all the traces.  


3. Hola
4. hola
5. 
NOTES:
1. All codes and file names are between quotes (' ') - example: 'make_egf_clusters.m' and 'eqinfo.mat'
2. File extension names are written with capital letters only to highlight them. Extension should be written in lower case for files. 
