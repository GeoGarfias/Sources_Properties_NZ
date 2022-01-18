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

Once you have the cluster file, you need to get individual station component SAC files for all your main-events. (I’m not sure why, but I decided to get NORDIC files first and then SAC files, the procedure to get the SAC files is explained ahead)

2. Getting SAC files.
   (I’m not sure why, but I decided to get NORDIC files first and then SAC files, that procedure is explained here)
   2.1. Run the python code 'Getnordic_EGF.py'




	a) Neet to run the python code call ‘Getnordic_EGF.py’.
	    This code was made to get data from VUW servers. 
    • Edits: Paths for .XML and  eqinfo.mat' files and id_list variable. 
              The id_list variable is the path of a directory that includes .CSV files for each MAIN event. The .CSV files only contain the *ID of the EGFs event asocciated to a MAIN event. 
		       		*You can find the ID of the events in the cluster_NETWORKNAME.mat file. 			              There is a  code name ‘get_idfiles.m’ that get the ID of the EGF events 					associated with a MAIN event and write the .CVS files for each MAIN event.

				You also need to edit the paths  of where the nordic files are  			           going to be saved in the def_nordic function are going to be 					saved.  There should be a directory for each MAIN event and 				subdirectories within them called ‘S_files’ and ‘Wav_files’. The 				output of the code is storage in those subdirectories. There 					is .csh code call ‘make_dir.csh’ that will make the directories for 				you.

    • Output: SEISAN files, s-files and wav-files for all MAIN and EGF events.		
	After getting all the SEISAN files, it is time to get the SAC files.

	b) Run the python code call ‘cat_to_EGF_sacfiles.py’
	This code will generate SAC traces per component per each single event. 
    • Edits: Only need to edit the paths. 
					
    • Output: Individual component stations SAC traces. 
      (Networkcode_station_component_eventname.sac_cut - AF_EOROX_Z_ev20180512145612.sac_cut)



3.
4. 
5. dfdf
6. Run EGF_decon.m:  If the EGF_decon.m presents a error saying that there are no picks, run the del_traces_nopicks.m code. this will delete the traces that not contains picks. We need P and S picks for all the traces.  


3. Hola
4. hola
5. 
NOTES:
1. All codes and file names are between quotes (' ') - example: 'make_egf_clusters.m' and 'eqinfo.mat'
2. File extension names are written with capital letters only to highlight them. Extension should be written in lower case for files. 
