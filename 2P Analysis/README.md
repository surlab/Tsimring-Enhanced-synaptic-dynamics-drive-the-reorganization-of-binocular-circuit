# 2P_analysis

Analyzes the time-series flourescence data collected during two photon imaging with visual stimulate and compiles analyzed properties for each spine into several .mat files
## Installation Instructions
Analysis run on MATLAB2021b, no specific software installation required besides downloading MATLAB (any version should suffice)
## Usage Instructions

### Overview
All code can be run through "scripts/meta_analysis_script.m", which contains 5 sections of code:

#### Setting directories/paths
Before running the code, be sure to change the following paths to your home directory and location of code: 
- homepath = '/Users/ktsimring/Documents/';
- path = '/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/2P_Analysis';


#### Section 1
 first section runs the code for "analysis_spine_2p_data_function.m", where for each FOV and viewing condition, it extracts and z-scores the flourescent values from each ROI (spine or soma), aligns the timestamps to the visual stimulus, calculates the meam and standard dev of traces aligned to the unique stimuli, and analyzes the ROI's tuning properties
##### Input file: data is in "demo_data/2P_data" folder 
1. ROIdata.mat (time series flourscent data from each ROI): Tx(M*3) array, where T represents timestamps and M represents spine ROI (for each spine, there is a corresponding dendritic rectangular ROI and background ROI drawn)
2. stim.mat (stimulus file): TxNx3 array, where T represents number of trials (10), N represents number of directions (8), the third dimension represents spatial frequencies testing (usually 1), and the fourth dimension is the brightness which is the same across stimuli
3. DendX*_VoltageRecording.csv (voltage recording to align timestamps: Tx2 array where first column is the timestamp, and second column is the voltage sent from the data acquisition device 
##### Output files: For each FOV, day, and viewing condition there will be two files: 
1. normalized_data_by_stim.mat (analyzed flourescent properties) 
2. mean_amp_ori_analysis.mat (orientation tuning properties for FOV)
##### Runtime: ~ 1 minute on all data in "demo_data/2P_data" folder

#### Section 2
Second section runs the code for "analyzed_processed_data_function.m", which compiles the analyzed properties for each ROI across FOVs and viewing conditions into a 1x3 array of structures
##### Input file: for each mouse, day, neuron, FOV, you need:
1. normalized_data_by_stim.mat
2. mean_amp_ori_analysis.mat
##### Output file: all_analyzed_spines.mat
- each field represents an array of all the spine's specific properties (i.e mean amplitude, OSI, DSI, etc), and each row is the properties for that specific viewing condition (binoc, contra, or ipsi)
##### Runtime: ~ 0.5 seconds on all data in "demo_data/2P_data" folder

#### Section 3
Third section runs the code for "identify_vis_responses.m" which identifies which spines are visually responsive for each visual stimulation condition
##### Input file: all_analyzed_spines.mat
##### Output file: visually_responsive_spines.mat
##### Runtime: ~ 0.05 seconds on all data in "demo_data/2P_data" folder

#### Section 4
Fourth sections runs the code for "create_table_of_all3sessions_function.m" which compiles all the previously analyzed properties of the spines, and runs additional analysis for each spine (i.e correlation to nearest neighbor)
##### Input file: 
1. all_analyzed_spines.mat
2. visually_responsive_spines.mat
3. dendritic_distance.csv: dendritic distance between spines for each FOV and viewing condition (data in "demo_data/2P_data")
4. spine_stats.csv: structural stats for each spine  (data in "demo_data/2P_data")
5. RoiSet.zip: roi identities (data in "demo_data/2P_data")
##### Output files:
1. spine_data_table.mat: N*M array, where N represents total number of spines recorded across viewing sessions, dendritic FOVs, neurons, and mice and M represents all the analyzed properties for each spine
2. soma_data_table.mat: N*M array, where N represents total number of somas recorded across viewing sessions, neurons, and mice and M represents all the analyzed properties for each soma
##### Runtime: ~ 7 seconds on all data in "demo_data/2P_data" folder

#### Section 5
 Fifth section runs the code for "create_table_tracked_spines_function" for spines that were tracked from D1 to D5 or D5 to D10 to combine analyzed properties from the tracked spine from the first and second timepoint into a single data table
##### Input file: 
1. spine_data_table.mat
2. variable "modify_mat" from alignment_based_on_fiducial.mat, which is a Nx2 array, providi the ROI identity of the tracked spine on timepoint 1 (first column) and 2 (second column). Spines that were missing on 1st or 2nd timepoint are labeled as 0, and spines that had unknown identity were labeled as -1
##### Output file: 
tracked_spines.mat, which has two variables that are N*2 arrays where N is the total number of spines, and the first column is data from D1 (or D5) and second column is data from D5 (or D10)
##### Runtime: ~ 0.2 seconds on all data in "demo_data/2P_data" folder
      
        
       
  
