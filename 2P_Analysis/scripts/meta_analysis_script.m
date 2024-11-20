% 1) Pre-process and analyze spine 2p data to obtain FOV's zscored data + OSI/Ori pref
% 2) collect all data into a structure array for each visual stim
% 3) identify visual responsive spines based on criteria
% 4) output all data into a table where each row is a spine, and the columns indicate different parameters measured for the spine for specific visual stim (binoc, contra, or ipsi)
% 5) output chronically tracked spines into two tables: spines tracked from D1 to D5 and spines tracked from D5 to D10

%% Set path
clear all
path = '/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/2P_Analysis';
addpath(genpath(path))
distance_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';

%% 1) run analysis spine 2p data function
disp('Running analysis spine 2p data function')
mouse_files = {'Mouse'};
input_path = fullfile(path, 'demo_data', '2P_data')
analysis_spine_2p_data_function(input_path, mouse_files)

%% 2) run analyze processed data function
disp('Running analysis processed data function')
mouse_files = {'Mouse'};
input_path = '/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/2P_Analysis/output_files';
savefile = 'all_analyzed_spines.mat';
analyzed_processed_data_function(path, savefile, mouse_files)

%% 3) run identify visual responses
disp('Running identify vis responses function')
analyzefile =  'all_analyzed_spines.mat';
savefile = 'visually_responsive_spines.mat';
mean_amp_thresh = 0.5;
p_val_thresh = 0.05;
identify_vis_responses(path,analyzefile, savefile, mean_amp_thresh,p_val_thresh);

%% 4) run create table for all 3 sessions per spines
disp('Running create table for all 3 sessions')
analyzefile = 'all_analyzed_spines.mat';
visfile = 'visually_responsive_spines.mat';
savefile = 'data_table';
savefile_mat = [savefile, '.mat'];
savefile_csv = [savefile, '.csv'];
savepath = fullfile(path, 'Analyzed Data');
load(fullfile(savepath,analyzefile));
load(fullfile(savepath,visfile));
create_table_of_all3sessions_function(savepath, distance_path, all_stims, all_vis_stims, savefile_mat, savefile_csv);

%% 5) run create table tracked spines
disp('Running create table tracked spines')
savefile = 'tracked_spines.mat';
spinefile = 'spine_data_table.mat';
create_table_tracked_spines_function(distance_path, path, spinefile, savefile)
