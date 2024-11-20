% 1) Pre-process and analyze spine 2p data to obtain FOV's zscored data + OSI/Ori pref
% 2) collect all data into a structure array for each visual stim
% 3) identify visual responsive spines based on criteria
% 4) output all data into a table where each row is a spine, and the columns indicate different parameters measured for the spine for specific visual stim (binoc, contra, or ipsi)
% 5) output chronically tracked spines into two tables: spines tracked from D1 to D5 and spines tracked from D5 to D10

%% Set path
clear all
homepath = '/Users/ktsimring/Documents/'; % change this path
path = '/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/2P_Analysis';
two_photon_data_folder = '2P_data'; %folder storing time-series data
FOV_alignment_data_folder = 'FOV_aligned_data'; %folder storing tracked spine identity data

output_folder = '2P_data_output'; %output folder (change for testing)

addpath(genpath(fullfile(homepath,path)))

%% 1) run analysis spine 2p data function
disp('Running analysis spine 2p data function')
mouse_files = {'Mouse'};
input_path = fullfile(homepath, path, 'demo_data', two_photon_data_folder);
tic
analysis_spine_2p_data_function(input_path, mouse_files)
toc
%% 2) run analyze processed data function
disp('Running analysis processed data function')
mouse_files = {'Mouse'};
input_path = fullfile(homepath, path, 'demo_data', two_photon_data_folder);
output_path = fullfile(homepath, path, 'demo_data', output_folder);
filename =  'all_analyzed_spines.mat';
tic
analyzed_processed_data_function(input_path, output_path, filename, mouse_files)
toc
%% 3) run identify visual responses
disp('Running identify vis responses function')
analyzefile =  'all_analyzed_spines.mat';
savefile = 'visually_responsive_spines.mat';
mean_amp_thresh = 0.5;
p_val_thresh = 0.05;
input_path = fullfile(homepath, path, 'demo_data', output_folder);
tic
identify_vis_responses(input_path,analyzefile, savefile, mean_amp_thresh,p_val_thresh);
toc
%% 4) run create table for all 3 sessions per spines
disp('Running create table for all 3 sessions')
analyzefile = 'all_analyzed_spines.mat';
visfile = 'visually_responsive_spines.mat';
savefile = 'data_table.mat';

input_path = fullfile(homepath, path, 'demo_data',  two_photon_data_folder);
savepath = fullfile(homepath, path, 'demo_data', output_folder);

load(fullfile(savepath,analyzefile));
load(fullfile(savepath,visfile));
tic
create_table_of_all3sessions_function(savepath, input_path, all_stims, all_vis_stims, savefile);
toc
%% 5) run create table tracked spines
disp('Running create table tracked spines')
savefile = 'tracked_spines.mat';
spinefile = 'spine_data_table.mat';
FOV_alignment_path = fullfile(homepath, path, 'demo_data', FOV_alignment_data_folder);
savepath = fullfile(homepath, path, 'demo_data', output_folder);
tic
create_table_tracked_spines_function(FOV_alignment_path, savepath, spinefile, savefile)
toc
