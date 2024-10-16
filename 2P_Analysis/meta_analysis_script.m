% 1) Pre-process and analyze spine 2p data to obtain FOV's zscored data + OSI/Ori pref
% 2) collect all data into a structure array for each visual stim
% 3) identify visual responsive spines based on criteria
% 4) output all data into a table where each row is a spine, and the columns indicate different parameters measured for the spine for specific visual stim (binoc, contra, or ipsi)
% 5) output chronically tracked spines into two tables: spines tracked from D1 to D5 and spines tracked from D5 to D10

%% Set path
clear all
input_path = 'F:\Binocular_matching\Spine_imaging\2P_imaging\';
stim_path =  '/Users/ktsimring/Dropbox (MIT)/Katya Tsimring/stim_events/Binocular Matching Katya';
path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';
distance_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';%distance_path = 'G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';
%% 1) run analysis spine 2p data function
disp('Running analysis spine 2p data function')
mouse_files = {'BM014' ,'BM015', 'BM016', 'BM017', 'BM018', 'BM019', 'BM020', 'BM021', 'BM023','BM024', 'BM025', 'BM026', 'BM027','BM029', 'BM030'};
analysis_spine_2p_data_function(input_path, stim_path, path, mouse_files)

%% 2) run analyze processed data longstim function
disp('Running analysis processed data function')
mouse_files = {'BM014' ,'BM015', 'BM016', 'BM017', 'BM018', 'BM019', 'BM020', 'BM021', 'BM023','BM024', 'BM025', 'BM026', 'BM027','BM029', 'BM030'};
%mouse_files = { 'BM021'};

savefile = 'all_analyzed_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_zscored_trace_active_trials.mat';
%savefile = 'all_analyzed_spines_BM021_zscored_trace.mat';

analyzed_processed_data_function(path, savefile, mouse_files)

%% 3) run identify visual responses
disp('Running identify vis responses function')
analyzefile =  'all_analyzed_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_zscored_trace_active_trials.mat';
savefile = 'visually_responsive_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_lax_criteria.mat';
mean_amp_thresh = 0.5;
p_val_thresh = 0.05;
identify_vis_responses(path,analyzefile, savefile, mean_amp_thresh,p_val_thresh);

%% 4) run create table for all 3 sessions per spines
disp('Running create table for all 3 sessions')
analyzefile = 'mean_tuning_all_analyzed_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_zscored_trace_active_trials.mat';
visfile = 'visually_responsive_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_lax_criteria.mat';
savefile = 'mean_tuning_table_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_lax_criteria_zscored_trace_active_trials';
savefile_mat = [savefile, '.mat'];
savefile_csv = [savefile, '.csv'];
savepath = fullfile(path, 'Analyzed Data');
load(fullfile(savepath,analyzefile));
load(fullfile(savepath,visfile));
create_table_of_all3sessions_function(savepath, distance_path, all_stims, all_vis_stims, savefile_mat, savefile_csv);

%% 5) run create table tracked spines
disp('Running create table tracked spines')
savefile = 'tracked_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_D1_D5_D10_lax_criteria_spine_area_dend_type_soma_props_trial_data.mat';
spinefile = 'mean_tuning_table_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_lax_criteria_zscored_trace_active_trials.mat';
spinefile = ['spine_',spinefile];
create_table_tracked_spines_function(distance_path, path, spinefile, savefile)
