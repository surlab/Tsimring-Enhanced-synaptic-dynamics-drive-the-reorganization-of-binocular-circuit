%% Load data (update paths)
clear all
chronic_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';
%chronic_path ='G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';
load(fullfile(chronic_path,"tracked_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_D1_D5_D10_lax_criteria_spine_area_dend_type_soma_props_trial_data.mat") );

vars = {'D1','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_vs_lost_D1_D5 = D1_D5_table(contains(D1_D5_table.structure_type, 'lost')|contains(D1_D5_table.structure_type, 'retained'),vars);
vars = {'D5','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_vs_formed_D1_D5 = D1_D5_table(contains(D1_D5_table.structure_type, 'retained')|contains(D1_D5_table.structure_type, 'formed'),vars);
vars = {'D5','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_vs_lost_D5_D10 = D5_D10_table(contains(D5_D10_table.structure_type, 'lost')|contains(D5_D10_table.structure_type, 'retained'),vars);
vars = {'D10','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_vs_formed_D5_D10 = D5_D10_table(contains(D5_D10_table.structure_type, 'retained')|contains(D5_D10_table.structure_type, 'formed'),vars);
%% get ori information
unique_orientations = [0:45:315];
oris_deg_psycopy = (180-unique_orientations); %% based on psycho conversion
oris_deg_psycopy(oris_deg_psycopy<0) = oris_deg_psycopy(oris_deg_psycopy<0)+360;
[inds1, inds] = sort(oris_deg_psycopy);
%% Plot fraction of trials that spine was active to each direction and soma's response to each direction during binoc viewing
% for lost spine on first time point
close all
mouse_cell = 'BM016-Cell4';
Dend = 'Dend3_0';
tbl = D5_D10_table;
d1 = 'D5';
d5 = 'D10';
close all
temp = tbl(strcmp(tbl.all_fovs,Dend)&strcmp(tbl.all_mice_cells, mouse_cell)&strcmp(tbl.session, "binoc"),:);
feature = 'all_active_spine_trials_smooth';
lost_spine = 16;
temp = temp(strcmp(temp.structure_type,'lost')&(temp.(d1).all_roi_inds)==lost_spine,:);
for i = 1:height(temp)
    figure('Position',[440,649,197,148])
    spine_trials = temp(i,:).(d1).(feature);
    inc_trials = temp(i,:).(d1).all_included_trial;
    soma_mean_resp = temp(i,:).(d1).soma_mean_resp;

    peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);
    sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
    spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
    fract_active_trials = cellfun(@(x,y,z) sum(squeeze(x),2)./z', spine_trials,peak_dir, sum_inc_trials,'UniformOutput',false);
    plot(soma_mean_resp{:}(inds)), hold on;
    ylabel('Soma mean amplitude')
    ylim([-1,6])
    yyaxis right
    plot(fract_active_trials{:}(inds))
    title(['lost:', num2str(temp(i,:).(d1).all_roi_inds)])
    ylim([-0.1,0.6])
    xlim([0,9])
    xticks([1:2:8])
    xticklabels([0:90:270])
end
%%
% for retained spine on first time point
temp = tbl(strcmp(tbl.all_fovs,Dend)&strcmp(tbl.all_mice_cells, mouse_cell)&strcmp(tbl.session, "binoc"),:);
feature = 'all_active_spine_trials_smooth';
retained_spine = 23;
temp = temp(strcmp(temp.structure_type,'retained')&(temp.(d1).all_roi_inds)==retained_spine,:);

for i = 1:(height(temp))
     figure('Position',[440,649,197,148])
    spine_trials = temp(i,:).(d1).(feature);
    inc_trials = temp(i,:).(d1).all_included_trial;
    soma_mean_resp = temp(i,:).(d1).soma_mean_resp;

    peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);
    sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
    spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
    fract_active_trials = cellfun(@(x,y,z) sum(squeeze(x),2)./z', spine_trials,peak_dir, sum_inc_trials,'UniformOutput',false);
    plot(soma_mean_resp{:}(inds)), hold on;
    ylabel('Soma mean amplitude')
    ylim([-1,6])
    yyaxis right
    plot(fract_active_trials{:}(inds))
    title(['retained:',num2str(temp(i,:).(d1).all_roi_inds)])
    ylim([-0.1,0.6])
    xlim([0,9])
    xticks([1:2:8])
    xticklabels([0:90:270])
end
%%
% for added spine on second time point
temp = tbl(strcmp(tbl.all_fovs,Dend)&strcmp(tbl.all_mice_cells, mouse_cell)&strcmp(tbl.session, "binoc"),:);
feature = 'all_active_spine_trials_smooth';
added_spine = 23;
temp = temp(strcmp(temp.structure_type,'formed')&(temp.(d5).all_roi_inds)==added_spine,:);

for i = 1:height(temp)
    figure('Position',[440,649,197,148])
    spine_trials = temp(i,:).(d5).(feature);
    inc_trials = temp(i,:).(d5).all_included_trial;
    soma_mean_resp = temp(i,:).(d5).soma_mean_resp;

    peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);
    sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
    spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
    fract_active_trials = cellfun(@(x,y,z) sum(squeeze(x),2)./z', spine_trials,peak_dir, sum_inc_trials,'UniformOutput',false);
    plot(soma_mean_resp{:}(inds)), hold on;
    ylabel('Soma mean amplitude')
    ylim([-1,6])
    yyaxis right
    plot(fract_active_trials{:}(inds))
    title(['added:',num2str(temp(i,:).(d5).all_roi_inds)])
    ylim([-0.1,0.6])
    xlim([0,9])
    xticks([1:2:8])
    xticklabels([0:90:270])
end

%%
% for retained spine on second timepoint (same spine as plotted on first
% timepoint)

temp = tbl(strcmp(tbl.all_fovs,Dend)&strcmp(tbl.all_mice_cells, mouse_cell)&strcmp(tbl.session, "binoc"),:);
feature = 'all_active_spine_trials_smooth';

retained_spine = 23;
temp = temp(strcmp(temp.structure_type,'retained')&(temp.(d1).all_roi_inds)==retained_spine,:);


for i = 1:height(temp)
     figure('Position',[440,649,197,148])
    spine_trials = temp(i,:).(d5).(feature);
    inc_trials = temp(i,:).(d5).all_included_trial;
    soma_mean_resp = temp(i,:).(d5).soma_mean_resp;

    peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);
    sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
    spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
    fract_active_trials = cellfun(@(x,y,z) sum(squeeze(x),2)./z', spine_trials,peak_dir, sum_inc_trials,'UniformOutput',false);
    plot(soma_mean_resp{:}(inds)), hold on;
    ylabel('Soma mean amplitude')
     ylim([-1,6])
    yyaxis right
    
    plot(fract_active_trials{:}(inds))
    %title(num2str(temp.D1.all_roi_inds))
    title(['retained:', num2str(temp(i,:).(d5).all_roi_inds)])
    ylim([-0.1,0.6])
    xlim([0,9])
    xticks([1:2:8])
    xticklabels([0:90:270])
end