% Get bootrapped distribution for 1 trial per direction and
% save for chronically tracked spines

%% Load data
clear all
chronic_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';
%chronic_path ='G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';
load(fullfile(chronic_path,"tracked_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_D1_D5_D10_lax_criteria_spine_area_dend_type_soma_props_trial_data.mat") );

vars = {'D1','D5','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_D1_D5 = D1_D5_table(contains(D1_D5_table.structure_type, 'retained'),vars);
retain_D1_D5.fovs_mouse_cell = strcat(retain_D1_D5.all_mice_cells, "_", retain_D1_D5.all_fovs);
retain_D1_D5.D5_roi = num2str(retain_D1_D5.D5.all_roi_inds);
retain_D1_D5.fovs_mouse_cell_roi = strcat(retain_D1_D5.all_mice_cells, "_", retain_D1_D5.all_fovs,"_", retain_D1_D5.D5_roi, '_', retain_D1_D5.session);
mice_cells_D1_D5 = unique(retain_D1_D5.fovs_mouse_cell);

vars = {'D5','D10','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_D5_D10 = D5_D10_table(contains(D5_D10_table.structure_type, 'retained'),vars);
retain_D5_D10.fovs_mouse_cell = strcat(retain_D5_D10.all_mice_cells, "_", retain_D5_D10.all_fovs);
mice_cells_D5_D10 = unique(retain_D5_D10.fovs_mouse_cell);
retain_D5_D10.D5_roi = num2str(retain_D5_D10.D5.all_roi_inds);
retain_D5_D10.fovs_mouse_cell_roi = strcat(retain_D5_D10.all_mice_cells, "_", retain_D5_D10.all_fovs,"_", retain_D5_D10.D5_roi, '_', retain_D5_D10.session);
mice_cells_D1_D10 = intersect(mice_cells_D5_D10,mice_cells_D1_D5);

%% Get bootrapped distribution for 1 trial per direction and save for chronically tracked spines
fract_trials = 0.1;
all_temps = {retain_D1_D5, retain_D5_D10};
save_temp = [];
bootstrap = 10000;
day = {'D1', 'D5', 'D10'};
for i = 1:length(all_temps)
    temp = all_temps{i};
    d1 = day{i};
    d5 = day{i+1};
    ds = {d1,d5};
    all_mean_amp_boot_d1 = [];
    all_mean_amp_boot_d5 = [];
    
    for roi = 1:height(temp)
        roi_temp = temp(roi,:);
        if mod(roi,10)==0
            disp(["Running roi: ", num2str(roi)])
        end
        
        if roi_temp.(d1).resp>0 & roi_temp.(d5).resp>0
            % get D1
            trial_amp_roi = squeeze(roi_temp.(d1).all_trial_amp{:});
            [num_oris,~] = size(trial_amp_roi);
            trials = squeeze(roi_temp.(d1).all_included_trial{:});
            boot_mean_amp = do_boostrap(trial_amp_roi,trials,fract_trials,num_oris,bootstrap);
            all_mean_amp_boot_d1 = [all_mean_amp_boot_d1; {boot_mean_amp}];

            % get D5
             trial_amp_roi = squeeze(roi_temp.(d5).all_trial_amp{:});
            [num_oris,~] = size(trial_amp_roi);
            trials = squeeze(roi_temp.(d5).all_included_trial{:});
            boot_mean_amp = do_boostrap(trial_amp_roi,trials,fract_trials,num_oris,bootstrap);
            all_mean_amp_boot_d5 = [all_mean_amp_boot_d5; {boot_mean_amp}];


        else
            all_mean_amp_boot_d1 = [all_mean_amp_boot_d1; {NaN}];
            all_mean_amp_boot_d5 = [all_mean_amp_boot_d5; {NaN}];
        end
    end
    temp.(d1).all_mean_amp_boot = all_mean_amp_boot_d1; 
    temp.(d5).all_mean_amp_boot = all_mean_amp_boot_d5;

    save_temp = [save_temp, {temp}];

end

%%
chronic_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';
%path = 'G:/My Drive/';
save(fullfile(chronic_path, 'tracked_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_D1_D5_D10_lax_criteria_spine_area_dend_type_soma_props_trial_data_boot.mat'), 'save_temp', '-v7.3')

%%
function boot_mean_amp = do_boostrap(trial_amp_roi,trials, fract_trials, num_oris, bootstrap)
        boot_mean_amp=zeros(num_oris,bootstrap);
        for ori = 1:num_oris
            trials_ori = trials(ori,:);
            inc_trials = find(trials_ori);
            num_trials = length(inc_trials);
            inds1 = arrayfun(@(x) randsample(inc_trials,ceil(fract_trials*num_trials)),[1:bootstrap], 'UniformOutput', false);
            boot_mean_amp(ori,:) = cellfun(@(x) mean(trial_amp_roi(ori,x)),inds1);
        end
end