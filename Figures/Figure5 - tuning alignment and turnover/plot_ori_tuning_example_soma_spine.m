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
%%
unique_orientations = [0:45:315];
oris_deg_psycopy = (180-unique_orientations); %% based on psycho conversion
oris_deg_psycopy(oris_deg_psycopy<0) = oris_deg_psycopy(oris_deg_psycopy<0)+360;
[oris_deg_sort,inds] = sort(oris_deg_psycopy);
oris_rad = (oris_deg_sort.*pi)./180;

%%
close all
mouse_cell = 'BM026-Cell2';
Dend = 'Dend1_0';
tbl = D5_D10_table;
d1 = 'D5';
d5 = 'D10';
close all
temp = tbl(strcmp(tbl.all_fovs,Dend)&strcmp(tbl.all_mice_cells, mouse_cell)&strcmp(tbl.session, "binoc"),:);
feature = 'all_active_spine_trials_smooth';

figure("Position", [283,290,176,168])
lost_temp = temp(strcmp(temp.structure_type,'lost'),:);
lost_roi = 12;
lost_temp = lost_temp(lost_temp.(d1).all_roi_inds == lost_roi,:);
soma_mean_resp = lost_temp.(d1).soma_mean_resp;
lost_spine_resp = lost_temp.(d1).spine_mean_resp;

ori_temp = [soma_mean_resp{:}];
plot_polarplots(ori_temp, oris_rad,'k', inds)

ori_temp = [lost_spine_resp{:}];
plot_polarplots(ori_temp, oris_rad,'g', inds)

retained_temp = temp(strcmp(temp.structure_type,'retained'),:);
retained_roi = 17;
retained_temp = retained_temp(retained_temp.(d1).all_roi_inds == retained_roi,:);
retained_spine_resp = retained_temp.(d1).spine_mean_resp;

ori_temp = [retained_spine_resp{:}];
plot_polarplots(ori_temp, oris_rad,'r', inds)
thetaticks([0:90:360])
thetaticklabels([90,135,0,45,180])

%%
figure("Position",[283,290,176,168])
added_temp = temp(strcmp(temp.structure_type,'formed'),:);
added_roi = 7;
added_temp = added_temp(added_temp.(d5).all_roi_inds == added_roi,:);
soma_mean_resp = added_temp.(d5).soma_mean_resp;
added_spine_resp = added_temp.(d5).spine_mean_resp;

ori_temp = [soma_mean_resp{:}];
plot_polarplots(ori_temp, oris_rad,'k', inds)

ori_temp = [added_spine_resp{:}];
plot_polarplots(ori_temp, oris_rad,'g', inds)

retained_temp = temp(strcmp(temp.structure_type,'retained'),:);
retained_roi = 17;
retained_temp = retained_temp(retained_temp.(d1).all_roi_inds == retained_roi,:);
retained_spine_resp = retained_temp.(d5).spine_mean_resp;

ori_temp = [retained_spine_resp{:}];
plot_polarplots(ori_temp, oris_rad,'r', inds)
thetaticks([0:90:360])
thetaticklabels([90,135,0,45,180])

%%
function plot_polarplots(ori_temp, oris_rad,c, inds)
    ori_temp(ori_temp<0)=0;
    ori_temp = ori_temp(inds);
    [OSI, Ori_pref, R_ori, mean_data_half, oris_rad_half] = get_orientation_tuning_vector(oris_rad, ori_temp);
    %polarplot([oris_rad_half*2,oris_rad_half(1)*2],[mean_data_half, mean_data_half(1)]/max(mean_data_half),c,'LineWidth',2); hold on;
    polarplot([0,mod(atan2(imag(R_ori),real(R_ori)),2*pi)], [0,OSI], c,'LineWidth',1); hold on;

    rlim([0, 0.7])
    disp(['OSI: ', num2str(OSI)])
    Ori_pref = Ori_pref-90;
    if Ori_pref<0
        Ori_pref = Ori_pref+180;
    end
    disp(['Ori pref: ', num2str(Ori_pref)])
    
end
