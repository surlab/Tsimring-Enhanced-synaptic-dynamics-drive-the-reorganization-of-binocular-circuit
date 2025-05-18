%% Supplementary Figure 7
% Panel A: plot the traces for 4 example spines (need to find the example
% figure I used), and the correlations for the trial amplitudes across the
% 10 trials vs the spines distances (heatmap of matrices)
% Panel B: Plot the distances between responsive spines (based on binoc
% response) for D1 and D10?
% Panel C: plot the trial amp correlation of spine pairs vs spine pair
% distance
% Panel D: plot the number of clusters and size of clusters by cell over
% development 
% Panel E: plot the trial-trial amp correlation for lost, added, and
% retained spines


%% Load data
clear all
path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
distance_path = fullfile(path, 'Chronic Imaging/FOV_alignment/'); %update path
load(fullfile(path, "spine_properties_table.mat"), "all_stim_table");

% merge table into contra and ipsi session per spine
all_stim_table.roi_fovs_mouse_cell_days = strcat(num2str(all_stim_table.all_roi_inds), '_', ...
                                                all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);
temp = all_stim_table;
contra = temp(strcmp([temp.session{:}], "contra"),:);
contra = contra(:, {'resp','roi_fovs_mouse_cell_days'});
ipsi = temp(strcmp([temp.session{:}], "ipsi"),:);
ipsi = ipsi(:, {'resp','roi_fovs_mouse_cell_days','days', 'all_roi_inds', 'mouse_cell', 'all_fovs', 'all_day'});

contra_ipsi = outerjoin(contra, ipsi, "Keys","roi_fovs_mouse_cell_days", "MergeKeys",true);

%remove rows with NaN values 
contra_ipsi=contra_ipsi(~any(ismissing(contra_ipsi),2),:);

%create one hot encoding for spines and soma eye pref
contra_ipsi.ci_spine = strcat(num2str(contra_ipsi.resp_contra), num2str(contra_ipsi.resp_ipsi));
%% Panel B&C: get distances of responsive spine pairs, unresponsive spine pairs, and resp-unresp spine pairs
feature = 'spine_mean_resp';
soma_get_resp = 0;
soma_feature = 'all_dir_corr';
unique_sessions = {'binoc'};
unique_day = {'D1'; 'D5';'D10'};
% get spine corrs by spine distance
[all_mean_soma_resp_corrs_by_session, all_resp_corrs_by_session,all_resp_max_corrs_by_session,all_resps_dists_by_session,...
    all_dists_by_session, ~,all_unresp_dists_by_session, all_resp_unresp_dists_by_session]...
    = get_spine_pair_correlations_by_dist(feature,soma_feature, soma_get_resp, all_stim_table, unique_sessions, unique_day, distance_path);
%%
% shuffle the spine corrs
bins = [0:5:20];
shuffle_times = 10000;
[CI_by_session, shuffled_means_by_sessions, shuffled_dists] = get_bootstrapped_CI_for_shuffled_corrs(bins, shuffle_times, all_resp_corrs_by_session, all_resps_dists_by_session, all_dists_by_session);
%% Panel C: plot the trial-trial amplitude corr of responsive spines pairs
close all
% plot the spine corrs
bins = [0:5:20];
f = figure("Name", "mean clustering by day")
set(gcf, 'Position',[1141.66666666667,435,800.666666666667,554.666666666667]); 
g = figure("Name", "all clustering by day")
set(gcf, 'Position',[1141.66666666667,435,800.666666666667,554.666666666667]);
h = figure("Name", "delta clustering by day")
set(gcf, 'Position',[1141.66666666667,435,800.666666666667,554.666666666667]);

plot_spine_corr_by_dist(f,g,h, all_resp_corrs_by_session,all_resps_dists_by_session,CI_by_session,shuffled_means_by_sessions,bins, unique_day,feature);

%% Panel B
unique_sessions = {'binoc'};
unique_day = {'D1'; 'D5'; 'D10'};
close all
for s = 1:length(unique_sessions)
figure
all_resp_dists = all_resps_dists_by_session.(unique_sessions{s});
all_unresp_dists = all_unresp_dists_by_session.(unique_sessions{s});
all_resp_unresp_dists = all_resp_unresp_dists_by_session.(unique_sessions{s});
all_dists = [];
all_day = [];
all_type = [];
for d = [1,3]
    resp_dists = [cell2mat(all_resp_dists{d})];
    unresp_dists = [cell2mat(all_unresp_dists{d})];
    resp_unresp_dists = [cell2mat(all_resp_unresp_dists{d})];
    temp = [resp_dists;unresp_dists;resp_unresp_dists];
    type = [repmat(1,length(resp_dists),1);repmat(2,length(unresp_dists),1);repmat(3,length(resp_unresp_dists),1)];
      all_dists = [all_dists;resp_dists;unresp_dists;resp_unresp_dists];
    all_day = [all_day; repmat(unique_day(d), length([resp_dists;unresp_dists;resp_unresp_dists]),1)];
    all_type = [all_type;repmat({'r'},length(resp_dists),1);repmat({'u'},length(unresp_dists),1);repmat({'ru'},length(resp_unresp_dists),1)];
end
end
boxchart(categorical(all_type), all_dists,'GroupByColor',all_day, 'Notch', 'on');
ylim([0,40])
%%
tbl = table(all_dists, all_type, all_day);
unique_type = unique(all_type);
for i = 1:length(unique_type)
    temp = tbl(strcmp(tbl.all_type, unique_type(i)),:);
    temp(any(isnan(temp.all_dists)),:) = [];
    [p,a] = ranksum(temp(strcmp(temp.all_day, 'D1'),:).all_dists,temp(strcmp(temp.all_day, 'D10'),:).all_dists)
end
%%
[p,stats,tbl] = anovan(all_dists, {all_day, all_type}, 'model', 'interaction');
[results,~,~,gnames]=multcompare(tbl, 'Dimension', [1,2])
tbl2 = array2table(results,"VariableNames", ...
["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A")=gnames(tbl2.("Group A"));
tbl2.("Group B")=gnames(tbl2.("Group B"))

%% Panel B: load data
clear all
chronic_path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
load(fullfile(chronic_path,"tracked_spines_properties_table.mat") );

vars = {'D1','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_vs_lost_D1_D5 = D1_D5_table(contains(D1_D5_table.structure_type, 'lost')|contains(D1_D5_table.structure_type, 'retained'),vars);
vars = {'D5','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_vs_formed_D1_D5 = D1_D5_table(contains(D1_D5_table.structure_type, 'retained')|contains(D1_D5_table.structure_type, 'formed'),vars);
vars = {'D5','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_vs_lost_D5_D10 = D5_D10_table(contains(D5_D10_table.structure_type, 'lost')|contains(D5_D10_table.structure_type, 'retained'),vars);
vars = {'D10','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_vs_formed_D5_D10 = D5_D10_table(contains(D5_D10_table.structure_type, 'retained')|contains(D5_D10_table.structure_type, 'formed'),vars);
%% Panel C: get distance to nearest responsive neighbor
close all
session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};
days = {'D1', 'D5', 'D5', 'D10'};
structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,20];

for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    all_spines = all_stim_table(strcmp([all_stim_table.session{:}],session{s}),:);
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        %temp = temp(temp.neighbor_spine_distance<=5&temp.resp>0&temp.neighbor_spine_resp>0,:);        
        temp = temp(temp.resp>=0&temp.neighbor_spine_resp>=0,:); 
        d = days{i};
        distances = find_nearest_unresp_resp_neighbor(temp,all_spines,d, distance_path);
        data_feature_not_retained = [data_feature_not_retained, {distances(ismember(distances.structure_type, structure{i}),:).distances}];
        data_feature_retained = [data_feature_retained, {distances(~ismember(distances.structure_type, structure{i}),:).distances}];
    end
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,['Distance to nearest resp spine (um)']);

end


%% Panel D: Plot corr to nearest spine neighbor data
close all
feature = 'co_activity_trial_amp';

session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,0.4];

for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        temp = temp(temp.resp>=0&temp.neighbor_spine_resp>=0,:);        
        %temp = temp(temp.neighbor_spine_resp>0,:);    
 
        data_feature_not_retained = [data_feature_not_retained, {temp(ismember(temp.structure_type, structure{i}),:).(feature)}];
        data_feature_retained = [data_feature_retained, {temp(~ismember(temp.structure_type, structure{i}),:).(feature)}];
    end
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature]);

end

%% Panel E: get distances of contra, ipsi, contra+ipsi spine pairs
elements = {0:1, 0:1, 0:1, 0:1}; %cell array with N vectors to combine
temp = get_combin(elements);
combinations = [];
for i = 1:size(temp,1)
    combinations = [combinations; {[num2str(temp(i,1)), num2str(temp(i,2)), '_', num2str(temp(i,3)), num2str(temp(i,4))]}];
end
%%
repeat_elems = {'10_01', '10_11', '01_11'};
%remove any unresponsive pair combos and repeated pairs
combinations = combinations(~contains(combinations,'00'));
for i = 1:length(repeat_elems)
    combinations = combinations(~contains(combinations,repeat_elems(i)));
end
unique_day = {'D1';'D10'};
all_combo_dists = get_eye_specific_spine_pair_dists(contra_ipsi, unique_day, distance_path, combinations);
%% Panel E: plot distribution between contra, ipsi, contra+ipsi spine pairs
close all
same_type_inds = [1,4,6]; % I am too lazy to figure out how split the combo pairs
mixed_type_inds = [2,3,5];

f = figure("Name",['Same Eye-Specific Identity Spine Pair'], "Position",[440,612,254,185]);
plot_pair_dists(all_combo_dists,combinations, same_type_inds, unique_day);

g = figure("Name",['Mixed Eye-Specific Identity Spine Pair'], "Position",[440,612,254,185]);
plot_pair_dists(all_combo_dists,combinations, mixed_type_inds, unique_day);
%%
%run ANOVA to see if there is any differences between days and eye-pairs
figure
tbl_anova = create_table_for_anova(all_combo_dists, combinations,unique_day);
[~,~, stats]= anovan(tbl_anova.all_data, {tbl_anova.all_day, tbl_anova.all_eye_combo},'varnames',{'all_day','all_eye_combo'});
[results,~,~,gnames]=multcompare(stats,"Dimension", [1,2])
tbl2 = array2table(results,"VariableNames", ...
["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A")=gnames(tbl2.("Group A"));
tbl2.("Group B")=gnames(tbl2.("Group B"))