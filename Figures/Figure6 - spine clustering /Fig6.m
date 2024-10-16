%% Figure 6
% Panel B: Plot the distances between responsive spines (based on binoc
% response) for D1 and D10
% Panel C: plot the trial amp correlation of responsive spine pairs vs spine pair
% distance
% Panel D: plot the trial-trial amp correlation for all retained spines to
% nearst neighbor
% Panel E: plot the number of clusters and size of clusters by cell over
% development 


%% Panel B&C: Load data
clear all
path = 'G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';%update this path for where data will be stored
path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';

savepath = fullfile(path, 'Analyzed Data');
distance_path = fullfile(path, 'Chronic Imaging/FOV_alignment/');
load(fullfile(savepath, "spine_mean_tuning_table_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_lax_criteria_zscored_trace_active_trials.mat"), "all_stim_table");

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
feature = 'all_trial_amp';
soma_get_resp = 0;
soma_feature = 'all_dir_corr';
unique_sessions = {'binoc'};
unique_day = {'D1'; 'D5';'D10'};
% get spine corrs by spine distance
[all_mean_soma_resp_corrs_by_session, all_resp_corrs_by_session,all_resp_max_corrs_by_session,all_resps_dists_by_session,...
    all_dists_by_session, ~,all_unresp_dists_by_session, all_resp_unresp_dists_by_session]...
    = get_spine_pair_correlations_by_dist(feature,soma_feature, soma_get_resp, all_stim_table, unique_sessions, unique_day, distance_path);
%% Panel B&C: shuffle the spine corrs
bins = [0:5:20];
shuffle_times = 10000;
[CI_by_session, shuffled_means_by_sessions, shuffled_dists] = get_bootstrapped_CI_for_shuffled_corrs(bins, shuffle_times, all_resp_corrs_by_session, all_resps_dists_by_session, all_dists_by_session);
%% Panel B: plot the shuffled and true distances of responsive spine pairs 
close all
f = figure("Name", "Spine Pair Distance by Day")
bins = [5:0.1:10];
unique_day = {'D1'; 'D5';'D10'};
histogram_shuffled_and_resp_pair_distances(bins,all_resps_dists_by_session,shuffled_dists, unique_day);

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
%%
bins = [5:20];
unique_day = {'D1', 'D5','D10'};
compare_corrs_between_days(all_resp_corrs_by_session,all_resps_dists_by_session, bins,unique_day)

%% Panel E: determine the number and size of clusters per day
feature = 'all_trial_amp';
soma_get_resp = 0;
soma_feature = 'all_dir_corr';
unique_sessions = {'binoc'};
unique_day = {'D1'; 'D10'};
all_mouse_cell_by_session = get_spine_pair_correlations_by_dist_matrix_version(feature,soma_feature, soma_get_resp, all_stim_table, unique_sessions, unique_day, distance_path);

%% Panel E: generate graph based on spine pair within dist_tresh and above spine's
% mean corr to other spine pairs
dist_thresh = 5;
close all;
all_clust_by_dend_by_day = [];
all_size_clust_by_day = [];
count = 1;
for i = 1:length(all_mouse_cell_by_session)
    f = figure( "Name",['Number of Clusters per Cell per Day, Session: ', unique_sessions{i}], "Position",[440,639,128,158])
    g= figure( "Name",['Mean Cluster Size per Cell per Day, Session: ', unique_sessions{i}], "Position",[440,639,128,158])


    mouse_cell_by_day = all_mouse_cell_by_session.(unique_sessions{i});
    for d = 1:length(mouse_cell_by_day)
        mouse_cell = mouse_cell_by_day{d};
        all_num_clust = [];
        all_size_clust = [];
        all_num_dends = [];
        all_mouse_cell = [];
        all_soma_resp = [];
 
        for m = 1:length(mouse_cell)
            
            sum_clust = 0;
            size_clust = 0;
            sum_dist = 0;
            if ~isempty(find(ismember(fieldnames(mouse_cell{m}),{'all_resp_dists'})))
                 resp_corrs_by_branch = mouse_cell{m}.all_corrs;
                 resp_dists_by_branch= mouse_cell{m}.all_resp_dists;
                 all_dists_by_branch= mouse_cell{m}.all_dists;
                 sum_dists = [];
                 for r = 1:length(resp_corrs_by_branch)
                    resp_corrs =  resp_corrs_by_branch{r};
                    resp_dists = resp_dists_by_branch{r};
                    dist_thresh_pairs = resp_dists<dist_thresh;
                   
                    [mean_resp_corrs] = get_mean_corrs_matrix_exclude_diagonal(resp_corrs);
      
                    
                    check_corrs = resp_corrs;
                    check_corrs(~(dist_thresh_pairs)) = 0;
                    corr_thresh_pairs = get_thresh_corrs(check_corrs, mean_resp_corrs);

                         
                    resp_dists(~(dist_thresh_pairs&corr_thresh_pairs)) = 0;
                    G = digraph(resp_dists);
                    [~, binsizes] = conncomp(G);
                    sum_clust = sum_clust+sum(binsizes>1); %number of clusters with 2 or more spines
                    sum_dist = sum_dist+max(max(all_dists_by_branch{r}));
                    size_clust = [size_clust,binsizes(binsizes>1)]; %number of spines within cluster (2 or more)

                 end
                all_num_clust = [all_num_clust, sum_clust*30/sum_dist];
                all_soma_resp = [all_soma_resp, mouse_cell{m}.soma_resp];
                all_size_clust = [all_size_clust, mean(size_clust(size_clust>0))];
                all_num_dends = [all_num_dends, mouse_cell{m}.num_dends];
                all_mouse_cell = [all_mouse_cell, mouse_cell{m}.mouse_cell];

            end
            
        end
        
        figure(f)

        clust_by_dend = all_num_clust;
        bar(count,mean(clust_by_dend)); hold on;
        errorbar(count,mean(clust_by_dend),std(clust_by_dend)./sqrt(length(clust_by_dend)), 'k');
        swarmchart(count*ones(size(clust_by_dend)), clust_by_dend,'.');
        all_clust_by_dend_by_day = [all_clust_by_dend_by_day;{clust_by_dend}];
        
        figure(g)
        bar(count,nanmean(all_size_clust)); hold on;
        errorbar(count,nanmean(all_size_clust),nanstd(all_size_clust)./sqrt(length(all_size_clust)), 'k');
        swarmchart(count*ones(size(all_size_clust)), all_size_clust, '.');
        all_size_clust_by_day = [all_size_clust_by_day;{all_size_clust}];
        count = count + 1.5;
    end
    figure(f)
    xticks([1,2.5])
    xticklabels(unique_day)
    ylabel('# of Clusters (per 30um)')
    [p,a] = ranksum(all_clust_by_dend_by_day{1},all_clust_by_dend_by_day{2})
    plot_sig(1.5,10, p)
    figure(g)
    xticks([1,2.5])
    xticklabels(unique_day)
    ylabel('Avg Size of Clusters')
    [p,a] = ranksum(all_size_clust_by_day{1},all_size_clust_by_day {2})
    plot_sig(1.5,10, p)
end
%% Panel D: load data
clear all
chronic_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';
%chronic_path ='G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/';
load(fullfile(chronic_path,"tracked_spines_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_D1_D5_D10_lax_criteria_spine_area_dend_type_soma_props_trial_data.mat") );

%% Panel D: Plot corr to nearest spine neighbor for retained spines
% Obtain retained spines
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

all_temps = {retain_D1_D5, retain_D5_D10};
unique_sessions = {'binoc'};
days = {'D1', 'D5', 'D10'};
g = figure;
h = figure;
count = 1;

for i = 1:length(all_temps)
    count2 = 1;
    for s = 1:length(unique_sessions)
        temp = all_temps{i};
        temp = temp(strcmp(temp.session,unique_sessions{s}), :);
        temp = splitvars(temp);

        d1 = days{i};
        d5 = days{i+1};
        
        temp = temp(temp.([d1, '_resp'])>=0 &temp.([d5, '_resp'])>=0,:);
        %temp = temp(temp.([d1, '_neighbor_spine_resp'])>0 &temp.([d5, '_neighbor_spine_resp'])>0,:);
        
       temp = temp(temp.([d1, '_neighbor_spine_distance'])<5 &temp.([d5, '_neighbor_spine_distance'])<5,:);

        figure(g)
        subplot(2,length(unique_sessions), count2)
        dir_corr_d1 = temp.([d1, '_co_activity_trial_amp']);
        dir_corr_d5 = temp.([d5, '_co_activity_trial_amp']);
        for ii = 1:length(dir_corr_d1)
            plot([count, count+1],[dir_corr_d1(ii), dir_corr_d5(ii)], 'k'); hold on;
            scatter([count, count+1],[dir_corr_d1(ii), dir_corr_d5(ii)], '.k'); hold on;
        end
        [p,a]= ranksum(dir_corr_d1,dir_corr_d5)
        plot_sig(count + 0.5, 1, p);
        plot([count, count+1],[mean(dir_corr_d1), mean(dir_corr_d5)],'r', 'LineWidth',2);
        ylabel(['trial-trial signal corr'])
        title(unique_sessions{s})
        xticks([1,2,3,4])
        xlim([0,5])
        xticklabels({days{1},days{2},days{2},days{3} })
        ylim([-0.5,1])
        
        

        count2 = count2 + 1;
    end
    count = count + 2;
    
end

