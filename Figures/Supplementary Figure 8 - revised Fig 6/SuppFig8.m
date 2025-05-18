%% reanalyze Figure 6B, 6C, and 6E with neurons that imaged D1 and D10


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
%% get indices of neurons that were imaged on D1 and D10 
all_stim_table.mouse_cell_days = strcat(all_stim_table.mouse_cell,'_', all_stim_table.days);
repeated_image = zeros(height(all_stim_table),1);
unique_mouse_cell = unique(all_stim_table.mouse_cell);
for i = 1:length(unique_mouse_cell)
    inds_D1 = find(ismember(all_stim_table.mouse_cell_days, [unique_mouse_cell{i},'_D1']));
    inds_D5 = find(ismember(all_stim_table.mouse_cell_days, [unique_mouse_cell{i},'_D5']));
    inds_D10 = find(ismember(all_stim_table.mouse_cell_days, [unique_mouse_cell{i},'_D10']));
    if ~isempty(inds_D1)&~isempty(inds_D10)
       repeated_image(inds_D1) = 1;
        repeated_image(inds_D5) = 1;
       repeated_image(inds_D10) = 1;
    end
end

all_stim_table_filter = all_stim_table(repeated_image==1,:);

%%
feature = 'all_trial_amp';
soma_get_resp = 0;
soma_feature = 'all_dir_corr';
unique_sessions = {'binoc'};
unique_day = {'D1'; 'D5';'D10'};
% get spine corrs by spine distance
[all_mean_soma_resp_corrs_by_session, all_resp_corrs_by_session,all_resp_max_corrs_by_session,all_resps_dists_by_session,...
    all_dists_by_session,all_mouse_cell_by_session,all_unresp_dists_by_session, all_resp_unresp_dists_by_session]...
    = get_spine_pair_correlations_by_dist(feature,soma_feature, soma_get_resp, all_stim_table_filter, unique_sessions, unique_day, distance_path);
%% shuffle the spine corrs
bins = [0:5:20];
shuffle_times = 10000;
[CI_by_session, shuffled_means_by_sessions, shuffled_dists] = get_bootstrapped_CI_for_shuffled_corrs(bins, shuffle_times, all_resp_corrs_by_session, all_resps_dists_by_session, all_dists_by_session);
%% Panel A: get the shuffled and true distances of responsive spine pairs 
close all
f = figure("Name", "Spine Pair Distance by Day")
bins = [5:0.1:10];
unique_day = {'D1'; 'D5';'D10'};
histogram_shuffled_and_resp_pair_distances(bins,all_resps_dists_by_session,shuffled_dists, unique_day);

%% Panel A: compare median distances of D1 to D10 neurons
mouse_cell_session = all_mouse_cell_by_session.binoc;
all_resps_dists = all_resps_dists_by_session.binoc;
all_resps_dists_mouse_day = [];
unique_mouse_cell = intersect(mouse_cell_session{1},mouse_cell_session{3});
for i = [1,3]
    resp_dists = all_resps_dists{i};
    mouse_cell = mouse_cell_session{i};

    resp_dists_all_mouse =  [];
    for ii = 1:length(unique_mouse_cell)
        inds = find(strcmp(mouse_cell, unique_mouse_cell(ii)));
        resp_dists_mouse = resp_dists(inds);
        median_resp_dists = median(cell2mat(resp_dists_mouse));
        resp_dists_all_mouse = [resp_dists_all_mouse; median_resp_dists];
    end
    all_resps_dists_mouse_day = [all_resps_dists_mouse_day,resp_dists_all_mouse];
end

bar([mean(all_resps_dists_mouse_day(:,1)), mean(all_resps_dists_mouse_day(:,2))]); hold on;
for i = 1:size(all_resps_dists_mouse_day,1)
    plot([1,2], [all_resps_dists_mouse_day(i,1),   all_resps_dists_mouse_day(i,2)], 'k');
    scatter([1,2], [all_resps_dists_mouse_day(i,1),   all_resps_dists_mouse_day(i,2)], 'black', 'filled');
end
ylabel('spine pair distance')
xticks([1,2])
xticklabels({'D1', 'D10'})
%% Panel B: plot the trial-trial amplitude corr of responsive spines pairs
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
bins = [0:5:20];
unique_day = {'D1', 'D5','D10'};
compare_corrs_between_days(all_resp_corrs_by_session,all_resps_dists_by_session, bins,unique_day)
%% Panel D: determine the number and size of clusters per day
feature = 'all_trial_amp';
soma_get_resp = 0;
soma_feature = 'all_dir_corr';
unique_sessions = {'binoc'};
unique_day = {'D1'; 'D10'};
all_mouse_cell_by_session = get_spine_pair_correlations_by_dist_matrix_version(feature,soma_feature, soma_get_resp, all_stim_table_filter, unique_sessions, unique_day, distance_path);
%% Panel C
% generate graph based on spine pair within dist_tresh and above spine's
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
%         scatter(d*ones(size(clust_by_dend(all_soma_resp ==0))),clust_by_dend(all_soma_resp ==0),'k','jitter', 'on','jitterAmount', 0.2, 'SizeData',4);
%         scatter(d*ones(size(clust_by_dend(all_soma_resp ==1))),clust_by_dend(all_soma_resp ==1),'.','k','jitter', 'on','jitterAmount', 0.2);
        swarmchart(count*ones(size(clust_by_dend)), clust_by_dend,'.');
        all_clust_by_dend_by_day = [all_clust_by_dend_by_day;{clust_by_dend}];
        
        figure(g)
        bar(count,nanmean(all_size_clust)); hold on;
        errorbar(count,nanmean(all_size_clust),nanstd(all_size_clust)./sqrt(length(all_size_clust)), 'k');
%         scatter(d*ones(size(all_size_clust(all_soma_resp ==0))),all_size_clust(all_soma_resp ==0),'k','jitter', 'on','jitterAmount', 0.2,'SizeData',4);
%         scatter(d*ones(size(all_size_clust(all_soma_resp ==1))),all_size_clust(all_soma_resp ==1),'.','k','jitter', 'on','jitterAmount', 0.2);
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