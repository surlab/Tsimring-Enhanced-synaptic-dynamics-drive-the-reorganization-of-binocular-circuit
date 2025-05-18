%% Figure 7
% Panel D (left): Plot the spine activity of pooled experimental data
% Panel E (left): plot the trial amp correlation of responsive spine pairs vs spine pair
% distance
% Panel F: Plot the delta mean corr between spine pairs by distance for
% experimental data and simulation data
%% Panel D (left): Plot the spine activity of pooled experimental data
load('ca_event_rate_retained_lost_added.mat')
data_retain = [data_feature_retained{1};data_feature_retained{3}];
data_lost = [data_feature_not_retained{1};data_feature_not_retained{3}];

figure
bar([mean(data_lost), mean(data_retain)]); hold on;
errorbar([mean(data_lost), mean(data_retain)],[std(data_lost)/sqrt(length(data_lost)), std(data_retain)/sqrt(length(data_retain))], '.k'); hold on;
ylim([0,0.02])
ranksum(data_lost,data_retain)
%% Panel E (left): Plot the fraction of spines active at soma's preferred direction for pooled experimental data
load('active_trials_retained_lost_added.mat')

data_retain = [data_feature_retained{1};data_feature_retained{3}];
data_lost = [data_feature_not_retained{1};data_feature_not_retained{3}];

figure

bar([mean(data_lost), mean(data_retain)]); hold on;
errorbar([mean(data_lost), mean(data_retain)],[std(data_lost)/sqrt(length(data_lost)), std(data_retain)/sqrt(length(data_retain))], '.k'); hold on;
set(gcf, 'Position',[191,408,98,141]);
ylim([0,40])
ylabel('% active trials (soma''s pref dir)')
ranksum(data_lost,data_retain)
%% Panel F: load expiremental data and simulation data 

% get shuffled correlations for corrs_aft in model
load('correlation_analysis.mat', 'all_resps_dists_by_session'); %experimental data 
load('pair_dists_corrs_correct.mat') %simulation data
load('pair_precisedists_corrs.mat', 'precise_dists') %simulation data

%% Panel F: Get delta mean corr per distance for experimental data (D10-D1) and simulation data (post-pre)
bins = [0:5:20];
dists_by_days = all_resps_dists_by_session.binoc;
corrs_by_days = all_resp_corrs_by_session.binoc;

[all_mean_grids] = compare_model_data_delta_corr(corrs_by_days,dists_by_days, bins);

data_delta_corrs_d1_d10 = (all_mean_grids(3,:)-all_mean_grids(1,:));

[all_mean_grids_model] = compare_model_data_delta_corr({corrs_bef,corrs_aft},{precise_dists,precise_dists}, bins);
model_delta_corrs_pre_post = all_mean_grids_model(2,:)-all_mean_grids_model(1,:);

%% Panel F: Plot delta mean corr per distance for experimental and simulation data
figure
yyaxis left
plot(bins(2:end), model_delta_corrs_pre_post), hold on;
ylabel('\Delta mean corr: post-pre model')
yyaxis right
scatter(bins(2:end),data_delta_corrs_d1_d10, 'k', 'filled')
ylim([0,0.12])
xlim([4,21])
xlabel('spine pair distance')

