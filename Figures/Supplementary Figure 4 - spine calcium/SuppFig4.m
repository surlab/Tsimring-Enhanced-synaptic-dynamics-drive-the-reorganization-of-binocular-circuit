%% Figure S4
% Panel A: plot the calcium activity for ipsi, contra, binoc traces
% Panel C: 
% Panel B: Swarm chart comparing spine calcium activity on D1, D5, and D10
% Panel C: bar graph comparing spine calcium activity for unresponsive and
% responsive spines
%% Parameters for smoothing
clear all
alpha = 0.4; % for smoothing
std_thresh = 3; %find peaks threshold above 3 std
conseq_frames = 3; % number of frames 
fs = 7.68; %I should change this to actual frame rate


%% Panel A
path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
%path = 'G:/My Drive/';

load(fullfile(path,'spine_properties_table.mat'));
all_stim_table.roi_fovs_mouse_cell_days = strcat(num2str(all_stim_table.all_roi_inds), '_', ...
                                                all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);
temp = all_stim_table;
features = {'spine_area','all_roi_inds','roi_fovs_mouse_cell_days','resp', 'soma_resp', ...
    'spine_mean_resp','days','soma_mean_resp', 'all_active_spine_trials',...
    'all_active_spine_trials_smooth', 'all_z_scored_trace','mouse_cell', 'session'};
contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, "roi_fovs_mouse_cell_days");

%Calculate overall frequency of spine events 
feature = 'all_z_scored_trace';
concatenate_z_scored_trace = cellfun(@(x,y,z) [x;y;z], contra_ipsi_binoc.([feature,'_contra']),...
contra_ipsi_binoc.([feature,'_ipsi']),contra_ipsi_binoc.(feature), 'UniformOutput',false);
   
%smooth trace using exponential weighted filter 
concatenated_smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

%count spine calcium event rate
[pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), concatenated_smoothed_z_scored_trace,'UniformOutput',false);
ca_rate = cellfun(@(x,y) length(x)*fs/length(y), pks, concatenate_z_scored_trace);
contra_ipsi_binoc.ca_rate = ca_rate;
%%
% plot ca rate by day
days = {'D1', 'D5', 'D10'};
ca_rate_by_day = [];
groups = [];
figure("Position",[440,377,242,420])

for i = 1:length(days)
    contra_ipsi_binoc_days = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.days, days{i}),:);
    %ca_rate = log10(contra_ipsi_binoc_days.ca_rate);
    ca_rate = contra_ipsi_binoc_days.ca_rate;
    ca_rate(isinf(ca_rate)) = NaN;
    subplot(2,1,1)
    swarmchart(i*ones(size(ca_rate)), ca_rate, '.'); hold on;
    %histogram(ca_rate, 'Normalization','probability'); hold on;
    ca_rate_by_day = [ca_rate_by_day; ca_rate];
    groups = [groups; i*ones(size(ca_rate))];
    if i == length(days)
      xlim([0,4]);
      xticks([1:length(days)]);
      xticklabels(days)
    end
    ylim([0,0.2]);
end
set(gca, 'YScale', 'log')
ylim([0,0.2])

%%
[p, tbl,stats] = kruskalwallis(ca_rate_by_day,groups)
multcompare(stats)
%% Panel B: load data get spine calcium of spines by turnover
chronic_path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
load(fullfile(chronic_path,"tracked_spines_properties_table.mat") );

vars = {'D1','D5','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_D1_D5 = D1_D5_table(contains(D1_D5_table.structure_type, 'retained'),vars);
retain_D1_D5.fovs_mouse_cell = strcat(retain_D1_D5.all_mice_cells, "_", retain_D1_D5.all_fovs);
retain_D1_D5.D5_roi = num2str(retain_D1_D5.D5.all_roi_inds);
retain_D1_D5.fovs_mouse_cell_roi_session = strcat(retain_D1_D5.all_mice_cells, "_", retain_D1_D5.all_fovs,"_", retain_D1_D5.D5_roi, '_', retain_D1_D5.session);
retain_D1_D5.fovs_mouse_cell_roi = strcat(retain_D1_D5.all_mice_cells, "_", retain_D1_D5.all_fovs,"_", retain_D1_D5.D5_roi);

mice_cells_D1_D5 = unique(retain_D1_D5.fovs_mouse_cell);

vars = {'D5','D10','session', 'all_mice_cells', 'all_fovs', 'structure_type'};
retain_D5_D10 = D5_D10_table(contains(D5_D10_table.structure_type, 'retained'),vars);
retain_D5_D10.fovs_mouse_cell = strcat(retain_D5_D10.all_mice_cells, "_", retain_D5_D10.all_fovs);
mice_cells_D5_D10 = unique(retain_D5_D10.fovs_mouse_cell);
retain_D5_D10.D5_roi = num2str(retain_D5_D10.D5.all_roi_inds);
retain_D5_D10.fovs_mouse_cell_roi_session = strcat(retain_D5_D10.all_mice_cells, "_", retain_D5_D10.all_fovs,"_", retain_D5_D10.D5_roi, '_', retain_D5_D10.session);
retain_D5_D10.fovs_mouse_cell_roi = strcat(retain_D5_D10.all_mice_cells, "_", retain_D5_D10.all_fovs,"_", retain_D5_D10.D5_roi);
mice_cells_D1_D10 = intersect(mice_cells_D5_D10,mice_cells_D1_D5);
all_temps = {retain_D1_D5, retain_D5_D10};
%% Panel B: retained spine delta ca rate
feature = 'all_z_scored_trace';
all_temps = {retain_D1_D5, retain_D5_D10};
ca_rate_by_day = [];
groups = [];
days = {'D1', 'D5', 'D10'};
count = 1;
close all
figure("Position",[440,646,146,151])
for i = 1:length(all_temps)
    temp = all_temps{i};
    d1 = days{i};
    d5 = days{i+1};
    ds = [{d1};{d5}];
    %extract each session and then combine
    features = {'fovs_mouse_cell_roi', 'structure_type', d1, d5};
    %concatenate trace across all 3 sessions
    contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, "fovs_mouse_cell_roi");
  
    for ii = 1:length(ds)
        contra_ipsi_binoc.(ds{ii}).cib = arrayfun(@(x,y,z) [num2str(x);num2str(y);num2str(z)], contra_ipsi_binoc.([ds{ii},'_contra']).('nonresp'),...
        contra_ipsi_binoc.([ds{ii},'_ipsi']).('nonresp'),contra_ipsi_binoc.(ds{ii}).nonresp, 'UniformOutput',false );
  
        concatenate_z_scored_trace = cellfun(@(x,y,z) [x;y;z], contra_ipsi_binoc.([ds{ii},'_contra']).(feature),...
        contra_ipsi_binoc.([ds{ii},'_ipsi']).(feature),contra_ipsi_binoc.(ds{ii}).(feature), 'UniformOutput',false);
        
        %smooth trace using exponential weighted filter 
        smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

        %count spine calcium event rate
        [pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), smoothed_z_scored_trace,'UniformOutput',false);
        ca_rate = cellfun(@(x,y) length(x)*fs/length(y), pks, concatenate_z_scored_trace);
        
        ca_rate(isinf(ca_rate))=NaN;
        contra_ipsi_binoc.(ds{ii}).ca_rate = ca_rate;
    end
    contra_ipsi_binoc = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.(d1).cib, '111')&strcmp(contra_ipsi_binoc.(d5).cib, '111'),:);
    ca_rate_d1 = contra_ipsi_binoc.(d1).ca_rate;
    ca_rate_d5 = contra_ipsi_binoc.(d5).ca_rate;
    ca_rate_by_day = [ca_rate_by_day; [ca_rate_d5-ca_rate_d1]];
    [p,a]=signrank(ca_rate_d5-ca_rate_d1)
    groups = [groups; i*ones(size(ca_rate))];
    swarmchart(count*ones(size(ca_rate_d1)),[ca_rate_d5-ca_rate_d1], '.'); hold on;
    
    
    %ylim([-.0,0.02])
    count = count + 2;
end

xticks([1,3])
xlim([0,4])
xticklabels({[days{1},'-',days{2}],[days{2},'-',days{3}]});

%% Panel C: load data
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

data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};
structure = {"lost","formed", "lost","formed"};

%% Panel C: get spine calcium of spines by turnover all spines by session
feature = 'all_z_scored_trace';
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};
use_log = 0;
plot_bar = 1;
sessions = {'', '_contra','_ipsi'}
for s = 1:length(sessions)
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
   
        temp = data{i};
        temp.roi_fovs_mouse_cell = strcat(num2str(temp.all_roi_inds),'_', temp.all_fovs, '_', temp.all_mice_cells);
    
        %extract each session and then combine
        features = {'nonresp','resp',feature,'roi_fovs_mouse_cell', 'structure_type'};
        contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, "roi_fovs_mouse_cell");
        contra_ipsi_binoc.cib = arrayfun(@(x,y,z) [num2str(x);num2str(y);num2str(z)], contra_ipsi_binoc.('nonresp_contra'),...
        contra_ipsi_binoc.('nonresp_ipsi'),contra_ipsi_binoc.nonresp, 'UniformOutput',false );
    
        %concatenate trace across all 3 sessions
        concatenate_z_scored_trace = contra_ipsi_binoc.([feature,sessions{s}]);
        
        %smooth trace using exponential weighted filter 
        smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);
    
        %count spine calcium event rate
        [pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), smoothed_z_scored_trace,'UniformOutput',false);
         if use_log
         ca_rate = cellfun(@(x,y) log10(length(x)*fs/length(y)), pks, concatenate_z_scored_trace);
         ylims = [-3,-1];
        else
         ca_rate = cellfun(@(x,y) (length(x)*fs/length(y)), pks, concatenate_z_scored_trace);
         ylims = [0,0.015];
        end
        ca_rate(isinf(ca_rate))=NaN;
        %ca_rate(ca_rate>0.02)=0.02;
        contra_ipsi_binoc.ca_rate = ca_rate;

        data_feature_not_retained = [data_feature_not_retained, { contra_ipsi_binoc(ismember( contra_ipsi_binoc.structure_type, structure{i}),:).ca_rate}];
        data_feature_retained = [data_feature_retained, { contra_ipsi_binoc(~ismember( contra_ipsi_binoc.structure_type, structure{i}),:).ca_rate}];
    
    end
    figure("Position", [195,243,246,118])
    if plot_bar
        plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylims,[sessions{s}, 'Calcium Event Rate']);
    else
        plot_violinchart_lost_added_retained(data_feature_not_retained,data_feature_retained, [sessions{s}, 'Log_{10} Calcium Event Rate'], ylims);
    end
end
%% Panel F-G: get population level ca activity of spines by day
path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
load(fullfile(path,'spine_properties_table.mat'));
all_stim_table.roi_fovs_mouse_cell_days = strcat(num2str(all_stim_table.all_roi_inds), '_', ...
                                                all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);
temp = all_stim_table;
features = {'all_roi_inds','roi_fovs_mouse_cell_days','resp', 'soma_resp', ...
    'spine_mean_resp','days','soma_mean_resp', 'all_active_spine_trials',...
    'all_active_spine_trials_smooth', 'all_z_scored_trace','mouse_cell', 'session'};
contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, "roi_fovs_mouse_cell_days");

%Calculate overall frequency of spine events 
feature = 'all_z_scored_trace';
concatenate_z_scored_trace = cellfun(@(x,y,z) [x;y;z], contra_ipsi_binoc.([feature,'_contra']),...
contra_ipsi_binoc.([feature,'_ipsi']),contra_ipsi_binoc.(feature), 'UniformOutput',false);
   
%smooth trace using exponential weighted filter 
concatenated_smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

%count spine calcium event rate
[pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), concatenated_smoothed_z_scored_trace,'UniformOutput',false);
ca_rate = cellfun(@(x,y) length(x)*fs/length(y), pks, concatenate_z_scored_trace);
contra_ipsi_binoc.ca_rate_all = ca_rate;


%Calculate overall frequency of spine events 
feature = 'all_z_scored_trace';
concatenate_z_scored_trace = cellfun(@(x,y,z) [x], contra_ipsi_binoc.([feature,'_contra']),...
contra_ipsi_binoc.([feature,'_ipsi']),contra_ipsi_binoc.(feature), 'UniformOutput',false);
   
%smooth trace using exponential weighted filter 
concatenated_smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

%count spine calcium event rate
[pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), concatenated_smoothed_z_scored_trace,'UniformOutput',false);
ca_rate = cellfun(@(x,y) length(x)*fs/length(y), pks, concatenate_z_scored_trace);
contra_ipsi_binoc.ca_rate_contra = ca_rate;

%Calculate overall frequency of spine events 
feature = 'all_z_scored_trace';
concatenate_z_scored_trace = cellfun(@(x,y,z) [y], contra_ipsi_binoc.([feature,'_contra']),...
contra_ipsi_binoc.([feature,'_ipsi']),contra_ipsi_binoc.(feature), 'UniformOutput',false);
   
%smooth trace using exponential weighted filter 
concatenated_smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

%count spine calcium event rate
[pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), concatenated_smoothed_z_scored_trace,'UniformOutput',false);
ca_rate = cellfun(@(x,y) length(x)*fs/length(y), pks, concatenate_z_scored_trace);
contra_ipsi_binoc.ca_rate_ipsi = ca_rate;

%Calculate overall frequency of spine events 
feature = 'all_z_scored_trace';
concatenate_z_scored_trace = cellfun(@(x,y,z) [z], contra_ipsi_binoc.([feature,'_contra']),...
contra_ipsi_binoc.([feature,'_ipsi']),contra_ipsi_binoc.(feature), 'UniformOutput',false);
   
%smooth trace using exponential weighted filter 
concatenated_smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

%count spine calcium event rate
[pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), concatenated_smoothed_z_scored_trace,'UniformOutput',false);
ca_rate = cellfun(@(x,y) length(x)*fs/length(y), pks, concatenate_z_scored_trace);
contra_ipsi_binoc.ca_rate = ca_rate;

 contra_ipsi_binoc.ci = arrayfun(@(x,y) [num2str(x);num2str(y)], contra_ipsi_binoc.('resp_contra'),...
        contra_ipsi_binoc.('resp_ipsi'), 'UniformOutput',false );
   contra_ipsi_binoc.cib = arrayfun(@(x,y,z) [num2str(x);num2str(y);num2str(z)], contra_ipsi_binoc.('resp_contra'),...
        contra_ipsi_binoc.('resp_ipsi'),contra_ipsi_binoc.('resp'), 'UniformOutput',false );
%% Panel F: plot contra eyes ca rate 
% plot ca rate by day
days = {'D1', 'D5', 'D10'};
%close all

figure("Position",[440,522,344,146])
eye_type = {'10', '11'};
all_temp = [];
 all_eye_type = [];
all_day = [];
use_log = 0;
plot_bar = 1;
for i = 1:length(days)

    subplot(1, length(days),i)
    contra_ipsi_binoc_days = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.days, days{i}),:);
    
    
    for ii = 1:length(eye_type)
         temp = contra_ipsi_binoc_days(strcmp(contra_ipsi_binoc_days.ci, eye_type{ii}),:).ca_rate_contra;
        if use_log
            temp = log10(temp);
            temp(isinf(temp)) = NaN;
            ylims = [-3,-0.5];
        else
            ylims = [0,0.2];
        end

        if plot_bar
            bar(ii, nanmean(temp)); hold on;
            errorbar(ii, nanmean(temp), nanstd(temp./length(temp)), '.k'); hold on;
            scatter(ii*ones(size(temp)), temp, '.k', 'SizeData',1,'jitter', 'on', 'jitterAmount', 0.1);
            title(days{i})
            ylim(ylims);
            xlim([0,5])
            if ii == length(eye_type)
             xticks(1:length(eye_type))
             xticklabels(eye_type)
            end
        end
         %swarmchart(ii*ones(size(temp)), temp, '.k'); hold on;

        all_temp = [all_temp; temp];
        all_day = [all_day; repmat(days(i), size(temp))];
        all_eye_type = [all_eye_type; repmat(eye_type(ii), size(temp))];
   
        
    end

    if plot_bar == 0
        inds = strcmp(all_day, days(i));
        v = violinplot(all_temp(inds), all_eye_type(inds));
        colors = [0.38,1.00,0.38; 0.30,0.75,0.93; 0.98,0.41,0.41; 0.99,0.99,0.38];
         for i = 1:length(v)
            v(i).ViolinColor = colors(i,:);
            v(i).ViolinAlpha = 0.8;
            v(i).MedianPlot.SizeData = 20;
            v(i).BoxWidth = 0.08;
            v(i).ScatterPlot.Marker = 'o';
            v(i).ScatterPlot.SizeData = 0.2;
            v(i).ScatterPlot.MarkerFaceColor = 'k';
            v(i).ScatterPlot.MarkerFaceAlpha = 0.2;
         end
         ylim([-3,-0.5])
         xlim([0,3])
         xtickangle(45)      
    end


end
%%
[p,tbl_p, stats] = anovan(all_temp, {all_day, all_eye_type}, 'model', 'interaction')
[results,~,~,gnames]=multcompare(stats, 'Dimension', [1,2])
tbl2 = array2table(results,"VariableNames", ...
["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A")=gnames(tbl2.("Group A"));
tbl2.("Group B")=gnames(tbl2.("Group B"))

groupsummary(table(all_temp, all_day, all_eye_type),["all_day", 'all_eye_type'])
%% Panel G: Plot ipsi eyes ca rate 
% plot ca rate by day
days = {'D1', 'D5', 'D10'};
%close all

figure("Position",[440,522,344,146])
eye_type = {'01', '11'};
all_temp = [];
 all_eye_type = [];
all_day = [];
use_log = 0;
plot_bar = 1;
for i = 1:length(days)

    subplot(1, length(days),i)
    contra_ipsi_binoc_days = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.days, days{i}),:);
    
    
    for ii = 1:length(eye_type)
         temp = contra_ipsi_binoc_days(strcmp(contra_ipsi_binoc_days.ci, eye_type{ii}),:).ca_rate_ipsi;
        if use_log
            temp = log10(temp);
            temp(isinf(temp)) = NaN;
            ylims = [-3,-0.5];
        else
            ylims = [0,0.2];
        end

        if plot_bar
            bar(ii, nanmean(temp)); hold on;
            errorbar(ii, nanmean(temp), nanstd(temp./length(temp)), '.k'); hold on;
            scatter(ii*ones(size(temp)), temp, '.k', 'SizeData',1,'jitter', 'on', 'jitterAmount', 0.1);
            title(days{i})
            ylim(ylims);
            xlim([0,5])
            if ii == length(eye_type)
             xticks(1:length(eye_type))
             xticklabels(eye_type)
            end
        end
         %swarmchart(ii*ones(size(temp)), temp, '.k'); hold on;

        all_temp = [all_temp; temp];
        all_day = [all_day; repmat(days(i), size(temp))];
        all_eye_type = [all_eye_type; repmat(eye_type(ii), size(temp))];
   
        
    end

    if plot_bar == 0
        inds = strcmp(all_day, days(i));
        v = violinplot(all_temp(inds), all_eye_type(inds));
        colors = [0.38,1.00,0.38; 0.30,0.75,0.93; 0.98,0.41,0.41; 0.99,0.99,0.38];
         for i = 1:length(v)
            v(i).ViolinColor = colors(i,:);
            v(i).ViolinAlpha = 0.8;
            v(i).MedianPlot.SizeData = 20;
            v(i).BoxWidth = 0.08;
            v(i).ScatterPlot.Marker = 'o';
            v(i).ScatterPlot.SizeData = 0.2;
            v(i).ScatterPlot.MarkerFaceColor = 'k';
            v(i).ScatterPlot.MarkerFaceAlpha = 0.2;
         end
         ylim([-3,-0.5])
         xlim([0,3])
         xtickangle(45)      
    end


end
%%
[p,tbl_p, stats] = anovan(all_temp, {all_day, all_eye_type}, 'model', 'interaction')
[results,~,~,gnames]=multcompare(stats, 'Dimension', [1,2])
tbl2 = array2table(results,"VariableNames", ...
["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A")=gnames(tbl2.("Group A"));
tbl2.("Group B")=gnames(tbl2.("Group B"))

groupsummary(table(all_temp, all_day, all_eye_type),["all_day", 'all_eye_type'])
