%% Figure 2 
% Panel A-B: plot the calcium activity for ipsi, contra, binoc traces
% Panel C: Bar chart of spine calcium activity for lost, retained, and
% added spines
% Panel D: CDF of spine calcium activity between responsive and
% unresponsive spines
% Panel E: Same as C but for responsive spines
% Panel F: Same as C but for unresponsive spines
% Panel G: Comparing binoc spine calcium activity by day for spines that are
% responsive to at least binoc 
%% Parameters for smoothing
clear all
alpha = 0.4; % for smoothing
std_thresh = 3; %find peaks threshold above 3 std
conseq_frames = 3; % number of frames 
fs = 7.68; 
%% Panel A and B
% get path details
close all
path = '/Volumes/GoogleDrive-108846495442099470486/';
fov_path = 'My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Processed Data/'; %update this path for where data will be stored

mouse = 'BM015';
day = 'p24_072922';
cell_name = 'Cell8';
dend_name = 'Dend1_0';
rois = [5,19];

dend_name_short = split(dend_name, '_');
dend_name_short = dend_name_short{1};
dend_binoc = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_binoc'],'normalized_data_by_stim.mat'), 'ROIdata_Z', 'timestamps');
dend_contra = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_contra'],'normalized_data_by_stim.mat'),'ROIdata_Z','timestamps');
dend_ipsi = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_ipsi'],'normalized_data_by_stim.mat'),'ROIdata_Z','timestamps');
sessions = {dend_binoc, dend_contra, dend_ipsi};
all_ca_rate = [];
for i = rois
    all_data = [];
    num_time = 0;
    for ii = 1:length(sessions)
        all_data = [all_data;sessions{ii}.ROIdata_Z(:,i)];
        num_time = num_time+sessions{ii}.timestamps(end);
    end
    [ca_rate, ~, ~, ~] = get_ca_rate(all_data, num_time,alpha, std_thresh, conseq_frames);
    all_ca_rate = [all_ca_rate, ca_rate]
end
h = heatmap(all_ca_rate); h.XDisplayLabels = rois;
colormap parula
h.ColorLimits = [0,0.02];

%% Panel B: 
close all
colors = {'g','r', 'b'};
rois = [5,19];
for i = rois
    figure("Position",[5,298,1162,420])
    for ii = 1:length(sessions)
        subplot(2,3,ii)
        data = sessions{ii};
        smooth_data = exp_weight_filter(data.ROIdata_Z(:,i), alpha);
        timestamps = data.timestamps;
        [pks,lks] = findpeaks(smooth_data, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames);  
        plot(timestamps,smooth_data, colors{ii}); hold on;
        scatter(timestamps(lks), pks, '*', colors{ii});
        plot(timestamps,std_thresh*ones(size(timestamps)), 'k');
        ylim([-2,14])

    end

end


%% Panel C: load data
chronic_path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files'; %update this path for where data will be stored
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

%% Panel C: get spine calcium of spines by turnover all spines
feature = 'all_z_scored_trace';
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};
structure = {"lost","formed", "lost","formed"};
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
    concatenate_z_scored_trace = cellfun(@(x,y,z) [x;y;z], contra_ipsi_binoc.([feature,'_contra']),...
        contra_ipsi_binoc.([feature,'_ipsi']),contra_ipsi_binoc.(feature), 'UniformOutput',false);
    
    %smooth trace using exponential weighted filter 
    smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

    %count spine calcium event rate

    [pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), smoothed_z_scored_trace,'UniformOutput',false);

    ca_rate = cellfun(@(x,y) (length(x)*fs/length(y)), pks, concatenate_z_scored_trace);
    ylims = [0,0.015];
    ca_rate(isinf(ca_rate))=NaN;
    %ca_rate(ca_rate>0.02)=0.02;
    contra_ipsi_binoc.ca_rate = ca_rate;

    data_feature_not_retained = [data_feature_not_retained, { contra_ipsi_binoc(ismember( contra_ipsi_binoc.structure_type, structure{i}),:).ca_rate}];
    data_feature_retained = [data_feature_retained, { contra_ipsi_binoc(~ismember( contra_ipsi_binoc.structure_type, structure{i}),:).ca_rate}];
end
figure("Position", [195,243,246,118])
plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylims,'Calcium Event Rate (Hz)');
%% Panel E-F: resp vs unresp spines turnover
resp_data_feature_not_retained = [];
resp_data_feature_retained = [];

unresp_data_feature_not_retained = [];
unresp_data_feature_retained = [];
feature = 'all_z_scored_trace';
resp_type = ['000'];
use_log = 0;
for i = 1:length(data)
    temp = data{i};
    temp.roi_fovs_mouse_cell = strcat(num2str(temp.all_roi_inds),'_', temp.all_fovs, '_', temp.all_mice_cells);

    %extract each session and then combine
    features = {'nonresp','resp',feature,'roi_fovs_mouse_cell', 'structure_type'};
    contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, "roi_fovs_mouse_cell");
    contra_ipsi_binoc.cib = arrayfun(@(x,y,z) [num2str(x);num2str(y);num2str(z)], contra_ipsi_binoc.('resp_contra'),...
        contra_ipsi_binoc.('resp_ipsi'),contra_ipsi_binoc.resp, 'UniformOutput',false );
    
    %concatenate trace across all 3 sessions
    concatenate_z_scored_trace = cellfun(@(x,y,z) [x;y;z], contra_ipsi_binoc.([feature,'_contra']),...
        contra_ipsi_binoc.([feature,'_ipsi']),contra_ipsi_binoc.(feature), 'UniformOutput',false);
    
    %smooth trace using exponential weighted filter 
    smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

    %count spine calcium event rate
    [pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), smoothed_z_scored_trace,'UniformOutput',false);
    if use_log
         ca_rate = cellfun(@(x,y) log10(length(x)*fs/length(y)), pks, concatenate_z_scored_trace);
         ylims = [-3.5,-1];
    else
         ca_rate = cellfun(@(x,y) (length(x)*fs/length(y)), pks, concatenate_z_scored_trace);
         ylims = [0,0.2];
    end
     ca_rate(isinf(ca_rate))=NaN;
    contra_ipsi_binoc.ca_rate = ca_rate;
    resp_cib = contra_ipsi_binoc(~strcmp(contra_ipsi_binoc.cib, resp_type),:);
    unresp_cib = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.cib, resp_type),:);
    
    resp_data_feature_not_retained = [resp_data_feature_not_retained, { resp_cib(ismember( resp_cib.structure_type, structure{i}),:).ca_rate}];
    resp_data_feature_retained = [resp_data_feature_retained, { resp_cib(~ismember( resp_cib.structure_type, structure{i}),:).ca_rate}];

    unresp_data_feature_not_retained = [unresp_data_feature_not_retained, { unresp_cib(ismember( unresp_cib.structure_type, structure{i}),:).ca_rate}];
    unresp_data_feature_retained = [unresp_data_feature_retained, { unresp_cib(~ismember( unresp_cib.structure_type, structure{i}),:).ca_rate}];

end

figure("Position",[345,394,290,190])
plot_bar_lost_added_retained(resp_data_feature_not_retained,resp_data_feature_retained, [0, 0.025],'Calcium Event Rate (Hz)');
figure("Position",[345,394,290,190])
plot_bar_lost_added_retained(unresp_data_feature_not_retained,unresp_data_feature_retained, [0, 0.025],'Calcium Event Rate (Hz)');

%% Panel D,G: load population data and get ca activity of spines by day

savepath = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files'; %update this path for where data will be stored
load(fullfile(savepath,'spine_properties_table.mat'));
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

%% Panel D: resp vs unresp overall calcium activity
contra_ipsi_binoc.cib = arrayfun(@(x,y,z) [num2str(x);num2str(y);num2str(z)], contra_ipsi_binoc.('resp_contra'),...
        contra_ipsi_binoc.('resp_ipsi'),contra_ipsi_binoc.resp, 'UniformOutput',false );
  
unresp = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.cib, ['000']),:);
resp = contra_ipsi_binoc(~strcmp(contra_ipsi_binoc.cib, ['000']),:);

figure
unresp = unresp.ca_rate;
unresp(isinf(unresp)) = NaN;
cdfplot(unresp); hold on;

resp = resp.ca_rate;
resp(isinf(unresp)) = NaN;
cdfplot(resp)
xlim([0,0.1]);

[p,a] = ranksum(unresp,resp)



%% Panel G: binoc spine calcium activity by day for spines that are responsive to at least binoc 
close all
days = {'D1', 'D5', 'D10'};
figure

eye_type = {'001', '101', '011', '111'};
all_temp = [];
all_eye_type = [];
all_day = [];

for i = 1:length(days)

    subplot(1, length(days),i)
    contra_ipsi_binoc_days = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.days, days{i}),:);
    
    
    for ii = 1:length(eye_type)
        temp = contra_ipsi_binoc_days(strcmp(contra_ipsi_binoc_days.cib, eye_type{ii}),:).ca_rate;
        
        ylims = [0,0.06];
        


        bar(ii, nanmean(temp)); hold on;
        errorbar(ii, nanmean(temp), nanstd(temp./length(temp)), '.k'); hold on;
                    title(days{i})
        ylim(ylims);
        xlim([0,5])
        if ii == length(eye_type)
             xticks(1:length(eye_type))
             xticklabels(eye_type)
        end

        all_temp = [all_temp; temp];
        all_day = [all_day; repmat(days(i), size(temp))];
        all_eye_type = [all_eye_type; repmat(eye_type(ii), size(temp))];
   
        
    end
end
%% Panel G stats
[p,tbl_p, stats] = anovan(all_temp, {all_day, all_eye_type}, 'model', 'interaction')
[results,~,~,gnames]=multcompare(stats, 'Dimension', [1,2])
tbl2 = array2table(results,"VariableNames", ...
["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A")=gnames(tbl2.("Group A"));
tbl2.("Group B")=gnames(tbl2.("Group B"))


