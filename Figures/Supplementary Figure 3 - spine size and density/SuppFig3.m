%% Supplementary Figure  3 
% % Panel A: Spine density over development (tracked neurons)
% Panel B: Spine Area over development
% Panel C: Spine Area Variance over development
% Panel D: Pie chart of retained, lost, added spines by response type 

%% Panel A spine density over development 
% use all the mice (even the untracked neurons to look at spine density,
% could be convinced otherwise 02/20/24)

% load data
path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files'; %update path
distance_path = fullfile(path, 'Chronic Imaging/FOV_alignment/'); %update path
load(fullfile(savepath, "spine_properties_table.mat"), "all_stim_table");
all_stim_table.fovs_mouse_cell_days = strcat(all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);
%% Plot spine density
% get spine density (using max distance between spines 02/20/24 -> might
% need to change to actual dendritic distance length)
temp = all_stim_table(strcmp([all_stim_table.session{:}], "binoc"), :);
unique_days = {'D1', 'D5', 'D10'};

get_resp = 0;
resp = 'resp';
[spine_density, unique_mouse_table] = plot_spine_density(distance_path, temp, unique_days, get_resp,1, resp);

close all
x = [1:size(spine_density,2)]; % number of days

% plot spine spine density
figure
for i = 1:size(spine_density ,1)
     plot(x,spine_density(i,:), '-o','Color','k', 'LineWidth',1,'MarkerSize',4,'MarkerFaceColor','k'); hold on;
end
xticklabels(unique_days)
xticks(x)
ylim([0,1.5])
xlim([0,x(end)+1])
ylabel('Spine density (per 10 um)')

days = [ones(size(spine_density,1),1); 2*ones(size(spine_density,1),1); 3*ones(size(spine_density,1),1)];
 x = spine_density(:);
x_nonan = x(~isnan(x));
days = days(~isnan(x));
[p, tbl,stats] = kruskalwallis(x_nonan, days)

%% Panel B-C: mean spine area and variance over development 

% load plot spine size
close all
temp = all_stim_table(strcmp([all_stim_table.session{:}], "binoc"), :);
unique_days = {'D1', 'D5', 'D10'};
f = figure("Name",'Spine Size by Day (pooled across neurons)',"Position",[353,483,187,148]);
g = figure("Name",'Mean Spine Size by Neuron by Day',"Position",[353,483,187,148]);
h = figure("Name",'Spine Size Variance by Neuron by Day', "Position",[353,483,187,148]);
plot_spine_size(temp, unique_days,f,g,h)
%% Load data 
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
%% Panel D: spine turnover by eye-specific category
g = figure("Position", [1000,1033,245.6666666666667,305])
h = figure("Position", [1000,1033,245.6666666666667,305])
data = {D1_D5_table, D5_D10_table};
for i = 1:length(data)
    if i == 1
        D1 = 'D1';
        D5 = 'D5';
    else
        D1 = 'D5';
        D5 = 'D10';
    end
    temp = data{i};
    
    lost = splitvars(temp(contains(temp.structure_type, 'lost'), {D1, 'all_mice_cells', 'all_fovs', 'session', 'structure_type'}));
    added = splitvars(temp(contains(temp.structure_type, 'formed'),{D5, 'all_mice_cells', 'all_fovs', 'session','structure_type'}));
    retained = splitvars(temp(contains(temp.structure_type, 'retained'),{D1, 'all_mice_cells', 'all_fovs', 'session','structure_type'}));
    temp = [lost;added; retained];
    temp.roi_fovs_mouse_cell = strcat(num2str(temp.all_roi_inds),'_', temp.all_fovs, '_', temp.all_mice_cells, '_', temp.structure_type);
    %extract each session and then combine
    features = {'resp','roi_fovs_mouse_cell', 'structure_type'};
    join_key = 'roi_fovs_mouse_cell';
    contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, join_key);
    
    contra_ipsi_binoc.ci = arrayfun(@(x,y) [num2str(x);num2str(y)], contra_ipsi_binoc.('resp_contra'),...
        contra_ipsi_binoc.('resp_ipsi'), 'UniformOutput',false );
   contra_ipsi_binoc.cib = arrayfun(@(x,y,z) [num2str(x);num2str(y);num2str(z)], contra_ipsi_binoc.('resp_contra'),...
        contra_ipsi_binoc.('resp_ipsi'),contra_ipsi_binoc.('resp'), 'UniformOutput',false );


    %plot pie chart for each eye-specific category
    resp_type = {'10', '01', '11'};
    figure(g)
    for ii = 1:length(resp_type)
        
        [l,f,r] = get_percents(contra_ipsi_binoc(strcmp(contra_ipsi_binoc.ci, resp_type{ii}),:));
        subplot(3,2,2*(ii-1) + i)
        pieHandle = pie([l,f,r]);
        text_handles = pieHandle(2:2:end);


        for iHandle = 1:length(text_handles)
            text_handles(iHandle).Position = 0.2*text_handles(iHandle).Position;
        end
        if ii == 1
            title([D1, ' to ', D5])
        end
    end
    figure(h)
    resp_type = {'000'};
    
    [l,f,r] = get_percents(contra_ipsi_binoc(strcmp(contra_ipsi_binoc.cib, resp_type),:));
    subplot(3,2,i)
    pie([l,f,r]);
    pieHandle = pie([l,f,r]);
    text_handles = pieHandle(2:2:end);


    for iHandle = 1:length(text_handles)
        text_handles(iHandle).Position = 0.2*text_handles(iHandle).Position;
    end
    title([D1, ' to ', D5])
    [l,f,r] = get_percents(contra_ipsi_binoc(contra_ipsi_binoc.resp==1,:));
    subplot(3,2,i+2)
    pie([l,f,r]);
    pieHandle = pie([l,f,r]);
    text_handles = pieHandle(2:2:end);


    for iHandle = 1:length(text_handles)
        text_handles(iHandle).Position = 0.2*text_handles(iHandle).Position;
    end
       
            
        
   
end

%% get percents lost, added, retained
function [lost, formed, retained] = get_percents(filter_temp_table)
    % retained spine properties
    retain_temp_table = filter_temp_table(contains...
        (filter_temp_table.structure_type, "retained"),:);
    retained= height(retain_temp_table)/height(filter_temp_table);        
                               
    %eliminated spine properties
    eliminated_temp_table = filter_temp_table(contains...
        (filter_temp_table.structure_type, "lost"),:);
    lost= height(eliminated_temp_table)/height(filter_temp_table);
      
    %formed spine properties   
    formed_temp_table = filter_temp_table(contains...
        (filter_temp_table.structure_type, "formed"),:);
    formed = height(formed_temp_table)/height(filter_temp_table);
end

