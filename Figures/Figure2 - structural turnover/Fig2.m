%% Figure 2 
% Panel B: Create bar chart show loss, addition and retention of spines
% from D1-D5 and D5-D10
% Panel C: Create flow chart of spines that were added, lost, retained across the 10 day period
% Panel D: loss, retention, addition of spines by response type 
% Panel E: binooc response type for retained spines and somas
% Panel F: eye-specific response type fore retained spines and somas
% Panel G: Fraction of binoc, contra, ipsi responsive spines on somas per
% day

%% Load data 
clear all
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

%% Panel B: Bar chart for spine loss, addition, retention
temp_structure.D1_D5_table = D1_D5_table;
temp_structure.D5_D10_table = D5_D10_table;
unique_mouse_cell = unique(union(D1_D5_table.all_mice_cells, D5_D10_table.all_mice_cells));
percents =get_percent_lost_formed_retained(temp_structure, unique_mouse_cell);

nansum(percents.all_spines)
all_percent_lost_formed_retained = percents.all_percent;
feature = all_percent_lost_formed_retained;
x_labels = {'lost', 'added', 'retained'};
y_labels = {'Fraction spines per cell'};
legend_labels = {'D1-D5', 'D5-D10'};
plot_bar_scatter_plot(feature, x_labels, y_labels, legend_labels)
title("Fraction lost, added, retained")
set(gca, 'FontSize',15);
ylim([0,0.8])
for i = 1:size(feature,2)
    d1_5 = feature(:,i,1);
    d5_d10 = feature(:,i,2);
    inds = any(isnan([d1_5, d5_d10]),2);
    
    [p,h] = ranksum(d1_5,d5_d10)
end


%% Panel C: Join D1-D5 and D5-D10 to get spines tracked across 10 days and create sandkey flow chart at http://sankeymatic.com
get_10day_spine_fates(D1_D5_table, D5_D10_table);


%% Panel D: Get fraction of contra, ipsi, contra+ipsi spines retained from D1 to D10 normalized by number of contra, ipsi, contra+ipsi spines 
close all
all_temps = {D1_D5_table, D5_D10_table};
get_fraction_retained_eye_preference_by_cell(all_temps)

%% Panel E-F (left): Load data and combine D1 to D5 and D5 to D10 table for retained spines
clear all
close all
chronic_path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files'; %update this path for where data will be stored
load(fullfile(chronic_path,"tracked_spines_properties_table.mat") );

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

%% Panel E (left): get eye specific identities of spines retained from D1 to D10 and create sandkey flow chart at http://sankeymatic.com
all_counts = get_10day_retained_spine_identity_contra_ipsi_only(retain_D1_D5, retain_D5_D10, mice_cells_D1_D10);
%% Panel F (left): get binoc response of spines retained from D1 to D10 and create sandkey flow chart at http://sankeymatic.com
all_counts = get_10day_retained_spine_identity_both_eyes(retain_D1_D5, retain_D5_D10, mice_cells_D1_D10);

%% Panel E-F (right): get eye-specific identities of neurons tracked from D1 to D10
path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files'; %update this path for where data will be stored
file = 'soma_properties_table.mat';
load(fullfile(path,file))
vars = {'resp', 'days', 'session', 'mouse_cell'}
D1_table = all_soma_stim_table(strcmp(all_soma_stim_table.days, 'D1'),vars);
D5_table = all_soma_stim_table(strcmp(all_soma_stim_table.days, 'D5'),vars);
D10_table = all_soma_stim_table(strcmp(all_soma_stim_table.days, 'D10'),vars);
unique_mouse_cell_D1_D5 = intersect(D1_table.mouse_cell,D5_table.mouse_cell);
unique_mouse_cell_D1_D10 = intersect(D10_table.mouse_cell,unique_mouse_cell_D1_D5);

data = {D1_table, D5_table, D10_table};
%% get soma eye specific identities from D1 to D10 and create sandkey flow chart at http://sankeymatic.com
all_counts = get_10day_soma_identity_contra_ipsi(data, unique_mouse_cell_D1_D10);
%% get soma binoc response from D1 to D10 and create sandkey flow chart at http://sankeymatic.com
all_counts = get_10day_soma_identity_both_eyes(data, unique_mouse_cell_D1_D10);
%% Panel G: fraction of spines that match the soma's pref 
% load data
path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files'; %update this path for where data will be stored
load(fullfile(path, "spine_properties_table.mat"), "all_stim_table");
all_stim_table.roi_fovs_mouse_cell_days = strcat(num2str(all_stim_table.all_roi_inds), '_',...
                                                all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);
temp = all_stim_table;
features = {'spine_mean_resp','resp','soma_resp','roi_fovs_mouse_cell_days', 'days', 'mouse_cell'};
join_key = 'roi_fovs_mouse_cell_days';
contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, join_key);

%create one hot encoding for spines and soma eye pref
contra_ipsi_binoc.ci_spine = strcat(num2str(contra_ipsi_binoc.resp_contra), num2str(contra_ipsi_binoc.resp_ipsi));
contra_ipsi_binoc.cib_spine = strcat(num2str(contra_ipsi_binoc.resp_contra), num2str(contra_ipsi_binoc.resp_ipsi),num2str(contra_ipsi_binoc.resp));

contra_ipsi_binoc.ci_soma = strcat(num2str(contra_ipsi_binoc.soma_resp_contra), num2str(contra_ipsi_binoc.soma_resp_ipsi));
contra_ipsi_binoc.cib_soma = strcat(num2str(contra_ipsi_binoc.soma_resp_contra), num2str(contra_ipsi_binoc.soma_resp_ipsi),num2str(contra_ipsi_binoc.soma_resp));


%% Fraction of Binoc Spines on Binoc Somas 
close all
figure("Position",[440,690,135,107], 'Name', 'Binoc')
contra_ipsi_binoc_resp = contra_ipsi_binoc(contra_ipsi_binoc.soma_resp==1,:);
days = {'D1', 'D5','D10'};
all_fraction_match = [];
all_days = [];
for i = 1:length(days)
    contra_ipsi_binoc_day = contra_ipsi_binoc_resp(strcmp(contra_ipsi_binoc_resp.days, days{i}),:);
    x = contra_ipsi_binoc_day.resp;
    g = findgroups( contra_ipsi_binoc_day.mouse_cell);
    fraction_match = splitapply(@mean, x, g);
    all_fraction_match = [all_fraction_match; fraction_match];
    all_days = [all_days; repmat(days(i), size(fraction_match))];
    bar(i, mean(fraction_match)); hold on;
       mean(fraction_match)
      std(fraction_match)/sqrt(length(fraction_match))

    errorbar(i,mean(fraction_match),std(fraction_match)/sqrt(length(fraction_match)), '.k' )
    scatter(i*ones(size(fraction_match)), fraction_match, '.k', 'jitter', 'on', 'jitterAmount', 0.05); 

end
xticks([1:3])
ylim([0,0.8])
xticklabels(days)
ylabel('Fraction spines/cell')

kruskalwallis(all_fraction_match, all_days)
%% Fraction of Contra Spines on Contra Somas 
close all
figure("Position",[440,690,135,107], 'Name', 'Contra')
contra_ipsi_binoc_resp = contra_ipsi_binoc(contra_ipsi_binoc.soma_resp_contra==1,:);
days = {'D1', 'D5','D10'};
all_fraction_match = [];
all_days = [];
for i = 1:length(days)
    contra_ipsi_binoc_day = contra_ipsi_binoc_resp(strcmp(contra_ipsi_binoc_resp.days, days{i}),:);
    x = contra_ipsi_binoc_day.resp_contra;
    g = findgroups( contra_ipsi_binoc_day.mouse_cell);
    fraction_match = splitapply(@mean, x, g);
    all_fraction_match = [all_fraction_match; fraction_match];
    all_days = [all_days; repmat(days(i), size(fraction_match))];
    bar(i, mean(fraction_match)); hold on;
      mean(fraction_match)
      std(fraction_match)/sqrt(length(fraction_match))
    errorbar(i,mean(fraction_match),std(fraction_match)/sqrt(length(fraction_match)), '.k' )
    scatter(i*ones(size(fraction_match)), fraction_match, '.k', 'jitter', 'on', 'jitterAmount', 0.05); 

end
xticks([1:3])
ylim([0,0.8])
xticklabels(days)
ylabel('Fraction spines/cell')
kruskalwallis(all_fraction_match, all_days)
%% Fraction of Ipsi Spines on Ipsi Somas 

close all
figure("Position",[440,690,135,107], 'Name', 'Ipsi')
contra_ipsi_binoc_resp = contra_ipsi_binoc(contra_ipsi_binoc.soma_resp_ipsi==1,:);
days = {'D1', 'D5','D10'};
all_fraction_match = [];
all_days = [];
for i = 1:length(days)
    contra_ipsi_binoc_day = contra_ipsi_binoc_resp(strcmp(contra_ipsi_binoc_resp.days, days{i}),:);
    x = contra_ipsi_binoc_day.resp_ipsi;
    g = findgroups( contra_ipsi_binoc_day.mouse_cell);
    fraction_match = splitapply(@mean, x, g);
    all_fraction_match = [all_fraction_match; fraction_match];
    all_days = [all_days; repmat(days(i), size(fraction_match))];
    bar(i, mean(fraction_match)); hold on;
        mean(fraction_match)
      std(fraction_match)/sqrt(length(fraction_match))

    errorbar(i,mean(fraction_match),std(fraction_match)/sqrt(length(fraction_match)), '.k' )
    scatter(i*ones(size(fraction_match)), fraction_match, '.k', 'jitter', 'on', 'jitterAmount', 0.05); 

end
xticks([1:3])
ylim([0,0.8])
xticklabels(days)
ylabel('Fraction spines/cell')
kruskalwallis(all_fraction_match, all_days)

%%
function [percents] = get_percent_lost_formed_retained(temp_structure, unique_mouse_cell)
    fn = fieldnames(temp_structure);
    all_percent_lost_formed_retained = NaN(length(unique_mouse_cell),3,numel(fn));
    all_spines = NaN(length(unique_mouse_cell),numel(fn));
    all_dends = NaN(length(unique_mouse_cell),numel(fn));
    for t = 1:numel(fn)
        temp = temp_structure.(fn{t});
        percent_lost_formed_retained = [];
        spines = [];
        dends = [];
        for m = 1:length(unique_mouse_cell) 
            %get spines from mouse
            filter_temp_table = temp(contains...
                (temp.all_mice_cells,unique_mouse_cell{m})&...
                contains(temp.session,"binoc"),:);
            num_spines = height(filter_temp_table);
            if num_spines == 0
                num_spines = NaN;
                
            end
            spines = [spines;num_spines];
            
            segments = filter_temp_table.all_fovs;
            dend_name = split(segments, '_');
            dend_name = dend_name(:,1);
            num_dends = length(unique(dend_name));
            if num_dends == 0
                num_dends = NaN;
            end
            dends = [dends; num_dends];
            [lost, formed, retained] = get_percents(filter_temp_table);
            %Find total percent retained, eliminated, formed
            percent_lost_formed_retained = [percent_lost_formed_retained;...
                [lost, formed, retained]];              
        end
        all_spines(:,t) = spines;
        all_percent_lost_formed_retained(:,:,t) = percent_lost_formed_retained;
        all_dends(:,t) = dends;  
    end
    percents.all_dends = all_dends;
    percents.all_spines = all_spines;
    percents.all_percent = all_percent_lost_formed_retained;
end
%%
function plot_bar_scatter_plot(feature, x_labels, y_labels, legend_labels)
figure
mean_percent = nanmean(feature,1);
mean_percent = squeeze(mean_percent);

bar(mean_percent); hold on
ngroups = size(mean_percent, 1);
nbars = size(mean_percent, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    n = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    y_by_d = feature(:,:,i);
    %y_by_d = reshape(y_by_d, size(y_by_d,2), size(y_by_d,1));
    hs = scatter(ones(length(y_by_d),1)*n,y_by_d,25,'k', "filled"); hold on;
    if i == 1
        n2 = (1:ngroups) - groupwidth/2 + (2*(i+1)-1) * groupwidth / (2*nbars);
        y_by_d2 = feature(:,:,i+1);
        for ii = 1:size(y_by_d2,2)
            for iii = 1:size(y_by_d2,1)
                plot([ones(length(y_by_d(iii,ii)),1)*n(ii),...
                    ones(length(y_by_d(iii,ii)),1)*n2(ii)],...
                    [y_by_d(iii,ii),y_by_d2(iii,ii)], 'k');
            end
        end
    end

end
clear nbars ngroups mean_percent groupwidth hs i ii iii n n2 y_by_d2 y_by_d
ylabel(y_labels);
legend(legend_labels);
xticks([1:length(x_labels)]);
xticklabels(x_labels)
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

