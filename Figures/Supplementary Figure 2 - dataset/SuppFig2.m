%% Supplementary Figure 2
% Panel A: plot the number of dendrites per cell by day
% Panel B: Plot the number of spines per cell by day
% Panel C: Plot the spine density per dendritic segment by day



%% load spine table across neurons and days and sessions
clear all 
close all
path_dir = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
file = 'spine_properties_table.mat'
filepath = fullfile(path_dir,file);
load(filepath);

all_stim_table.roi_fovs_mouse_cell_days = strcat(num2str(all_stim_table.all_roi_inds), '_', ...
                                                all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);
%%
temp = all_stim_table;
contra = temp(strcmp([temp.session{:}], "contra"),:);
contra = contra(:, {'resp','spine_mean_resp','soma_mean_resp','roi_fovs_mouse_cell_days','days'});
ipsi = temp(strcmp([temp.session{:}], "ipsi"),:);
ipsi = ipsi(:, {'resp','spine_mean_resp','soma_mean_resp','roi_fovs_mouse_cell_days','days'});
binoc = temp(strcmp([temp.session{:}], "binoc"),:);
binoc = binoc(:, {'resp','spine_mean_resp','soma_mean_resp','roi_fovs_mouse_cell_days','days', 'mouse_cell', 'all_fovs'});

contra_ipsi = outerjoin(contra, ipsi, "Keys","roi_fovs_mouse_cell_days", "MergeKeys",true);
contra_ipsi_binoc = outerjoin(contra_ipsi, binoc, "Keys","roi_fovs_mouse_cell_days", "MergeKeys",true);

%remove rows with NaN values 
contra_ipsi_binoc=contra_ipsi_binoc(~any(ismissing(contra_ipsi_binoc),2),:);

%% number of spines, dendrites, spine by fov per cell per day
all_day_table = [];
all_spines_table = [];
all_spine_density_table = [];
all_dend_table = [];
all_mouse_table = [];

temp = contra_ipsi_binoc;
temp = all_stim_table(strcmp([all_stim_table.session{:}], "binoc"), :);
unique_mouse_cell = unique(temp.mouse_cell);
unique_day = {'D1', 'D5', 'D10'};
for i = 1:length(unique_mouse_cell)
    dends_table = NaN(1,3);
    spines_table = NaN(1,3);
    spine_density_table = NaN(1,3);
    mouse_cell = unique_mouse_cell{i};
    %days_table = zeros(3,1);
    vals = temp(ismember(temp.mouse_cell,mouse_cell),:);
    for d = 1:length(unique_day)
    
        rois = vals(ismember(vals.days,unique_day{d}),:); 
        if ~isempty(rois)
            
            dends = table2array(rois(:,"all_fovs"));
            temp_dends = cellfun(@(x) strsplit(x, '_'), dends, 'UniformOutput', false);
            temp_dends = cellfun(@(x) x{:,1}, temp_dends, 'UniformOutput', false); 
            temp_dends = string(temp_dends);
                
            unique_dends = unique(temp_dends);
            spine_density = [];
            for dd = 1:length(unique_dends)
                spines_by_dend = rois(contains...
                    (rois.all_fovs,unique_dends{dd}),:);
                num_spines = height(spines_by_dend);
                spine_density = [spine_density,num_spines];
            end
    
            dends_table(d) = length(unique_dends);
        
            spines_table(d) = length(temp_dends);
            spine_density_table(d) = mean(spine_density);
        %             OSI_table = [OSI_table, table2array(rois_binoc(1,"all_OSI"))];
        %             mean_amp_table = [mean_amp_table, table2array(rois_binoc(1,"all_mean_Z_day"))];   
            %days_table(ind) = day};
            mouse_table = {mouse_cell};
        end
end
all_dend_table = [all_dend_table; dends_table];
all_spines_table = [all_spines_table; spines_table];
all_spine_density_table = [all_spine_density_table; spine_density_table];
%all_day_table = [all_day_table; days_table];
all_mouse_table = [all_mouse_table; mouse_table];
end

%% Panel A 

feature = all_spines_table;
title_name = 'Spines per Cell';
y_label = 'Count';
x_label = 'Days';
x_ticks = {'D1', 'D5', 'D10'};

num_cells = size(feature,1);
%colors = rand(num_cells,3);
figure("Position",[500,352,354,445])
line_plot_cells(x_ticks,feature,num_cells,all_mouse_table,x_label, y_label, title_name);
nansum(feature)
nanmean(feature)
nanstd(feature)
day1 = sum(~isnan(feature(:,1)))
day2 = sum(~isnan(feature(:,2)))
day3 = sum(~isnan(feature(:,3)))
%% Panel B: Plot dendrites per cell per day
feature = all_dend_table;
title_name = 'Dendritic Segments per Cell';
y_label = 'Count';
x_label = 'Days';
x_ticks = {'D1', 'D5', 'D10'};

figure("Position",[500,352,354,445])
line_plot_cells(x_ticks,feature,num_cells,all_mouse_table,x_label, y_label, title_name);
nansum(feature)
nanmean(feature)
nanstd(feature)

%% Panel C: plot average spines per dendrite per cell per day
feature = all_spine_density_table;
title_name = 'Spines per Dendritic Segment';
y_label = 'Count';
x_label = 'Days';
x_ticks = {'D1', 'D5', 'D10'};
num_cells = size(feature,1);

figure("Position",[500,352,354,445])
line_plot_cells(x_ticks,feature,num_cells,all_mouse_table, x_label, y_label, title_name);
nansum(feature)
nanmean(feature)
nanstd(feature)
size(sum(~isnan(feature(:,1))))
%% Line plot of features by cell by day
function line_plot_cells(x_ticks,feature, num_cells, all_mouse_table,x_label, y_label, title_name)
x = 1:length(x_ticks);    
for c = 1:num_cells
        if sum(isnan(feature(c,:)))==0
            plot(x,feature(c,:), '-o','Color','g', 'LineWidth',3,'MarkerSize',8, 'MarkerFaceColor','g'); hold on;
        elseif sum(~isnan(feature(c,1:2)))==2
             plot(x,feature(c,:),'-o','Color','r', 'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','r'); hold on;
        elseif sum(~isnan(feature(c,2:3)))==2
            plot(x,feature(c,:), '-o','Color','b', 'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','b'); hold on;
        elseif sum(~isnan(feature(c,[1,3])))==2
            plot(x,feature(c,:), 'o','Color','m', 'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','m'); hold on;
        else
            plot(x,feature(c,:), '-o','Color','k', 'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','k'); hold on;
        end
    
    end
    xticks(x);
    xticklabels(x_ticks);
    xtickangle(45);
    xlabel(x_label);
    ylabel(y_label);
    title(title_name);
    xlim([1-0.5, x(end)+0.5]);
    legend(all_mouse_table, "Location","bestoutside");
    set(gca,"FontSize",10)
end