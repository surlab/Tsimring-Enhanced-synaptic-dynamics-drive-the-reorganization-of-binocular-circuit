%% Figure 4
% Panel A: plot the % active trials to soma's preferred and orthogonal directions for lost, added, retained spine example
% Panel B: plot the % active trials to soma's preferred direction for lost, added, retained spines
% Panel C: plot the % active trials to soma's orthogonal directions for lost, added, retained spines
% Panel D: plot the % active trials to soma's preferred direction for the same retained spines
% Panel E: plot the % active trials to soma's orthogonal directions for the same retained spines

%% Panel A
plot_active_trials_soma_spine
%% Panel B&E: load data 
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
%% Panel B: Plot fract of active trials to soma's pref direction among lost, added, and retained spines during binoc viewing
% both soma and spine are responsive to binoc viewing
close all
feature = 'all_active_spine_trials_smooth';
feature2 = 'fract_active_trials_pref_dir';
session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,50];

for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        temp = temp(temp.resp>0,:);        
        temp = temp(temp.soma_resp>0,:); 
        soma_mean_resp = temp.soma_mean_resp;
        peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);
 
        spine_trials = temp.(feature);
        inc_trials = temp.all_included_trial;
        sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
        spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
        fract_active_trials_pref_dir = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,peak_dir, sum_inc_trials);
        temp.(feature2) = fract_active_trials_pref_dir;
 
        data_feature_not_retained = [data_feature_not_retained, {temp(ismember(temp.structure_type, structure{i}),:).(feature2)}];
        data_feature_retained = [data_feature_retained, {temp(~ismember(temp.structure_type, structure{i}),:).(feature2)}];
    end
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature]);

end
%% Panel C: fract of trials active to soma's orthogonal direction among lost, added, retained spines during binoc viewing
% both soma and spine are responsive to binoc viewing
close all
feature = 'all_active_spine_trials_smooth';
feature2 = 'fract_active_trials_pref_dir';
session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,10];
num_dir = 8;
for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        temp = temp(temp.resp>0,:);        
        temp = temp(temp.soma_resp>0,:); 
        soma_mean_resp = temp.soma_mean_resp;
        peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);
        ortho_dir1 = cellfun(@(x) mod(x+2,num_dir), peak_dir, 'UniformOutput',false);
        ortho_dir2 = cellfun(@(x)(x-2), peak_dir, 'UniformOutput',false);
        ortho_dir1([ortho_dir1{:}]==0)={8};
        ortho_dir2([ortho_dir2{:}]<0)=cellfun(@(x) num_dir+x, ortho_dir2([ortho_dir2{:}]<0), 'UniformOutput',false);
        ortho_dir2([ortho_dir2{:}]==0)={8};
        spine_trials = temp.(feature);
        inc_trials = temp.all_included_trial;
        sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
        spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
        fract_active_trials_ortho1 = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,ortho_dir1, sum_inc_trials);
        fract_active_trials_ortho2 = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,ortho_dir2, sum_inc_trials);
        mean_fract_active_trials = mean([fract_active_trials_ortho1,fract_active_trials_ortho2],2);
        temp.(feature2) = mean_fract_active_trials;
 
        data_feature_not_retained = [data_feature_not_retained, {temp(ismember(temp.structure_type, structure{i}),:).(feature2)}];
        data_feature_retained = [data_feature_retained, {temp(~ismember(temp.structure_type, structure{i}),:).(feature2)}];
    end
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature]);

end

%% Panel D-E: Load data for retained spines 
close all
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

%% Panel D: plot fraction of trials spine was active for soma's preferred and orthogonal direction during binoc viewing
% retained spines must be visually responsive for both time points
all_temps = {retain_D1_D5, retain_D5_D10};
feature = 'all_active_spine_trials_smooth';
unique_sessions = {'binoc'};
days = {'D1', 'D5', 'D10'};
g = figure("Position",[440,460,170,337]); % preferred direction
h = figure("Position",[440,460,170,337]); % orthogonal direction

count = 1;

for i = 1:length(all_temps)
    count2 = 1;
    for s = 1:length(unique_sessions)
        temp = all_temps{i};
        temp = temp(strcmp(temp.session,unique_sessions{s}), :);
        temp = splitvars(temp);

        d1 = days{i};
        d5 = days{i+1};
        
        %temp = temp(temp.([d1, '_soma_resp'])>0 &temp.([d5, '_soma_resp'])>0,:);
        temp = temp(temp.([d1, '_soma_resp'])>0&temp.([d5, '_soma_resp'])>0&...
            temp.([d1, '_resp'])>0&temp.([d5, '_resp'])>0,:);

        soma_mean_resp = temp.([d1, '_soma_mean_resp']);
        peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);
        spine_trials = temp.([d1, '_',feature]);
        inc_trials = temp.([d1, '_all_included_trial']);
        sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
        spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
        
        fract_active_trials_pref_dir_d1 = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,peak_dir, sum_inc_trials);
        ortho_dir1 = cellfun(@(x) mod(x+2,num_dir), peak_dir, 'UniformOutput',false);
        ortho_dir2 = cellfun(@(x)(x-2), peak_dir, 'UniformOutput',false);
        ortho_dir1([ortho_dir1{:}]==0)={8};
        ortho_dir2([ortho_dir2{:}]<0)=cellfun(@(x) num_dir+x, ortho_dir2([ortho_dir2{:}]<0), 'UniformOutput',false);
        ortho_dir2([ortho_dir2{:}]==0)={8};

        fract_active_trials_ortho1 = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,ortho_dir1, sum_inc_trials);
        fract_active_trials_ortho2 = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,ortho_dir2, sum_inc_trials);
        fract_active_trials_ortho_d1 = mean([ fract_active_trials_ortho1, fract_active_trials_ortho2 ],2);


    
        soma_mean_resp =temp.([d5, '_soma_mean_resp']);
        peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);
        spine_trials = temp.([d5, '_',feature]);
        inc_trials = temp.([d5, '_all_included_trial']);
        sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
        spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
        fract_active_trials_pref_dir_d5 = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,peak_dir, sum_inc_trials);
      
        ortho_dir1 = cellfun(@(x) mod(x+2,num_dir), peak_dir, 'UniformOutput',false);
        ortho_dir2 = cellfun(@(x)(x-2), peak_dir, 'UniformOutput',false);
        ortho_dir1([ortho_dir1{:}]==0)={8};
        ortho_dir2([ortho_dir2{:}]<0)=cellfun(@(x) num_dir+x, ortho_dir2([ortho_dir2{:}]<0), 'UniformOutput',false);
        ortho_dir2([ortho_dir2{:}]==0)={8};       
        
        fract_active_trials_ortho1 = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,ortho_dir1, sum_inc_trials);
        fract_active_trials_ortho2 = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))*100/z(y), spine_trials,ortho_dir2, sum_inc_trials);
        fract_active_trials_ortho_d5 = mean([ fract_active_trials_ortho1, fract_active_trials_ortho2 ],2);

        figure(g)
        subplot(2,length(unique_sessions), count2)
        for ii = 1:length(fract_active_trials_pref_dir_d5)
            plot([count, count+1],[fract_active_trials_pref_dir_d1(ii), fract_active_trials_pref_dir_d5(ii)], 'k'); hold on;
            scatter([count, count+1],[fract_active_trials_pref_dir_d1(ii), fract_active_trials_pref_dir_d5(ii)], '.k'); hold on;
        end
        [p,a]= signrank(fract_active_trials_pref_dir_d1,fract_active_trials_pref_dir_d5)
        plot_sig(count + 0.5, 100, p);
        plot([count, count+1],[mean(fract_active_trials_pref_dir_d1), mean(fract_active_trials_pref_dir_d5)],'r', 'LineWidth',2);
        ylabel(['% active trials (soma'' pref dir)'])
        title(unique_sessions{s})
        xticks([1,2,3,4])
        xlim([0,5])
        xticklabels({days{1},days{2},days{2},days{3} })
        ylim([0,100])

        
        figure(h)
        subplot(2,length(unique_sessions), count2)
        for ii = 1:length(fract_active_trials_pref_dir_d5)
            plot([count, count+1],[fract_active_trials_ortho_d1(ii), fract_active_trials_ortho_d5(ii)], 'k'); hold on;
            scatter([count, count+1],[fract_active_trials_ortho_d1(ii), fract_active_trials_ortho_d5(ii)], '.k'); hold on;
        end
        [p,a]= signrank(fract_active_trials_ortho_d1,fract_active_trials_ortho_d5)
        plot_sig(count + 0.5, 40, p);
        plot([count, count+1],[mean(fract_active_trials_ortho_d1), mean(fract_active_trials_ortho_d5)],'r', 'LineWidth',2);
        ylabel(['% active trials (soma'' ortho dirs)'])
        title(unique_sessions{s})
        xticks([1,2,3,4])
        xlim([0,5])
        xticklabels({days{1},days{2},days{2},days{3} })
        ylim([0,40])
        

        count2 = count2 + 1;
    end
    count = count + 2;
    
end
figure(g)
title('Preferred direction')
figure(h)
title('Orthogonal directions')