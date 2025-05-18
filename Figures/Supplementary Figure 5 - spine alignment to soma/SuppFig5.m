%% Panel A&B: load data 
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
%% Panel A: Plot fract of active trials to soma's pref direction among lost, added, and retained spines
close all
feature = 'all_active_spine_trials_smooth';
feature2 = 'fract_active_trials_pref_dir';
session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,5];

for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        temp = temp(temp.resp<=0,:);        
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
    %plot_violinchart_lost_added_retained(data_feature_not_retained,data_feature_retained, '% active trials', [0,100]);
    %plot_boxchart_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature])
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature]);

end
%% Panel B: fract of trials active to soma's orthogonal direction among lost, added, retained spines
close all
feature = 'all_active_spine_trials_smooth';
feature2 = 'fract_active_trials_pref_dir';
session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,20];
num_dir = 8;
for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        temp = temp(temp.resp<=0,:);        
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
    %plot_boxchart_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature])


end