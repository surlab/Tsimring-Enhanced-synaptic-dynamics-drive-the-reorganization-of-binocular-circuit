%% Figure 5 
% Panel A: Plot the polarplots for lost, added and retained spine examples for two time points
% Panel B: Plot OSI for lost, retained, added spines (binoc viewing)
% Panel C: plot delta OSI for retained spines (binoc viewing)
% Panel D: plot delta Ori Pref for retained spines (binoc viewing)
% Panel E: plot orientation offset to soma between lost, retained and added spines
% (binoc viewing)
% Panel F: CDF of orientation offset to soma for all spines from D1 to D10
%% Panel A: Plot the polarplots for lost, added and retained spine examples for two time points
plot_ori_tuning_example_soma_spine
%% Panel B: load data
clear all
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

%% Panel B: Plot OSI for binoc session lost, added, retained spines (binoc viewing)
% only include spines that are visually responsive to binoc viewing
close all
feature = 'all_OSI_vector';

session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,1];

for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        temp = temp(temp.resp>0,:);        
 
        data_feature_not_retained = [data_feature_not_retained, {temp(ismember(temp.structure_type, structure{i}),:).(feature)}];
        data_feature_retained = [data_feature_retained, {temp(~ismember(temp.structure_type, structure{i}),:).(feature)}];
    end
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature]);

end



%% Panel E: Plot ori offset among lost, added, and retained spines vector based approach
close all
feature = 'all_ori_alignment_vector';

session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,90];

for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        %temp = temp(temp.neighbor_spine_distance<=5&temp.resp>0,:);        
        temp = temp(temp.resp>0&temp.soma_resp>0&temp.all_OSI_vector>0.3,:);    
 
        data_feature_not_retained = [data_feature_not_retained, {temp(ismember(temp.structure_type, structure{i}),:).(feature)}];
        data_feature_retained = [data_feature_retained, {temp(~ismember(temp.structure_type, structure{i}),:).(feature)}];
    end
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature]);

end
%% Panel C-D: Load bootstrapped data 
clear all
chronic_path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
load(fullfile(chronic_path, 'tracked_spines_ori_bootstrap.mat'),'save_temp');

%% Panel C-D: plot bootstrapped OSI and orientation preference
close all
sessions = {'binoc', 'contra', 'ipsi'};
days = {'D1', 'D5', 'D10'};
colors = {'g', 'r', 'b'}
close all
g = figure("Name",['Orientation Difference'], "Position",[440,479,437,318]);
h = figure("Name",['Direction Difference'],"Position",[440,479,437,318]);
a = figure("Name",['OSI Difference'],"Position",[440,479,437,318]);
b = figure("Name",['DSI Difference'],"Position",[440,479,437,318]);

for i = 1:1
    [temp] =compare_boot_strap_distribution_ori_offset_spines(g,h,a,b,save_temp, days, length(sessions), i, sessions{i});
    dist_temp.(sessions{i}) = temp;
end

%% Panel F: Load data
clear all
save_path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
load(fullfile(save_path, "spine_properties_table.mat"), "all_stim_table");

% merge table into contra and ipsi session per spine
all_stim_table.roi_fovs_mouse_cell_days = strcat(num2str(all_stim_table.all_roi_inds), '_', ...
                                                all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);

%% Panel F: Plot orientation offset for binoc session
sessions = {'binoc'};
resp_spine_stim = all_stim_table(all_stim_table.resp>0&all_stim_table.all_OSI_vector>0.3,:);
resp_spine_stim = resp_spine_stim(resp_spine_stim.soma_resp>0,:);
offset_type = 'all_ori_alignment_vector';
figure
%bins = [-1:0.2:1];
% all visually responsive spines
for i = 1:length(sessions)
      all_ori = [];
      all_days = [];
      session_spine = resp_spine_stim(strcmp([resp_spine_stim.session{:}], sessions(i))',:);
      days = {'D1','D5', 'D10'};
      for ii = 1:length(days)
          subplot(2,3,i);
          session_spine_day = session_spine(strcmp(session_spine.days,days{ii}),:);
          ori = session_spine_day.(offset_type);
          %swarmchart(ii*ones(size(ori)), ori); hold on;
          cdfplot(ori); hold on;
          disp(['Mean ori: ', num2str(mean(ori))])
          disp(['SEM ori: ', num2str(std(ori)/sqrt(length(ori)))])
          all_ori= [all_ori;ori];
          all_days = [all_days; repmat(days(ii),size(ori))];

      end

      if i == 1
          
          xlabel('Orientation Offset (soma-spine)')
          
      end
      title(sessions(i))      
end
[p,a, stats] = kruskalwallis(all_ori, all_days);
[results,~,~,gnames]=multcompare(stats, 'Dimension', [1,2])
tbl2 = array2table(results,"VariableNames", ...
["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A")=gnames(tbl2.("Group A"));
tbl2.("Group B")=gnames(tbl2.("Group B"))




