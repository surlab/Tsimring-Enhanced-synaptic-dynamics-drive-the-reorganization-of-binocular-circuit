%% Supp Fig Figure 6 
% Panel A: show example FOV with tuning properties for spine and soma
% Panel B: Plot DSI for lost, retained, added spines
% Panel C: plot delta OSI, DSI, Ori Pref 
% Panel F: Bar plot of the mean change in DSI grouped by timepoint and
% session
%% Panel C: Load bootstrapped data 
clear all
chronic_path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
load(fullfile(chronic_path, 'tracked_spines_ori_bootstrap.mat'),'save_temp');

%% Panel C: plot bootstrapped OSI, DSI, and orientation preference
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
%%
all_temp = [];
for i = 1:1
    temp = dist_temp.(sessions{i});
    unique_days = unique(temp.all_days);
    for ii = 1:length(unique_days)
        %[p,a] = ttest(temp(strcmp(temp.all_days, unique_days{ii}),:).all_osi_difference);
        %disp(['P-value for ', unique_days{ii}, ' and session ',sessions{i}, ': ', num2str(a) ]);
        disp(['Mean for ', unique_days{ii}, ' and session ',sessions{i}, ': ', num2str(nanmean(temp(strcmp(temp.all_days, unique_days{ii}),:).all_ori_difference)) ])
        disp(['SEM for ', unique_days{ii}, ' and session ',sessions{i}, ': ', num2str(helper.nansem(temp(strcmp(temp.all_days, unique_days{ii}),:).all_ori_difference,1)) ])

    end
     all_temp = [all_temp; temp];
end
[p,stats,tbl] = anovan(all_temp.all_dsi_difference, {all_temp.all_days, all_temp.all_sessions}, 'model', 'interaction')
multcompare(tbl,'Dimension',[1,2])



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



%% Panel C: Plot DSI for binoc session lost, added, retained spines
close all
feature = 'all_DSI_vector';

session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,0.6];

for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        temp = temp(temp.resp>0,:);        
        %temp = temp(temp.neighbor_spine_resp>0,:);    
 
        data_feature_not_retained = [data_feature_not_retained, {temp(ismember(temp.structure_type, structure{i}),:).(feature)}];
        data_feature_retained = [data_feature_retained, {temp(~ismember(temp.structure_type, structure{i}),:).(feature)}];
    end
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature]);

end

%% Panel E: Plot direction offset among lost, added, and retained spines vector based approach
close all

feature = 'Dir offset (soma, spine)'
session = {'binoc'};
data = {splitvars(retain_vs_lost_D1_D5),splitvars(retain_vs_formed_D1_D5),splitvars(retain_vs_lost_D5_D10),splitvars(retain_vs_formed_D5_D10)};

structure = {"lost",  "formed", "lost", "formed"};
ylim_axis = [0,180];

for s = 1:length(session)
    figure("Position",[345,394,290,190])
    data_feature_not_retained = [];
    data_feature_retained = [];
    for i = 1:length(data)
        temp = data{i};
        temp = temp(ismember(temp.session, session{s}),:);
        %temp = temp(temp.neighbor_spine_distance<=5&temp.resp>0,:);        
        temp = temp(temp.resp>0&temp.soma_resp>0&temp.all_DSI_vector>0.3,:);    
       
        spine = temp(ismember(temp.structure_type, structure{i}),:).all_Dir_pref_vector;
        soma=  temp(ismember(temp.structure_type, structure{i}),:).soma_dir_pref_vector;
        offset = arrayfun(@(x,y) get_dir_offset(x,y), soma,spine);
        data_feature_not_retained = [data_feature_not_retained, {offset}];
        
        spine = temp(~ismember(temp.structure_type, structure{i}),:).all_Dir_pref_vector;
        soma=  temp(~ismember(temp.structure_type, structure{i}),:).soma_dir_pref_vector;
        offset = arrayfun(@(x,y) get_dir_offset(x,y), soma,spine);
        data_feature_retained = [data_feature_retained, {offset}];

        
    end
    plot_bar_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,[session{s}, ' ' feature]);

end

%% Panel F
clear all
save_path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
load(fullfile(save_path, "spine_properties_table.mat"), "all_stim_table");


% merge table into contra and ipsi session per spine
all_stim_table.roi_fovs_mouse_cell_days = strcat(num2str(all_stim_table.all_roi_inds), '_', ...
                                                all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);



%% Panel F: Plot direction offset for binoc session
sessions = {'binoc'};
resp_spine_stim = all_stim_table(all_stim_table.resp>0&all_stim_table.all_DSI_vector>0.3,:);
resp_spine_stim = resp_spine_stim(resp_spine_stim.soma_resp>0,:);

g = figure
%bins = [-1:0.2:1];
% all visually responsive spines
for i = 1:length(sessions)
      all_dir = [];
      all_days = [];
      session_spine = resp_spine_stim(strcmp([resp_spine_stim.session{:}], sessions(i))',:);
      days = {'D1','D5', 'D10'};
      for ii = 1:length(days)
          subplot(2,3,i);
          session_spine_day = session_spine(strcmp(session_spine.days,days{ii}),:);
          spine = session_spine_day.all_Dir_pref_vector;
          soma=  session_spine_day.soma_dir_pref_vector;
          offset = arrayfun(@(x,y) get_dir_offset(x,y), soma,spine);

          cdfplot(offset); hold on;
          temp_plot = offset;
          bar(ii, mean(temp_plot )); hold on;
          errorbar(ii,mean(temp_plot ), std(temp_plot )/sqrt(length(temp_plot)), '.k');

       disp(['Mean dir: ', num2str(mean(offset))])
          disp(['SEM dir: ', num2str(std(offset)/sqrt(length(offset)))])
          all_dir= [all_dir;offset];
          all_days = [all_days; repmat(days(ii),size(offset))];

      end

      if i == 1
          
          xlabel('Orientation Offset (soma-spine)')
          
      end
      title(sessions(i))      
end
[p,a, stats] = kruskalwallis(all_dir, all_days);
[results,~,~,gnames]=multcompare(stats, 'Dimension', [1,2])
tbl2 = array2table(results,"VariableNames", ...
["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl2.("Group A")=gnames(tbl2.("Group A"));
tbl2.("Group B")=gnames(tbl2.("Group B"))





     
%% Old method for bootstrapp:

%% Panel B: get bootstrapped OSI difference for retained spines and plot
all_temps = {retain_D1_D5, retain_D5_D10};
close all
sessions = {'binoc'};
days = {'D1', 'D5', 'D10'};
colors = {'g', 'r', 'b'}
all_delta_OSI = [];
all_delta_Ori = [];
for i = 1:length(sessions)
    %[dist_temp] = boot_strap_Ori_and_OSI_offset(all_temps, sessions{i}, days, 10000);

    figure("Name", sessions{i}, "Position", [440,656,560,140]);
    ori_delta = plot_Ori_offset(dist_temp,days,colors{i})
    all_delta_Ori.(sessions{i}) = ori_delta;

    figure("Name", sessions{i}, "Position", [440,656,560,200]);
    [osi_delta,all_sig] = plot_OSI_offset(dist_temp, days,colors{i});
    all_delta_OSI.(sessions{i}) = osi_delta;
    
end
%% Panel F: plot the bar chart mean OSI difference grouped by timepoint and session (maybe make into supplement)
close all
figure
count = 1;
feature = all_delta_Ori;
session = {'binoc'};
eye_pref = [];
day = [];
data = [];
for i = 1:length(sessions)
    temp = feature.(sessions{i});
    
    for ii = 1:length(temp)
        if ii == 1
            c = colors{i};
        else
            c = 'none';
        end
        plot_temp = temp{ii};
        plot_sig = all_sig{ii};
        bar(count,mean(plot_temp(plot_sig>=0)),'EdgeColor', colors{i}, 'FaceColor',c , 'EdgeColor', colors{i}); hold on;
        errorbar(count, mean(plot_temp(plot_sig>=0)),std(plot_temp(plot_sig>=0))/sqrt(length(plot_temp(plot_sig>=0))), '.k');
        %scatter(count*ones(length(temp.(unique_session{ii})),1), temp.(unique_session{ii}), 'jitter', 'on', 'jitterAmount',0.05)
        %swarmchart(count*ones(length(temp.(unique_session{ii})),1), temp.(unique_session{ii})); hold on;
        %plot([count-0.5, count+0.5],[mean(temp.(unique_session{ii})), mean(temp.(unique_session{ii}))], 'k')
        count = count + 1;
        eye_pref = [eye_pref;repmat(sessions(i),size(plot_temp(plot_sig>=0)))];
        day = [day; repmat(ii,size(plot_temp(plot_sig>=0)))];
        data = [data; {plot_temp(plot_sig>=0)}];
        [p,a]=signrank(plot_temp(plot_sig>=0));
        disp(['Sign rank test p-value for session ', sessions{i}, 'timepoint ', num2str(ii), ' : ', num2str(p)]);
    end
    count = count+1;
end
xticks([1,2])
xticklabels({'D5-D1', 'D10-D5'})
ylabel('\Delta OSI')
set(gca, 'FontSize',10)
ylim([0,40])

%%
function offset = get_dir_offset(spine, soma)
    offset = abs(spine-soma);
    if offset > 180
        offset = 360-offset;
    end
end

%%
function offset = get_ori_offset(spine, soma)
    spine = mod(spine,180);
    soma = mod(soma,180);
    offset = abs(spine-soma);
    if offset > 90
        offset = 180-offset;
    end
end
%%
function [ori_offset, dir_offset] = get_offset_maxresp(spine, soma,oris)
    [~,spine_ind] = find(spine == max(spine));
    [~,soma_ind] = find(soma == max(soma));
    spine_ori = oris(spine_ind);
    soma_ori = oris(soma_ind);
    dir_offset = get_dir_offset(spine_ori, soma_ori);
    ori_offset = get_ori_offset(spine_ori, soma_ori);
    
end




