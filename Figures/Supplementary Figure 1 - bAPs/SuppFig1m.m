%% Supplementary Figure 1
% Panel B: compare the calcium traces of the dendrite and spine
% Panel C: Plot the dendrite, spine, and spine substracted trace
% Plot D: Plot the OSI, DSI, Ori Pref, Dir Pref for spines with or without
% bAPs
%% Panel B and C: compare the calcium traces of the dendrite and spine
clear all 
close all
path = '/Volumes/GoogleDrive-108846495442099470486/';
fov_path = 'My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Processed Data/';
mouse = 'BM030';
day = 'p24_030623';
cell_name = 'Cell2';
dend_name = 'Dend2_0';
session = 'contra';

dend_name_short = split(dend_name, '_');
dend_name_short = dend_name_short{1};
load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_', session],'normalized_data_by_stim.mat'));
roi = 17;

%% find F for each ROI 
%calculate number of ROIs
[TimeLength,numROIs]=size(ROIdata);
%preallocate matrix for blank frame data (period b/t stims, ISI)
F_data=zeros(1,(numROIs));
%loop to collect ISI 
for ii=1:length(blankframes)
    Chunk = (ROIdata(timestamps(:,1) >=(blankframes(ii)+isi_dur/2)&timestamps(:,1) <=(blankframes(ii)+isi_dur),:));
    F_data=[F_data;Chunk];
end
%exclude first empty row
F_data=F_data(2:end,:);
%calculate mean of each ROI's ISI
meanF=nanmean(F_data(1:end,:),1);
clear F_data chunk ii

%% Delta F/F
ROIdata_F=zeros(TimeLength,numROIs);
for ii=1:numROIs
ROIdata_F(:,ii)=(ROIdata(:,ii)-meanF(1,ii))./meanF(1,ii);
end
clear ii meanF

%% Subtract dendrite contamination
%identify dendritic shaft ROIs, and for each spine subtract its closest
%shaft, use robust regression to subtract shaft signal, and calculate
%remaining shaft correlation
roi_temp = (roi*2-1);

b = robustfit(ROIdata_F(:,shaft(roi_temp,2)), ROIdata_F(:,roi_temp));
ROIdata_fit = ROIdata_F(:,roi_temp) - b(2).*ROIdata_F(:,shaft(roi_temp,2));
x_space = linspace(-1,max(ROIdata_F(:,shaft(roi_temp,2))));
y_space = b(1) + b(2)*x_space;
figure("Position", [420,258,206,539])
subplot(4,1,1)
scatter(ROIdata_F(:,shaft(roi_temp,2)),ROIdata_F(:,roi_temp)); hold on;
plot(x_space,y_space, '--k', 'LineWidth',1);
ylim([-1,4])
xlim([-1,4])
subplot(4,1,2)
plot(ROIdata_F(:,shaft(roi_temp,2)), 'b','LineWidth',2), hold on;
xlim([0,1000])
ylim([-1,4])
subplot(4,1,3) 
plot(ROIdata_F(:,roi_temp),'LineWidth',2)
xlim([000,1000])
ylim([-1,4])
subplot(4,1,4) 
plot(ROIdata_fit,'LineWidth',2)
xlim([0,1000])
ylim([-1,4])

%% Panel D-F: bAP Comparisons
path_dir = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
file = 'spine_properties_table.mat';
filepath = fullfile(path_dir,file);
load(filepath);

all_stim_table = all_stim_table(~isnan(all_stim_table.all_OSI_vector),:);
all_spines_resp_cells = all_stim_table(all_stim_table.soma_resp==1,:);
%% Panel D: determine number of bAPs per cell per day
close all
unique_mice = unique(all_stim_table.mouse_cell);
days = {'D1', 'D5', 'D10'};
session = {"binoc";"contra"; "ipsi"};
colors = {'g', 'r', 'b'};
bAPs_by_session = [];
for s = 1:length(session)
    num_bAPs_per_cell_per_day = NaN(length(unique_mice), length(days));
    soma_resp_by_mouse_by_day = NaN(length(unique_mice), length(days));
    figure("Position", [794,312,122,331])
    subplot(2,1,1)
    for i = 1:length(unique_mice)
        for ii = 1:length(days)
            temp = all_stim_table(strcmp(all_stim_table.mouse_cell, unique_mice{i})...
                &strcmp(all_stim_table.days, days{ii})&strcmp([all_stim_table.session{:}]', session{s}),:);
            bAPs = temp.all_bAP_trial;
            included_trials = temp.all_included_trial;
            bAPs = find_number_bAPs_pref_dir(temp);
            num_bAPs_per_cell_per_day(i,ii) = mean(unique(bAPs));
            if ~isempty(temp)
                soma_resp_by_mouse_by_day(i,ii) = temp.soma_resp(1);
            end
        end
    end
    plot_resp_by_mouse(num_bAPs_per_cell_per_day,soma_resp_by_mouse_by_day,colors{s});
    xticks([1:3])
    xticklabels(days)
    xlim([0,4])
    ylabel('% bAP trials')
    ylim([0,100])
    bAPs_by_session.(session{s}) =  num_bAPs_per_cell_per_day;
end
%% Panel D: Compare fraction of active trials to soma's pref direction w/ and w/o bAPs
close all
bAPs_by_session = [];
temp = all_stim_table(all_stim_table.soma_resp==1,:);
num_bAPs = find_number_bAPs_pref_dir(temp);

[active_trials_with_bAPs, active_trials_bAPs_removed] = find_active_trials_pref_dir(temp);
figure
subplot(2,2,1)
scatter(active_trials_with_bAPs, active_trials_bAPs_removed,[],num_bAPs, 'filled', 'MarkerFaceAlpha',0.3); hold on;
h = colorbar;
h.Label.String = "#bAPs (soma's pref dir)";
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";

[p,S] = polyfit(active_trials_with_bAPs,active_trials_bAPs_removed,1); 
f = polyval(p,active_trials_with_bAPs); 
plot(active_trials_with_bAPs,f,'k', 'LineWidth',2)
rsq = 1 - (S.normr/norm(active_trials_bAPs_removed - mean(active_trials_bAPs_removed)))^2;
text(20,100,{['r^2: ', num2str(rsq)]});
ylabel('without bAPs')
xlabel('with bAPs')
title('% active trials (soma''s pref dir)')
%%

tbl = table(active_trials_with_bAPs,active_trials_bAPs_removed, num_bAPs);
mdl = fitlm(tbl,'active_trials_bAPs_removed ~ active_trials_with_bAPs + num_bAPs')

plot(mdl, '.k', LineWidth=2)



%% Panel E: Compare orientation tuning properties
figure
num_bAPs_resp_cells = find_number_bAPs_pref_dir(all_spines_resp_cells);
resp_spines_bAPs_resp_cell = num_bAPs_resp_cells(all_spines_resp_cells.resp==1);

%responsive somas
resp_spine_table = all_spines_resp_cells(all_spines_resp_cells.resp==1,:);

%transform ori preference for bAP remove to orientation (did not do this in
%create_table_of_all3sessions)
Oris =resp_spine_table.all_Ori_pref_vector_bAPs_removed;
Oris = Oris-90;
Oris(Oris<0) = Oris(Oris<0)+180; 
resp_spine_table.all_Ori_pref_vector_bAPs_removed = Oris;

subplot(2,2,1)
scatter(resp_spine_table.all_Ori_pref_vector, resp_spine_table.all_Ori_pref_vector_bAPs_removed, [],resp_spines_bAPs_resp_cell ); hold on;
mdl = fitlm(resp_spine_table.all_Ori_pref_vector,resp_spine_table.all_Ori_pref_vector_bAPs_removed)
rsq = mdl.Rsquared.Ordinary
plot([0,180], [0,180], 'k', LineWidth=2)
text(1,200,{rsq});
title('Orientation preference')
ylabel('without bAPs');
xlabel('with bAPs');
h = colorbar;
h.Label.String = "#bAPs for preferred direction";
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";

subplot(2,2,2)
scatter(resp_spine_table.all_Dir_pref_vector, resp_spine_table.all_Dir_pref_vector_bAPs_removed, [],resp_spines_bAPs_resp_cell ); hold on;
mdl = fitlm(resp_spine_table.all_Dir_pref_vector,resp_spine_table.all_Dir_pref_vector_bAPs_removed)
rsq = mdl.Rsquared.Ordinary
plot([0,360], [0,360], 'k', LineWidth=2)
text(1,360,{rsq});
title('Direction preference')
ylabel('without bAPs');
xlabel('with bAPs');
h = colorbar;
h.Label.String = "#bAPs for preferred direction";
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";

subplot(2,2,3)
scatter(resp_spine_table.all_OSI_vector, resp_spine_table.all_OSI_vector_bAPs_removed, [],resp_spines_bAPs_resp_cell ); hold on;
mdl = fitlm(resp_spine_table.all_OSI_vector,resp_spine_table.all_OSI_vector_bAPs_removed)
rsq = mdl.Rsquared.Ordinary;
plot([0,1], [0,1], 'k', LineWidth=2)
text(0,1,{rsq});
title('OSI')
ylabel('without bAPs');
xlabel('with bAPs');
h = colorbar;
h.Label.String = "#bAPs for preferred direction";
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";

subplot(2,2,4)
scatter(resp_spine_table.all_DSI_vector, resp_spine_table.all_DSI_vector_bAPs_removed, [],resp_spines_bAPs_resp_cell ); hold on;
plot([0,1], [0,1], 'k', LineWidth=2)
mdl = fitlm(resp_spine_table.all_DSI_vector,resp_spine_table.all_DSI_vector_bAPs_removed)
rsq = mdl.Rsquared.Ordinary
text(0,1,{rsq});
title('DSI')
ylabel('without bAPs');
xlabel('with bAPs');
h = colorbar;
h.Label.String = "#bAPs for preferred direction";
h.Label.Rotation = 270;
h.Label.VerticalAlignment = "bottom";

%%
function plot_resp_by_mouse(fraction_resp_by_mouse_by_day,soma_resp_by_mouse_by_day,c)
mean_fract = squeeze(nanmean(fraction_resp_by_mouse_by_day,1));
num_cells_by_day = sum(~isnan(fraction_resp_by_mouse_by_day));
sem_fract = squeeze(nanstd(fraction_resp_by_mouse_by_day,1))./sqrt(num_cells_by_day);
for i = 1:size(fraction_resp_by_mouse_by_day,1)
    plot([1:3], fraction_resp_by_mouse_by_day(i,:),'Color', [0.5,0.5,0.5],'LineWidth',0.5); hold on;

    for ii = 1:size(fraction_resp_by_mouse_by_day,2)
         if soma_resp_by_mouse_by_day(i,ii)==1
            scatter(ii, fraction_resp_by_mouse_by_day(i,ii), 100,'.',c, 'MarkerFaceAlpha',0.6); hold on;
         else
            scatter(ii, fraction_resp_by_mouse_by_day(i,ii), 100,'.','k','MarkerFaceAlpha',0.6); hold on;
         end
    end
    %errorbar([1:size(mean_fract,2)], mean_fract, sem_fract, '-o', 'MarkerFaceColor', 'auto', 'MarkerSize',5, 'LineWidth',2); hold on;
end
disp(['mean ', num2str(mean_fract)]);
disp(['sem ', num2str(sem_fract)]);

end

%%
function num_bAPs = find_number_bAPs_pref_dir(spines)
mean_resp = spines.soma_mean_resp;
bAP_trials = spines.all_bAP_trial;
%included_trials =  spines.all_included_trial;
pref_dir = cellfun(@(x) find(x == max(x)), mean_resp);
num_bAPs = cellfun(@(x,y) sum(x(:,y,:))*100/numel(x(:,y,:)), bAP_trials,num2cell(pref_dir));
end
%%
function [active_trials_with_bAPs, active_trials_bAPs_removed] = find_active_trials_pref_dir(spines)
mean_resp = spines.spine_mean_resp;
bAP_trials = spines.all_bAP_trial;
included_trials =  spines.all_included_trial;
bAP_included_trials = cellfun(@(x,y) squeeze(x)'&~y, included_trials, bAP_trials, 'UniformOutput',false);
active_trials = spines.all_active_spine_trials;
pref_dir = cellfun(@(x) find(x == max(x)), mean_resp);
active_trials_with_bAPs = cellfun(@(x,z,y) sum(x(squeeze(z(:,y,:))==1))*100/numel(x(squeeze(z(:,y,:))==1)), active_trials, included_trials,num2cell(pref_dir));
active_trials_bAPs_removed = cellfun(@(x,z,y,a) sum(x(squeeze(z(:,y,:))==1))*100/numel(x(squeeze(a(:,y,:))==1)), active_trials,bAP_included_trials,num2cell(pref_dir),included_trials);
end