%% Figure 1 - Characterizing eye-specific and binocular properties of dendritic spines
% Panel B: OI done for BM015 on 081022 (separate animal than for C)
% Panel C: Trial averaged and mean amplitude responses
% Panel D-F: Plot the praction of eye-specific responsive spines per day
% for each soma
%% Panel C: Trial averaged and mean amplitude responses
% BM018 p23 Cell 3 Dend 2 
close all
clear all
% get path details
path = '/Volumes/GoogleDrive-108846495442099470486/'; %update this path for where data will be stored
fov_path = 'My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Processed Data/'; %update this path for where data will be stored
mouse = 'BM018';
day = 'p23_081522';
cell_name = 'Cell3';
dend_name = 'Dend2_0';
soma_name = 'Soma_0';
dend_name_short = split(dend_name, '_');
dend_name_short = dend_name_short{1};
soma_name_short = split(soma_name, '_');
soma_name_short = soma_name_short{1};

% load dendritic fov
dend_binoc = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_binoc'],'normalized_data_by_stim.mat'));
dend_binoc_ori = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_binoc'],'mean_amp_ori_analysis_vonMises.mat'));
dend_contra = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_contra'],'normalized_data_by_stim.mat'));
dend_contra_ori = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_contra'],'mean_amp_ori_analysis_vonMises.mat'));
dend_ipsi = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_ipsi'],'normalized_data_by_stim.mat'));
dend_ipsi_ori = load(fullfile(path,fov_path,mouse,day,cell_name,dend_name,[dend_name_short, '_ipsi'],'mean_amp_ori_analysis_vonMises.mat'));
all_dend_sessions = {dend_binoc, dend_contra, dend_ipsi};
all_dend_vonMises = {dend_binoc_ori, dend_contra_ori, dend_ipsi_ori};
dend_rois = [24,14,16];

%load somatic fov
soma_binoc = load(fullfile(path,fov_path,mouse,day,cell_name,soma_name,[soma_name_short, '_binoc'],'normalized_data_by_stim.mat'));
soma_binoc_ori = load(fullfile(path,fov_path,mouse,day,cell_name,soma_name,[soma_name_short, '_binoc'],'mean_amp_ori_analysis_vonMises.mat'));
soma_contra = load(fullfile(path,fov_path,mouse,day,cell_name,soma_name,[soma_name_short, '_contra'],'normalized_data_by_stim.mat'));
soma_contra_ori = load(fullfile(path,fov_path,mouse,day,cell_name,soma_name,[soma_name_short, '_contra'],'mean_amp_ori_analysis_vonMises.mat'));
soma_ipsi = load(fullfile(path,fov_path,mouse,day,cell_name,soma_name,[soma_name_short, '_ipsi'],'normalized_data_by_stim.mat'));
soma_ipsi_ori = load(fullfile(path,fov_path,mouse,day,cell_name,soma_name,[soma_name_short, '_ipsi'],'mean_amp_ori_analysis_vonMises.mat'));
all_soma_sessions = {soma_binoc, soma_contra, soma_ipsi};
all_soma_vonMises = {soma_binoc_ori, soma_contra_ori, soma_ipsi_ori};

soma_rois = [1];

%set orientation details
oris_temp = 180-[0:45:315]; %assuming 8 directions, and PsychToolBox
oris_temp(oris_temp<0) = oris_temp(oris_temp <0)+360;
[oris_deg_sort,inds] = sort(oris_temp);

% plot dendritic roi and soma trial-averaged traces per session
c = {'g', 'r', 'b'};
g = figure("Name", 'Trial-averaged traces',"Position",[440,486,685,311]);
h = figure("Name", 'Mean responses per direction',"Position",[440,486,137,311]);
func = @(x,y,z) max(int8(x' > 0.5) + int8(y' < 0.05) + int8(z > 5));
for i = 1:length(all_soma_sessions)
    
    % plot soma's traces
    soma = all_soma_sessions{i};
    soma_ori = all_soma_vonMises{i};
    included_trials = soma.included_trials;
    indices = soma.indices;
    mean_amplitude = soma.mean_amplitude;
    std_amplitude = soma.std_amplitude;
    isi_dur = soma.isi_dur;
    stim_dur = soma.stim_dur;
    t_test = soma.t_test;

    num_oris = size(mean_amplitude,2);
    Datacell = soma.Datacell;
    [~,~,num_dims] = size(mean_amplitude);
    if num_dims>1
        ind = find_preferred_SF(mean_amplitude, t_test);
        mean_amplitude = squeeze(mean_amplitude(:,:,ind));
        std_amplitude = squeeze(std_amplitude(:,:,ind));
        t_test = squeeze(t_test(:,:,ind));
    else
        ind = 1;
    end
    std_amplitude = squeeze(std_amplitude(:,:,1)); %02/06/24 work around

    figure(g)
    subplot(length(dend_rois)+1, length(all_soma_sessions),i);
    plot_traces(1, indices, num_oris, isi_dur,stim_dur,ind, included_trials, Datacell,c{i});
   ylim([-1,6])
    figure(h)
    subplot(length(dend_rois)+1, 1,1);
    
    num_trials = sum(squeeze(included_trials(1,:,ind,:)),2);
    e = errorbar(oris_deg_sort,mean_amplitude(1,inds), std_amplitude(1,inds)./sqrt(num_trials'),'Color',c{i}); hold on;
    criteria = func(mean_amplitude(1,:),t_test(1,:),num_trials);
    if criteria==3
        scatter(0+i*10,4.5,'*',c{i});
    end
    e.CapSize = 0;
    e.MarkerSize = 10;
    e.Marker = '.';
    ylim([-1,4.5]);
    % plot spine's traces
    dend = all_dend_sessions{i};
    dend_ori = all_dend_vonMises{i};
    indices = dend.indices;
    included_trials = dend.included_trials;
    mean_amplitude = dend.mean_amplitude;
    std_amplitude = dend.std_amplitude;
    curve = dend_ori.Curve_Dir;
    gof = dend_ori.GOF_vonMises_Dir;
    t_test = dend.t_test;
    isi_dur = dend.isi_dur;
    stim_dur = dend.stim_dur;
    num_oris = size(mean_amplitude,2);
    Datacell = dend.Datacell;
    ind = 1;
    for ii = 1:length(dend_rois)
        figure(g)
        subplot(length(dend_rois)+1, length(all_soma_sessions),i+length(all_soma_sessions)*ii);
        plot_traces(dend_rois(ii), indices, num_oris, isi_dur,stim_dur,ind, included_trials, Datacell,c{i});
        ylim([-1,4])
        figure(h)
        subplot(length(dend_rois)+1, 1,ii+1);
        num_trials = sum(squeeze(included_trials(1,:,ind,:)),2);
        e = errorbar(oris_deg_sort,mean_amplitude(dend_rois(ii),inds), std_amplitude(dend_rois(ii),inds)./sqrt(num_trials'),'Color',c{i}); hold on;
        %plot(curve(dend_rois(ii),:), c{i});
       
        criteria = func(mean_amplitude(dend_rois(ii),:),t_test(dend_rois(ii),:),num_trials);
        e.CapSize = 0;
        e.Marker = '.';
        e.MarkerSize = 10;
        ylim([-1,2.5]);
        if criteria==3
        scatter(0+i*10,2.5,'*',c{i});
        end
    end
end
%% Panel D-F: Load Spine's and Soma's eye-specific & binoc preference
clear all 
close all
% load spine table across neurons and days and sessions
home_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/';
path_dir = 'Sur Lab/Development project/Binocular_Matching/Spine_imaging/Analyzed Data/';
spine_file = 'spine_mean_tuning_table_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_lax_criteria_zscored_trace_active_trials.mat';
soma_file = 'soma_mean_tuning_table_BM014_15_16_17_19_18_20_21_23_24_25_26_27_29_30_lax_criteria_zscored_trace_active_trials.mat';
filepath = fullfile(home_path, path_dir,spine_file);
load(filepath);
filepath = fullfile(home_path, path_dir,soma_file);
load(filepath)
% merge table into contra and ipsi session per spine
all_stim_table.roi_fovs_mouse_cell_days = strcat(num2str(all_stim_table.all_roi_inds), '_', ...
                                                all_stim_table.mouse_cell, '_',...
                                                all_stim_table.all_fovs, '_',...
                                                all_stim_table.days);


all_soma_stim_table.mouse_cell_days = strcat(all_soma_stim_table.mouse_cell, '_',all_soma_stim_table.days);
%% Panel D-F: Plot responsive spine by session by mouse by day
temp = all_stim_table;
features = {'spine_mean_resp','resp','soma_resp','roi_fovs_mouse_cell_days', 'days', 'mouse_cell'};
join_key = 'roi_fovs_mouse_cell_days';
contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, join_key);

%create one hot encoding for spines and soma eye pref
contra_ipsi_binoc.ci_spine = strcat(num2str(contra_ipsi_binoc.resp_contra), num2str(contra_ipsi_binoc.resp_ipsi));
contra_ipsi_binoc.cib_spine = strcat(num2str(contra_ipsi_binoc.resp_contra), num2str(contra_ipsi_binoc.resp_ipsi),num2str(contra_ipsi_binoc.resp));

contra_ipsi_binoc.ci_soma = strcat(num2str(contra_ipsi_binoc.soma_resp_contra), num2str(contra_ipsi_binoc.soma_resp_ipsi));
contra_ipsi_binoc.cib_soma = strcat(num2str(contra_ipsi_binoc.soma_resp_contra), num2str(contra_ipsi_binoc.soma_resp_ipsi),num2str(contra_ipsi_binoc.soma_resp));

%plot fraction of binoc, contra, ipsi responsive spines for each cell by
%day
spine_resp_by_session_by_neuron_by_day



%% plot trial-averaged traces
function plot_traces(roi, indices, num_oris, isi_dur,stim_dur,jj, included_trials, Datacell,c)
    % indices: Cell array stating how many trials there were unique stim
    % included trials: how many trials were included in unique stim
    % DataCell: Cell array with the stim traces aligned to the stim onset
    % num oris: number of unique directions 
    % isi_dur: duration of ISI
    % stim_dur: duration of stimulus
    % jj: default is 1 (but changes on number of contrasts, spatial freq
    % c: color 
    time_axis = -round(isi_dur):.1:round(stim_dur+isi_dur);
    index_isi = find(time_axis < 0 & time_axis >= -isi_dur/3); 
    index_stim = find(time_axis >=0 & time_axis <= 1.5); % plot from stim onset to 1.5 sec after 
    end_stim = find(time_axis==1);
    time_period = time_axis + isi_dur; %used for plotting the traces
    n=0; 
    for ii = 1:num_oris
        n = n+1;        
        %create empty matrix to contain interpolated stim traces, 
        % 10 frames before stim onset, 10 frames during stim onset, 10
        % frames after stim 
        
        stim_traces=zeros(round(20*isi_dur + 10*stim_dur),length(indices{n}));
        for kk = 1:length(indices{n})
            %interpolate each trial on a 0.1 S scale  
            stim_traces(:,kk)=interp1(Datacell{ii,jj,kk}(:,1),Datacell{ii,jj,kk}(:,roi+1),(round(-isi_dur*10):1:floor(stim_dur*10+isi_dur*10))/10);
            
            %plot stim trace
            %plot(time_period([index_isi,index_stim]),stim_traces([index_isi,index_stim],kk),'Color',[0.5,0.5,0.5],'LineWidth',0.25); hold on;
            
        end
        include_trials = find(squeeze(included_trials(roi,ii,jj,:)));
        
        %Take mean of traces 
        mean_trace=nanmean(stim_traces(:,include_trials),2);
        SEM_trace = nanstd(stim_traces(:,include_trials),[],2)./sqrt(length(include_trials));
        
        %plot mean trace
        plot_SEM(time_period([index_isi,index_stim]),mean_trace([index_isi,index_stim]),SEM_trace([index_isi,index_stim]),c,c);

        %plot mean trace
        plot(time_period([index_isi,index_stim]),mean_trace([index_isi,index_stim]),'Color',c,'LineWidth',1);
        patch([time_period(index_stim(1)),time_period(end_stim),time_period(end_stim),time_period(index_stim(1))], [-5, -5,10, 10],[0.5,0.5,0.5], 'FaceAlpha', 0.3,'LineStyle', 'none')          
        time_period = time_period + round(stim_dur*10+(isi_dur/2)*10)/10;
        
    end
end

