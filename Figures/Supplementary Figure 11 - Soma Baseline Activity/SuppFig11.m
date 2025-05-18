%% Longitudinal spines tuning by day
clear all
close all
path = '/Users/ktsimring/Documents/GitHub/Tsimring-Enhanced-synaptic-dynamics-drive-the-reorganization-of-binocular-circuit/Mat Files';
distance_path = fullfile(path, 'Chronic Imaging/FOV_alignment/'); %update path
load(fullfile(path, "soma_properties_table.mat"), "all_soma_stim_table");


%%
path = '/Volumes/GoogleDrive-108846495442099470486/';
fov_path = 'My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Processed Data/';
mouse = [{'BM015'};{'BM029'}];
days = [{'p24_072922'}, {'p29_080322'},{'p34_080822'};{'p22_030523'}, {'p28_031123'},{'p33_031623'}];
cell_name = [{'Cell8'};{'Cell3'}];
soma_name = 'Soma_0';
soma_name_short = split(soma_name, '_');
soma_name_short = soma_name_short{1};

for i = 1:length(mouse)
    day = days(i,:);
    figure
    for ii = 1:length(day)
        mouse_x = mouse{i};
        day_x = day{ii};
        cell_x = cell_name{i};
        soma_contra = load(fullfile(path,fov_path,mouse_x,day_x,cell_x,soma_name,[soma_name_short, '_contra'],'normalized_data_by_stim.mat'));
        subplot(length(day),1,ii)
        plot(soma_contra.timestamps,soma_contra.ROIdata_Z(:,1),'Color','k', 'LineWidth',1); hold on;
        

         for iii = 1:length(soma_contra.startevent)
            patch([soma_contra.startevent(iii),soma_contra.blankframes(iii),soma_contra.blankframes(iii),soma_contra.startevent(iii)], [-5, -5,15, 15],[0.5,0.5,0.5], 'FaceAlpha', 0.3,'LineStyle', 'none') 
         end
         %xlim([soma_contra.startevent(5),soma_contra.startevent(40)])
         ylim([-2,15])
         xlim([soma_contra.startevent(1), soma_contra.startevent(end)])
    end
end
%%
close all
path = '/Volumes/GoogleDrive-108846495442099470486/';
fov_path = 'My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Processed Data/';
mouse = [{'BM015'};{'BM029'}];
days = [{'p24_072922'}, {'p29_080322'},{'p34_080822'};{'p22_030523'}, {'p28_031123'},{'p33_031623'}];
cell_name = [{'Cell8'};{'Cell3'}];
soma_name = 'Soma_0';
soma_name_short = split(soma_name, '_');
soma_name_short = soma_name_short{1};
oris_deg_psycopy = (180-[0:45:315]); %% based on psycho conversion
oris_deg_psycopy(oris_deg_psycopy<0) = oris_deg_psycopy(oris_deg_psycopy<0)+360;
[~,inds] = sort(oris_deg_psycopy);

for i = 1:length(mouse)
    day = days(i,:);
    figure
     set(gcf,'Position', [667,377,286,420]);
    for ii = 1:length(day)
        mouse_x = mouse{i};
        day_x = day{ii};
        cell_x = cell_name{i};
        soma_contra = load(fullfile(path,fov_path,mouse_x,day_x,cell_x,soma_name,[soma_name_short, '_contra'],'mean_amp_ori_analysis_vonMises.mat'));
        mean_data = load(fullfile(path,fov_path,mouse_x,day_x,cell_x,soma_name,[soma_name_short, '_contra'],'normalized_data_by_stim.mat'), 'mean_amplitude', 'std_amplitude');
         kappa = soma_contra.coeffs_dir(1,4);
        theta_HWHM = acos(1 - log(2) / kappa);
        theta_HWHM_deg = rad2deg(theta_HWHM);
        subplot(length(day),1,ii)
        von_vector = soma_contra.Curve_Dir(1,:);
       % mean_data = mean_data.mean_amplitude(1,inds);
        plot(von_vector);hold on; scatter(oris_deg_psycopy(inds),mean_data.mean_amplitude(1,inds), 'filled')
       errorbar(oris_deg_psycopy(inds),mean_data.mean_amplitude(1,inds), mean_data.std_amplitude(1,inds), "LineStyle",'none','CapSize',0)
        title('HWHM: ', num2str(theta_HWHM_deg));
        %plot(soma_contra.timestamps,soma_contra.ROIdata_Z(:,1),'Color','k', 'LineWidth',1); hold on;
        if i == 1
            ylim([-1,3])
        else
            ylim([-1,10])
        end
        xlim([0,360])
        xticks([0:45:315])
        xticklabels([0:45:315])
    end
   
end

%%
all_soma_stim_table.mouse_cell_days = strcat(all_soma_stim_table.mouse_cell,'_', all_soma_stim_table.days);
repeated_image = zeros(height(all_soma_stim_table),1);
unique_mouse_cell = unique(all_soma_stim_table.mouse_cell);
for i = 1:length(unique_mouse_cell)
    inds_D1 = find(ismember(all_soma_stim_table.mouse_cell_days, [unique_mouse_cell{i},'_D1']));
    inds_D5 = find(ismember(all_soma_stim_table.mouse_cell_days, [unique_mouse_cell{i},'_D5']));
    inds_D10 = find(ismember(all_soma_stim_table.mouse_cell_days, [unique_mouse_cell{i},'_D10']));
    if ~isempty(inds_D1)&~isempty(inds_D5)&~isempty(inds_D10)
       repeated_image(inds_D1) = 1;
       repeated_image(inds_D5) = 1;
       repeated_image(inds_D10) = 1;
    end
end

all_stim_table_filter = all_soma_stim_table(repeated_image==1,:);
%%
feature = 'all_ISI_data';
features = {'resp','all_z_scored_trace', 'all_trial_data','all_ISI_data','all_OSI_vector','all_coeffs_dir','mouse_cell_days', 'days', 'mouse_cell'};
join_key = 'mouse_cell_days';
alpha = 0.4; % for smoothing
std_thresh = 3; %find peaks threshold above 3 std
conseq_frames = 3; % number of frames 
fs = 7.68; %I should change this to actual frame rate

contra_ipsi_binoc = concate_contra_ipsi_binoc_v2(all_stim_table_filter, features, join_key);
    
concatenate_z_scored_trace = cellfun(@(x,y,z) [x;y;z], contra_ipsi_binoc.([feature,'_contra']),...
contra_ipsi_binoc.([feature,'_ipsi']),contra_ipsi_binoc.(feature), 'UniformOutput',false);
   
%smooth trace using exponential weighted filter 
concatenated_smoothed_z_scored_trace = cellfun(@(x) exp_weight_filter(x, alpha), concatenate_z_scored_trace,'UniformOutput',false);

%count spine calcium event rate
[pks,locs] = cellfun(@(x) findpeaks(x, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames), concatenated_smoothed_z_scored_trace,'UniformOutput',false);
ca_rate = cellfun(@(x,y) length(x)*fs/length(y), pks, concatenate_z_scored_trace);
contra_ipsi_binoc.ca_rate = ca_rate;
contra_ipsi_binoc.cib_resp = arrayfun(@(x,y,z) sum([x,y,z]),contra_ipsi_binoc.resp_contra,contra_ipsi_binoc.resp_ipsi,contra_ipsi_binoc.resp);

%% Panel B: Plot somas that lose response 

mouse_cell = unique(contra_ipsi_binoc.mouse_cell);
days = {'D1', 'D5', 'D10'};
all_ca_rate_by_day = [];
soma_resp = [];
figure('Position',[440,384,279,208])
hold on;
for i = 1:length(mouse_cell)
    ca_rate_by_day = [];
    resp_by_day = [];
    if i == 14
        k = 12;
    end
    for ii = 1:length(days)
        temp = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.mouse_cell_days, [mouse_cell{i}, '_', days{ii}]),:);
        ca_rate = temp.ca_rate;
        ca_rate_by_day = [ca_rate_by_day, ca_rate];
        resp_by_day = [resp_by_day,temp.cib_resp];
    end
    if resp_by_day(1)>0 && resp_by_day(3)==0
        plot([1:3],ca_rate_by_day, 'r');
        scatter([1:3], ca_rate_by_day, 'r', 'filled')
        resp = 0;
    else
        plot([1:3],ca_rate_by_day, 'k');
        scatter([1:3], ca_rate_by_day, 'k', 'filled')
        resp = 1;
    end
    all_ca_rate_by_day = [all_ca_rate_by_day; ca_rate_by_day];
    soma_resp = [soma_resp;resp_by_day>0]

end
xticks([1:3])
xticklabels(days)
xlim([0,4])
ylabel('Calcium Rate During ISI (Hz)')

%% Stats (using chatgpt)
calciumRates =  all_ca_rate_by_day;
responseType = soma_resp(:,3);

% Create neuron IDs
numNeurons = size(calciumRates, 1);
numTimePoints = size(calciumRates, 2);
neuronIDs = repelem((1:numNeurons)', numTimePoints);

% Create a time variable (1,2,3 for days)
time = repmat((1:numTimePoints)', numNeurons, 1);

% Flatten calcium rates into a vector
calciumRatesVec = calciumRates(:);

% Expand response type for each neuron across time points
responseTypeExpanded = repelem(responseType, numTimePoints);

% Convert response type to categorical
responseTypeExpanded = categorical(responseTypeExpanded);
%%
tbl = table(neuronIDs, time, calciumRatesVec, responseTypeExpanded, ...
    'VariableNames', {'Neuron', 'Time', 'CalciumRate', 'ResponseType'});
%%
% Display results
disp(ranovatbl);
% Fit a linear mixed model with random intercepts for neurons
lme = fitlme(tbl, 'CalciumRate ~ ResponseType * Time + (1|Neuron)');

% Display results
disp(lme);
%%
% Create a table in wide format
tbl = array2table(calciumRates, 'VariableNames', {'Day1', 'Day2', 'Day3'});
tbl.ResponseType = categorical(responseType); % Add Response Type

%% Panel B: Plot half-width-half-max for somas that lose response 

mouse_cell = unique(contra_ipsi_binoc.mouse_cell);
days = {'D1', 'D5', 'D10'};
all_theta_HWHM_deg_by_day = [];
soma_resp = [];
figure('Position',[440,384,279,208])
hold on;
for i = 1:length(mouse_cell)
     theta_HWHM_deg_by_day = [];
    resp_by_day = [];
    if i == 14
        k = 12;
    end
    for ii = 1:length(days)
        temp = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.mouse_cell_days, [mouse_cell{i}, '_', days{ii}]),:);
        von_mises_coeffs = temp.all_coeffs_dir;
        kappa = cellfun(@(x) x(4), von_mises_coeffs);
        theta_HWHM = acos(1 - log(2) / kappa);
        theta_HWHM_deg = rad2deg(theta_HWHM);
        theta_HWHM_deg_by_day = [ theta_HWHM_deg_by_day,  theta_HWHM_deg];
        resp_by_day = [resp_by_day,temp.resp];
    end
    if resp_by_day(1)>0 && resp_by_day(3)==0
        plot([1:3],  theta_HWHM_deg_by_day, 'r');
        scatter([1:3],   theta_HWHM_deg_by_day, 'r', 'filled')
        resp = 0;
    else
        plot([1:3],  theta_HWHM_deg_by_day, 'k');
        scatter([1:3],   theta_HWHM_deg_by_day, 'k', 'filled')
        resp = 1;
    end
    all_theta_HWHM_deg_by_day = [all_theta_HWHM_deg_by_day;  theta_HWHM_deg_by_day];
    soma_resp = [soma_resp;resp_by_day>0]

end
xticks([1:3])
xticklabels(days)
xlim([0,4])
ylabel('Binoc Half-Width Half-Max')
%% Stats (using chatgpt)
calciumRates =  all_theta_HWHM_deg_by_day;
responseType = soma_resp(:,3);

% Create neuron IDs
numNeurons = size(calciumRates, 1);
numTimePoints = size(calciumRates, 2);
neuronIDs = repelem((1:numNeurons)', numTimePoints);

% Create a time variable (1,2,3 for days)
time = repmat((1:numTimePoints)', numNeurons, 1);

% Flatten calcium rates into a vector
calciumRatesVec = calciumRates(:);

% Expand response type for each neuron across time points
responseTypeExpanded = repelem(responseType, numTimePoints);

% Convert response type to categorical
responseTypeExpanded = categorical(responseTypeExpanded);
%%
tbl = table(neuronIDs, time, calciumRatesVec, responseTypeExpanded, ...
    'VariableNames', {'Neuron', 'Time', 'CalciumRate', 'ResponseType'});
%%
% Create a table in wide format
tbl = array2table(calciumRates, 'VariableNames', {'Day1', 'Day2', 'Day3'});
tbl.ResponseType = categorical(responseType); % Add Response Type

% Fit the repeated measures model
rm = fitrm(tbl, 'Day1-Day3 ~ ResponseType', 'WithinDesign', [1 2 3]);

% Run repeated measures ANOVA
ranovatbl = ranova(rm);

% Display results
disp(ranovatbl);

%%
% Display results
disp(ranovatbl);
% Fit a linear mixed model with random intercepts for neurons
lme = fitlme(tbl, 'CalciumRate ~ ResponseType * Time + (1|Neuron)');

% Display results
disp(lme);

%% Panel C: 
close all
% load spine table across neurons and days and sessions
clear all
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
%% responsive spine by session by mouse by day
temp = all_stim_table;
features = {'spine_mean_resp','resp','soma_resp','roi_fovs_mouse_cell_days', 'days', 'mouse_cell'};
join_key = 'roi_fovs_mouse_cell_days';
contra_ipsi_binoc = concate_contra_ipsi_binoc_v2(temp, features, join_key);

%create one hot encoding for spines and soma eye pref
contra_ipsi_binoc.ci_spine = strcat(num2str(contra_ipsi_binoc.resp_contra), num2str(contra_ipsi_binoc.resp_ipsi));
contra_ipsi_binoc.cib_spine = strcat(num2str(contra_ipsi_binoc.resp_contra), num2str(contra_ipsi_binoc.resp_ipsi),num2str(contra_ipsi_binoc.resp));

contra_ipsi_binoc.ci_soma = strcat(num2str(contra_ipsi_binoc.soma_resp_contra), num2str(contra_ipsi_binoc.soma_resp_ipsi));
contra_ipsi_binoc.cib_soma = strcat(num2str(contra_ipsi_binoc.soma_resp_contra), num2str(contra_ipsi_binoc.soma_resp_ipsi),num2str(contra_ipsi_binoc.soma_resp));

%% Plot % responsive spines responsive to sessions by mouse by day (pooled across sessions)
% color code by whether soma was responsive or not
close all
unique_mouse = unique(contra_ipsi_binoc.mouse_cell);
session = { '','_contra', '_ipsi'}; %empty is the binoc 
unique_days = {'D1', 'D5', 'D10'};
colors = {'g','r', 'b' };
all_soma_resp_by_session = [];
all_session = [];
figure("Position", [283,132,150,344])
  all_fraction_resp = [];
  all_soma_resp = [];
for i = 1:length(session)
    all_resp = 0;
    resp_type = ['resp', session{i}];

    [fraction_resp_by_mouse_by_day,soma_resp_by_mouse_by_day] = get_eye_resp_mouse_day(unique_mouse,unique_days, contra_ipsi_binoc, all_resp,resp_type);
    all_soma_resp = [all_soma_resp;soma_resp_by_mouse_by_day];
    all_fraction_resp = [all_fraction_resp;fraction_resp_by_mouse_by_day];
end
    D1_D5_fraction_resp = [all_fraction_resp(:,[1,2]);all_fraction_resp(:,[2,3])];
    D1_D5_soma_resp = [all_soma_resp(:,[1,2]); all_soma_resp(:,[2,3])];
    
    D1_D5_soma_resp(any(isnan(D1_D5_soma_resp ), 2), :) = [];
    D1_D5_fraction_resp(any(isnan(D1_D5_fraction_resp), 2), :) = [];

    keep_resp_inds = find(D1_D5_soma_resp(:,1)==1&D1_D5_soma_resp(:,2)==1);
    lose_resp_inds = find(D1_D5_soma_resp(:,1)==1&D1_D5_soma_resp(:,2)==0);
    gain_resp_inds = find(D1_D5_soma_resp(:,1)==0&D1_D5_soma_resp(:,2)==1);

    delta_D1_D5_fraction_resp = (D1_D5_fraction_resp(:,2)-D1_D5_fraction_resp(:,1));
    soma_resp_fract_spines = delta_D1_D5_fraction_resp(keep_resp_inds);
    soma_unresp_fract_spines = delta_D1_D5_fraction_resp(lose_resp_inds);
    soma_gain_fract_spines = delta_D1_D5_fraction_resp(gain_resp_inds);
    
    mean_soma_resp_fract_spines = mean(soma_resp_fract_spines);
    mean_soma_unresp_fract_spines = mean(soma_unresp_fract_spines);
    sem_soma_resp_fract_spines = std(soma_resp_fract_spines)/sqrt(length(soma_resp_fract_spines));
    sem_soma_unresp_fract_spines = std(soma_resp_fract_spines)/sqrt(length(soma_unresp_fract_spines));
    
   
    bar([mean_soma_unresp_fract_spines,mean_soma_resp_fract_spines]); hold on;
    errorbar([mean_soma_unresp_fract_spines,mean_soma_resp_fract_spines],[sem_soma_unresp_fract_spines, sem_soma_resp_fract_spines], '.k');
    scatter(ones(size(soma_unresp_fract_spines)), soma_unresp_fract_spines, '.k','jitter', 'on', 'jitterAmount', 0.05);
    scatter(2*ones(size(soma_resp_fract_spines)), soma_resp_fract_spines,'.k', 'jitter', 'on', 'jitterAmount', 0.05);

    xticklabels({'soma lose response', 'soma maintain response'})
    ylabel('\Delta % responsive spines')
    [p,a] = ranksum(soma_resp_fract_spines, soma_unresp_fract_spines)

  



%% Plot % responsive spines responsive to sessions by mouse by day
% color code by whether soma was responsive or not
close all
unique_mouse = unique(contra_ipsi_binoc.mouse_cell);
session = { '','_contra', '_ipsi'}; %empty is the binoc 
unique_days = {'D1', 'D5', 'D10'};
colors = {'g','r', 'b' };
all_soma_resp_by_session = [];
all_session = [];
figure("Position", [283,132,150,344])
   
for i = 1:length(session)
    all_resp = 0;
    resp_type = ['resp', session{i}];

    [fraction_resp_by_mouse_by_day,soma_resp_by_mouse_by_day] = get_eye_resp_mouse_day(unique_mouse,unique_days, contra_ipsi_binoc, all_resp,resp_type);
    D1_D5_fraction_resp = [fraction_resp_by_mouse_by_day(:,[1,2]);fraction_resp_by_mouse_by_day(:,[2,3])];
    D1_D5_soma_resp = [soma_resp_by_mouse_by_day(:,[1,2]);soma_resp_by_mouse_by_day(:,[2,3])];
    
    D1_D5_soma_resp(any(isnan(D1_D5_soma_resp ), 2), :) = [];
    D1_D5_fraction_resp(any(isnan(D1_D5_fraction_resp), 2), :) = [];

    keep_resp_inds = find(D1_D5_soma_resp(:,1)==1&D1_D5_soma_resp(:,2)==1);
    lose_resp_inds = find(D1_D5_soma_resp(:,1)==1&D1_D5_soma_resp(:,2)==0);
    gain_resp_inds = find(D1_D5_soma_resp(:,0)==1&D1_D5_soma_resp(:,2)==1);

    delta_D1_D5_fraction_resp = (D1_D5_fraction_resp(:,2)-D1_D5_fraction_resp(:,1));
    soma_resp_fract_spines = delta_D1_D5_fraction_resp(keep_resp_inds);
    soma_unresp_fract_spines = delta_D1_D5_fraction_resp(lose_resp_inds);
    soma_gain_fract_spines = delta_D1_D5_fraction_resp(gain_resp_inds);
    
    mean_soma_resp_fract_spines = mean(soma_resp_fract_spines);
    mean_soma_unresp_fract_spines = mean(soma_unresp_fract_spines);
    sem_soma_resp_fract_spines = std(soma_resp_fract_spines)/sqrt(length(soma_resp_fract_spines));
    sem_soma_unresp_fract_spines = std(soma_resp_fract_spines)/sqrt(length(soma_unresp_fract_spines));
    
    subplot(3,1,i)
    bar([mean_soma_unresp_fract_spines,mean_soma_resp_fract_spines]); hold on;
    errorbar([mean_soma_unresp_fract_spines,mean_soma_resp_fract_spines],[sem_soma_unresp_fract_spines, sem_soma_resp_fract_spines], '.k');
    scatter(ones(size(soma_unresp_fract_spines)), soma_unresp_fract_spines, '.k','jitter', 'on', 'jitterAmount', 0.05);
    scatter(2*ones(size(soma_resp_fract_spines)), soma_resp_fract_spines,'.k', 'jitter', 'on', 'jitterAmount', 0.05);
    title(session{i})
    xticklabels({'soma lose response', 'soma maintain response'})
    ylabel('\Delta % responsive spines')
    [p,a] = ranksum(soma_resp_fract_spines, soma_unresp_fract_spines)

  

end

%%
function [fraction_resp_by_mouse_by_day,soma_resp_by_mouse_by_day] = get_eye_resp_mouse_day(unique_mouse,unique_days, contra_ipsi_binoc, all_resp,resp_type)
    fraction_resp_by_mouse_by_day = NaN(length(unique_mouse), length(unique_days));
    soma_resp_by_mouse_by_day = NaN(length(unique_mouse), length(unique_days));
    for i = 1:length(unique_mouse)
         for iii = 1:length(unique_days)
             temp = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.mouse_cell, unique_mouse{i})&...
                 strcmp(contra_ipsi_binoc.days, unique_days{iii}),:);
             if ~isempty(temp)
                 if all_resp
                     resp = height(temp(~strcmp(temp.cib_spine, resp_type),:));
                     soma_resp = max(~strcmp(temp.cib_soma, resp_type));
                 else
                     resp = height(temp(temp.(resp_type)==1,:));
                     soma_resp = max(temp.(['soma_',resp_type]));
                 end
                 fraction_resp_by_mouse_by_day(i,iii) = resp*100/height(temp);
                 soma_resp_by_mouse_by_day(i,iii) = soma_resp;
             end
         end
    end
end
%%
function [contra_ipsi_binoc] = concate_contra_ipsi_binoc_v2(temp, features, join_key)
    if ~isempty(find(strcmp([temp.session(:)], 'contra')))
        contra = temp(strcmp([temp.session(:)], 'contra'),features);
        ipsi = temp(strcmp([temp.session(:)], 'ipsi'),features);
        binoc = temp(strcmp([temp.session(:)], 'binoc'),features);
    else
        contra = temp(strcmp([temp.session{:}], "contra"),features);
        ipsi = temp(strcmp([temp.session{:}], "ipsi"),features);
        binoc = temp(strcmp([temp.session{:}], "binoc"),features);
    end
    contra_ipsi = outerjoin(contra, ipsi, "Keys",join_key, "MergeKeys",true);
    contra_ipsi_binoc = outerjoin(contra_ipsi, binoc, "Keys",join_key, "MergeKeys",true);
    
end
