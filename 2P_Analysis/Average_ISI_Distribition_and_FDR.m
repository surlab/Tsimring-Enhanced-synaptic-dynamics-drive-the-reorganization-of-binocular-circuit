%% Noise Distribition and FDR analysis for spines
% 1. Calculate the distribution of the trial-averaged ISI periods for spines across
% all sessions, days, and cells to determine the 95% percentile
% 2. Calculate the FDR rate using the mean amplitiude and t-test criteria

mouse_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Processed Data/';
%mouse_path = 'G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Processed Data/';
mouse_files = { 'BM014','BM015', 'BM016', 'BM017', 'BM018', 'BM019', 'BM020', 'BM021', 'BM023','BM024', 'BM025', 'BM026', 'BM027','BM029', 'BM030'};

%% Find the noise distribition of the trial-averaged ISI periods for all pooled spines
all_ISI = [];
save_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';
save_file = [ 'soma_ISI_distribtion',strjoin(mouse_files, '_'),'.mat'];
figure
for i = 1:length(mouse_files)
    ISI = []
    all_norm_files = dir(fullfile(mouse_path,mouse_files{i},'**/normalized_data_by_stim.mat'));
    all_norm_files = all_norm_files(contains({all_norm_files.folder}, 'Soma'));
    for ii = 1:length(all_norm_files)
        load(fullfile(all_norm_files(ii).folder, all_norm_files(ii).name),'pre_mean');
        mean_pre_mean = mean(pre_mean,4);
        ISI = [ISI; mean_pre_mean(:)]; %find the ISI by mouse to see if there's any differences
    end
    title(mouse_files{i})
    histogram(ISI), shg
    all_ISI = [all_ISI; ISI];
end
y=prctile(all_ISI,99.8)
save(fullfile(save_path,save_file), 'y', 'all_ISI');
%% FDR analysis for spines: permutate among the start event times
% this keeps the structure of the trials the same
% 01/23/24 update: this doesn't seem to actually work (need to figure out
% the math...)

fraction_resp = [];
boot_strap = 10;
mean_amp = 0.5;
p_thresh = 0.05;
all_mean_amplitude = [];
all_t_test = [];
all_large_pre_mean = [];
all_artifact_trial = [];
for i = 1:length(mouse_files)
    all_norm_files = dir(fullfile(mouse_path,mouse_files{i},'**/normalized_data_by_stim.mat'));
    all_norm_files = all_norm_files(~contains({all_norm_files.folder}, 'Soma'));
    for ii = 1:length(all_norm_files)
        load(fullfile(all_norm_files(ii).folder, all_norm_files(ii).name),'zscore_avg_shaft','ROIdata_Z', 'startevent',...
            'stimevents','timestamps', 'max_dim','max_index', 'stim_dur', 'isi_dur', 'numROIs_fit', 'stim_pairs');
        [~,~,num_dims] = size(stimevents);
        if num_dims>1
            % this is for the multiple sfs
        else
            reshape_events = stimevents';
            reshape_events = reshape_events(:);
        end
        artifact_trial = create_artifact_trials(zscore_avg_shaft, startevent, timestamps,stim_dur);
        perm_inds = arrayfun(@(x) randperm(length(startevent)), [1:boot_strap]', 'UniformOutput',false);
        all_startevents = cellfun(@(x) startevent(x), perm_inds, 'UniformOutput',false);
        all_artifact_trial_shuffle = cellfun(@(x) artifact_trial(x), perm_inds, 'UniformOutput',false);
        unique_orientations = unique(reshape_events);
        all_mean_amplitude_by_shuffle = [];
        all_t_test_by_shuffle = [];
        all_large_pre_mean_by_shuffle = [];
        all_artifact_trial_by_shuffle = [];
        for s = 1:length(all_startevents)
            artifact_trial_temp = reshape(all_artifact_trial_shuffle{s},[max_index,length(unique_orientations)]);
            
            stim_pairs_temp = [reshape_events,all_startevents{s}];  
            stim_pairs_temp = sortrows(stim_pairs_temp);
            [Datacell, indices, unique_rows] = stim_align_traces(stim_pairs_temp,ROIdata_Z, timestamps,stim_dur,isi_dur);
            [mean_amplitude, t_test,large_pre_mean] = mean_response_p_values(Datacell,indices,...
                artifact_trial_temp, numROIs_fit, max_dim, max_index, unique_orientations, unique_rows);

            total_artifacts = sum(artifact_trial_temp,1);
            all_mean_amplitude_by_shuffle = [all_mean_amplitude_by_shuffle, {mean_amplitude}];
            all_t_test_by_shuffle = [all_t_test_by_shuffle, {t_test}];
            all_large_pre_mean_by_shuffle = [all_large_pre_mean_by_shuffle, {large_pre_mean}];
            all_artifact_trial_by_shuffle = [all_artifact_trial_by_shuffle , {repmat(total_artifacts,size(t_test,1),1)}];
        end
        all_mean_amplitude = [all_mean_amplitude;all_mean_amplitude_by_shuffle];
        all_t_test = [all_t_test; all_t_test_by_shuffle ];
        all_large_pre_mean = [all_large_pre_mean;all_large_pre_mean_by_shuffle];
        all_artifact_trial = [all_artifact_trial; all_artifact_trial_by_shuffle];
    end
    
end
%% FDR analysis for spines: randomly pick 80 start times from timestamps
boot_strap = 1000;
all_mean_amplitude = [];
all_t_test = [];
all_large_pre_mean = [];
all_artifact_trial = [];
for i = 1:length(mouse_files)
    all_norm_files = dir(fullfile(mouse_path,mouse_files{i},'**/normalized_data_by_stim.mat'));
    all_norm_files = all_norm_files(~contains({all_norm_files.folder}, 'Soma'));
    for ii = 1:length(all_norm_files)
        load(fullfile(all_norm_files(ii).folder, all_norm_files(ii).name),'zscore_avg_shaft','ROIdata_Z',...
            'stimevents','timestamps', 'max_dim','max_index', 'stim_dur', 'isi_dur', 'numROIs_fit', 'stim_pairs');
        [~,~,num_dims] = size(stimevents);
        if num_dims>1
            % this is for the multiple sfs
        else
            reshape_events = stimevents';
            reshape_events = reshape_events(:);
        end
        unique_orientations = unique(reshape_events);
        last_timepoint = timestamps(end);
        first_timepoint = timestamps(1);
        trim_timestamps = timestamps(timestamps < (last_timepoint-(isi_dur+stim_dur)) & ...
            timestamps > (first_timepoint+(isi_dur+stim_dur)));
        all_mean_amplitude_by_shuffle = [];
        all_t_test_by_shuffle = [];
        all_large_pre_mean_by_shuffle = [];
        all_artifact_trial_by_shuffle = [];
        for s = 1:boot_strap
            
            startevent_temp = trim_timestamps(randperm(length(trim_timestamps),length(reshape_events)));
            artifact_trial_temp = create_artifact_trials(zscore_avg_shaft, startevent_temp, timestamps,stim_dur);
            artifact_trial_temp = reshape(artifact_trial_temp ,[max_index,length(unique_orientations)]);
           
            stim_pairs_temp = [reshape_events,startevent_temp];  
            stim_pairs_temp = sortrows(stim_pairs_temp);
            [Datacell, indices, unique_rows] = stim_align_traces(stim_pairs_temp,ROIdata_Z, timestamps,stim_dur,isi_dur);
            [mean_amplitude, t_test,large_pre_mean] = mean_response_p_values(Datacell,indices,...
                artifact_trial_temp, numROIs_fit, max_dim, max_index, unique_orientations, unique_rows);

            total_artifacts = sum(artifact_trial_temp,1);
            all_mean_amplitude_by_shuffle = [all_mean_amplitude_by_shuffle, {mean_amplitude}];
            all_t_test_by_shuffle = [all_t_test_by_shuffle, {t_test}];
            all_large_pre_mean_by_shuffle = [all_large_pre_mean_by_shuffle, {large_pre_mean}];
            all_artifact_trial_by_shuffle = [all_artifact_trial_by_shuffle , {repmat(total_artifacts,size(t_test,1),1)}];
        end
        all_mean_amplitude = [all_mean_amplitude;all_mean_amplitude_by_shuffle];
        all_t_test = [all_t_test; all_t_test_by_shuffle ];
        all_large_pre_mean = [all_large_pre_mean;all_large_pre_mean_by_shuffle];
        all_artifact_trial = [all_artifact_trial; all_artifact_trial_by_shuffle];
    end
    
end
%% FDR analysis for spines: find the fraction of responsive spines per permutation
% check the FDR for different mean amplitude and p value thresholds
save_path = 'G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';

mouse_files = { 'BM014','BM015', 'BM016', 'BM017', 'BM018', 'BM019', 'BM020', 'BM021', 'BM023','BM024', 'BM025', 'BM026', 'BM027','BM029', 'BM030'};
save_file = [ 'boot_strap_FDR',strjoin(mouse_files, '_'),'.mat'];
boot_strap = 1000;
load(fullfile(save_path,save_file))
mean_amp_thresh = [0.5,0.75,1];
p_thresh = [0.01,0.05];
func = @(m,p,a,b,c,d) max(int8(a>m) + int8(b<p) +int8((c+d)<=5));
all_m = [];
for m = mean_amp_thresh
    all_p =[];
    for p = p_thresh
        all_resp = [];
        for s = 1:boot_strap
            artifact_trial_temp = cell2mat(all_artifact_trial(:,s));
            large_pre_mean = cell2mat(all_large_pre_mean(:,s));
            mean_amp = cell2mat(all_mean_amplitude(:,s));
            t_test = cell2mat(all_t_test(:,s));
            criteria = [];
            for roi = 1:size(t_test,1)
                criteria = [criteria;func(m,p,mean_amp(roi,:),t_test(roi,:),large_pre_mean(roi,:), artifact_trial_temp(roi,:))];
            end
            resp = criteria==3;
            all_resp = [all_resp,resp];

        end
        all_p = [all_p,{all_resp}];
    end
    all_m = [all_m; {all_p}];
end
%% FDR plot
% check the FDR for different mean amplitude and p value thresholds
save_path = 'G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';
save_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';

mouse_files = { 'BM014','BM015', 'BM016', 'BM017', 'BM018', 'BM019', 'BM020', 'BM021', 'BM023','BM024', 'BM025', 'BM026', 'BM027','BM029', 'BM030'};
save_file = [ 'boot_strap_FDR',strjoin(mouse_files, '_'),'.mat'];
boot_strap = 1000;
load(fullfile(save_path,save_file))
%%
mean_amp_thresh = [0.5];
p_thresh = [0.01];
func = @(m,p,a,b,c,d) int8(a>m) + int8(b<p) +int8((c+d)<=5);
all_m = [];
for m = mean_amp_thresh
    all_p =[];
    for p = p_thresh
        all_resp = [];
        for s = 1:boot_strap
            artifact_trial_temp = cell2mat(all_artifact_trial(:,s));
            large_pre_mean = cell2mat(all_large_pre_mean(:,s));
            mean_amp = cell2mat(all_mean_amplitude(:,s));
            t_test = cell2mat(all_t_test(:,s));
            criteria = [];
            for roi = 1:size(t_test,1)
                ind = find(mean_amp(roi,:)==max(mean_amp(roi,:)));
                
                
                criteria = [criteria;func(m,p,mean_amp(roi,ind),t_test(roi,ind),large_pre_mean(roi,ind), artifact_trial_temp(roi,ind))];
            end
            resp = criteria==3;
            all_resp = [all_resp,resp];

        end
        all_p = [all_p,{all_resp}];
    end
    all_m = [all_m; {all_p}];
end
%% Plot swarmchart for all mean amplitude and p values
count = 1;
figure
mean_amp = [0.5,0.75,1];
p_thresh = [0.01,0.05];
for m = 1:length(all_m)
    all_p = all_m{m,:};
    for p = 1:length(all_p)
        subplot(length(all_m),length(all_p), count)
        fraction_resp = nanmean(all_p{p},2);
        mean(fraction_resp)
        swarmchart(ones(size(fraction_resp)),fraction_resp);
        count = count + 1;
        ylim([0,1])
    end
end

        

%%
function artifact_trial = create_artifact_trials(zscore_avg_shaft, startevent, timestamps, stim_dur)
    artifact_trial = arrayfun(@(x) max(zscore_avg_shaft(timestamps>=startevent(x) & timestamps<=startevent(x)+stim_dur)<=-2), 1:length(startevent));
end


%% align traces to traces to unique stimulus trials
function [Datacell,indices, unique_rows] = stim_align_traces(stim_pairs,ROIdata_Z, timestamps,stim_dur,isi_dur)
    %find every combination of unique values from each column (excluding time)
    [unique_rows,~,stim_index] = unique(stim_pairs(:,1:end-1), 'rows');
    %Find the row indices of each unique combination
    indices = accumarray(stim_index, find(stim_index), [], @(rows){rows});
    %find unique orientations/1st dimension
    unique_orientations = unique(stim_pairs(:,1));
    numROIs_fit = size(ROIdata_Z,2);
    n = 0;
    for ii = 1:length(unique_orientations) % number of unique orietantions or first dimensions
        for jj = 1:sum(unique_rows==unique_orientations(ii)) % number of unique stims for that orientation/first dimension
            n=n+1;  
            for kk = 1:length(indices{n})
                
                Datacell{ii,jj,kk}=timestamps(timestamps >=(stim_pairs(indices{n}(kk),end)-isi_dur)&timestamps <=(stim_pairs(indices{n}(kk),end)+stim_dur+isi_dur));
                % ROI data for all ROIs for that time range around the stim
                Datacell{ii,jj,kk}(:,2:numROIs_fit+1)=ROIdata_Z(timestamps >=(stim_pairs(indices{n}(kk),end)-isi_dur)&timestamps <=(stim_pairs(indices{n}(kk),end)+stim_dur+isi_dur),:);
                %normalize times of each event to be on a scale with 0 as the event
                %by subtracting the stim time
                Datacell{ii,jj,kk}(:,1)=Datacell{ii,jj,kk}(:,1)-stim_pairs(indices{n}(kk),end);   
            end
        end
    end
end

%%
function [mean_amplitude, t_test,large_pre_mean] = mean_response_p_values(Datacell,indices, artifact_trial, numROIs_fit, max_dim,max_index, unique_orientations, unique_rows)
    %mean of 0 to 1.5 post stim
    stim_mean = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);
    %mean of pre period -1.5 to 0 seconds
    pre_mean = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);
    %preallocate matrix for each individual pre trial period
    pre_peak=zeros(numROIs_fit,length(unique_orientations),max_dim,max_index); 
    %T test across trial
    t_test = nan(numROIs_fit,length(unique_orientations),max_dim);
    %Take mean of trial amplitudes
    mean_amplitude =  nan(numROIs_fit,length(unique_orientations),max_dim);
    %pre means above 3 STD
    large_pre_mean = zeros(numROIs_fit,length(unique_orientations),max_dim);
    num_orientations = length(unique_orientations);
    peak_threshold = 3;
    for l=1:numROIs_fit
    n = 0;
        for ii = 1:num_orientations 
            for jj = 1:1:sum(unique_rows==unique_orientations(ii)) 
               n = n+1;
                %create empty matrix to contain interpolated stim traces
                for kk = 1:length(indices{n})
                    %interpolate each trial on a 0.1 S scale  
                    stim_trace = Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=0&Datacell{ii,jj,kk}(:,1)<=1.5,l+1);
                    pre_stim_trace = Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=-1.5 & Datacell{ii,jj,kk}(:,1)<0,l+1);                
                    %Take mean of 0 to 1.5 seconds post stim
                    stim_mean(l,ii,jj,kk) = mean(stim_trace);
                    %take mean of pre-stim periods -1.5 to 0
                    pre_mean(l,ii,jj,kk) = mean(pre_stim_trace);
                    %take peak of pre-stim periods 
                    pre_peak(l,ii,jj,kk)=max(pre_stim_trace);
                    %take amplitude of individual trials over pre periods
                    trial_amp(l,ii,jj,kk) = stim_mean(l,ii,kk) - pre_mean(l,ii,kk);
    
                end
                %logical array for trials with pre_mean above 2 std
                large_pre_mean(l,ii,jj) = sum(permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) >= peak_threshold);
              
                %Calculate excluding pre periods above 2 STD and artifact
                %trials
                include_trials = permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) < peak_threshold & artifact_trial(:,ii) == 0;
                mean_amplitude(l,ii,jj) = mean(trial_amp(l,ii,jj,include_trials));
                [~,t_test(l,ii,jj)] = ttest(stim_mean(l,ii,jj,include_trials),pre_mean(l,ii,jj,include_trials),'Tail','right');
            end
        end
    end
end
