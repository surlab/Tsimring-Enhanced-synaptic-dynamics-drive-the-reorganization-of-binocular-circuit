function normalize_stimevents_spine_longstim(timestamps,ROIdata_scan,graph,files,framelength)
neuropil_sub_flag = false;
Shaft_corr_flag = true;

%% subtract neuropil from corresponding spine and reduce size of matrix.
[numrows,numcolumns] = size(ROIdata_scan);
if neuropil_sub_flag == true
    disp('subtracting neuropil...');
    ROIdata = NaN(numrows,numcolumns-numcolumns/3);
    ROIdata(:,1:2:end) = ROIdata_scan(:,1:3:end) - ROIdata_scan(:,3:3:end); %subtract background from spine
    ROIdata(:,2:2:end) = ROIdata_scan(:,2:3:end); %keep shafts the same
    %clear ROIdata_scan
else
    ROIdata_scan(:,3:3:end) = [];
    ROIdata = ROIdata_scan;
    %clear ROIdata_scan
end

%set up shaft matrix. Every shaft labeled in column 1 by a 1, every spine
%has the index of the column of its corresponding shaft in column 2. 
shaft = zeros(numcolumns-numcolumns/3,2);
shaft(2:2:end,1) = 1;
shaft(1:2:end,2) = [2:2:numcolumns-numcolumns/3];

%% .mat file that has the stimulus order saved. 
load(fullfile(files(1).folder, 'stim', 'stim.mat')); %Your stim file from disp_grating or disp_movie
if ~isempty(stim)
[~,~,num_dims] = size(stim);
stimevents = stim;
%remove redundant dimensions with only one variable
if num_dims>1
    numunique = arrayfun(@(x) length(unique(stimevents(:,:,x))), [1:num_dims]); 
    numunique_catch = ones(length(numunique),1);
    numunique_catch(numunique>1) = 0;
    stimevents(:,:,numunique_catch == 1) = [];
    numunique(numunique_catch == 1) = [];
    clear uniquemm numunique_catch
    else
    numunique = length(unique(stimevents));
end
[~,~,num_dims] = size(stimevents);
  
%find # stimulations total
stimnumber=numel(stimevents(:,:,1));

%Resphape matrix into 2D ordered from 1st row to last row by time.
stim_sets = arrayfun(@(x) reshape(stimevents(:,:,x)',1,stimnumber), [1:num_dims], 'UniformOutput', false);
stim_setsmat = zeros(stimnumber,num_dims);
for mm = 1:num_dims
    stim_setsmat(:,mm) = stim_sets{mm};
end
clear mm

%set KT_flag if Katya stimuli used
KT_flag = ismember(888,stim_setsmat);

%% Calculate Ephys
voltage_file = files(contains({files.name}, '.csv') & contains({files.name}, 'VoltageRecording'));
file = fullfile(voltage_file(1).folder, voltage_file(1).name);
electrophystime = csvread(file,2);
   
%round to 0 or 5, whichever voltage is closest
roundedtimes=electrophystime;
round_targets = [0,5];
[~,Index1] = histc(electrophystime(:,2),[-Inf interp1(1:numel(round_targets),round_targets,0.5 + (1:numel(round_targets)-1)) Inf]);
roundedtimes(:,2) = round_targets(Index1);
%take difference roundedtimes(2,2) - roundedtime(1,2)...roundedtimes(n+1,2) - roundedtime(n,2)
diff_voltage = diff(roundedtimes(:,2));
%add 0 to start to get same length as timestamps
diff_voltage = [0;diff_voltage];
%5 is start of a square wave, -5 is the end of a square wave. Take these
%timestamps. Divide by 1000 to be seconds. 
startevent = roundedtimes(diff_voltage==5,1)/1000;
blankframes = roundedtimes(diff_voltage==-5,1)/1000;
%If first event is less than 0.8 seconds, remove and assume it is a pupil
%trigger. Inform user. 
if length(blankframes) == length(startevent)
if blankframes(1) - startevent(1) <0.8
    startevent = startevent(2:end);
    blankframes = blankframes(2:end);
    disp('First event excluded. If this was not a trigger event comment out the if statement in the code');
end

%Calculate stim duration and interstimulus interval duration
stim_dur = mean(arrayfun(@(x) blankframes(x) - startevent(x), 1:length(startevent)-1));
isi_dur = mean(arrayfun(@(x)startevent(x+1) - blankframes(x), 1:length(startevent)-1));   

%If KT_flag==1, remove trials longer than 2 seconds
if KT_flag == 1
    disp('KT flag active, removing KT trials');
    KT_trials = stim_setsmat==888; %find KT trials
    stim_setsmat(KT_trials,:) = []; %remove from stims
    numunique = numunique-1;
    %Calculate stim duration and interstimulus interval duration
    stim_dur = mode(arrayfun(@(x) blankframes(x) - startevent(x), 1:length(startevent)-1));
    isi_dur = mode(arrayfun(@(x)startevent(x+1) - blankframes(x), 1:length(startevent)-1)); 
    startevent(KT_trials) = []; %remove from startevents
else
    %Calculate stim duration and interstimulus interval duration
    stim_dur = mean(arrayfun(@(x) blankframes(x) - startevent(x), 1:length(startevent)-1));
    isi_dur = mean(arrayfun(@(x)startevent(x+1) - blankframes(x), 1:length(startevent)-1));        
end
clear mm

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

%% Identify trials with motion artifacts and bAPs
avg_shaft = mean(ROIdata(:,logical(shaft(:,1))),2); %average shaft ROIs
zscore_avg_shaft = zscore(avg_shaft);
%Report 1 if trial has a value below -2
artifact_trial = arrayfun(@(x) max(zscore_avg_shaft(timestamps>=startevent(x) & timestamps<=startevent(x)+stim_dur)<=-2), 1:length(startevent));
bAP_trial = arrayfun(@(x) mean(zscore_avg_shaft(timestamps>=startevent(x) & timestamps<=startevent(x)+stim_dur))>=2, 1:length(startevent));

%% Subtract dendrite contamination
%identify dendritic shaft ROIs, and for each spine subtract its closest
%shaft, use robust regression to subtract shaft signal, and calculate
%remaining shaft correlation
if Shaft_corr_flag==true
    [ROIdata_fit,shaft_corr] = sub_chen_KT(ROIdata_F, framelength, numROIs,shaft);
    [~,numROIs_fit] = size(ROIdata_fit);
    clear ROIdata_F numROIs
else
    ROIdata_fit = ROIdata_F(:,1:2:end);
    [~,numROIs_fit] = size(ROIdata_fit);
    clear ROIdata_F numROI
end
 
%% Recollect ISIs and calculate Z-Score
%preallocate matrix for blank frame data (period b/t stims, ISI)
F_data=zeros(1,(numROIs_fit));
%loop to collect ISI 
for ii=1:length(blankframes)
    Chunk = (ROIdata_fit(timestamps(:,1) >=(blankframes(ii)+isi_dur/2)&timestamps(:,1) <=(blankframes(ii)+isi_dur),:));
    F_data = [F_data; Chunk];
end
%exclude first empty row
F_data=F_data(2:end,:);
%calculate stdev
stdevF=nanstd(F_data(1:end,:),[],1);
%calculate mean of each ROI's ISIs
meanF = nanmean(F_data(1:end,:),1);
clear F_data Chunk ii
%Z-score
ROIdata_Z = zeros(TimeLength,numROIs_fit);
for ii=1:numROIs_fit
ROIdata_Z(:,ii)=(ROIdata_fit(:,ii)-meanF(1,ii))./stdevF(1,ii);
end
clear ii ROIdata_fit

% if there are more stims than startevents, calculate startevents based on
% isi and stim durations
if length(startevent)<length(stim_setsmat)
    for ii = length(startevent):length(stim_setsmat)
        startevent(ii) = startevent(ii-1) + isi_dur + stim_dur;
    end
    disp('more stims reported than show up in electrophys recording. This may be because the electrophys recording cut out early');
end
%append startevent and artifact_trials as final columns
stim_pairs = [stim_setsmat,startevent,artifact_trial', bAP_trial'];
%stim_pairs = [stim_setsmat,startevent,artifact_trial'];
clear stim_sets stim_setsmat

%Sort by values of 1st dimension
stim_pairs = sortrows(stim_pairs);
%take artifact_trial out
artifact_trial = stim_pairs(:,end-1);
bAP_trial = stim_pairs(:,end);
stim_pairs = stim_pairs(:,1:end-2);

%% extract time series data. 
%Because the number of dimensions cannot be known a priori, treat every
%unique value of the 1st dimension (1st column) as a category (usually
%orientation) and then every combination of the other dimensions as a
%variation of that category. For example, 90 degrees with 0.04 CPD, 90
%degrees with 0.1 CPD, etc. Organize cell output as 1st dim x variations x
%repetitions.

%find every combination of unique values from each column (excluding time)
[unique_rows,~,stim_index] = unique(stim_pairs(:,1:end-1), 'rows');
%Find the row indices of each unique combination
indices = accumarray(stim_index, find(stim_index), [], @(rows){rows});
%find unique orientations/1st dimension
unique_orientations = unique(stim_pairs(:,1));

n = 0;
for ii = 1:numunique(1) % number of unique orietantions or first dimensions
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
clear ii jj kk n

%% Calculate active spine events smoothed trace over entire recording 
% Bin active spine events by stim pairs
max_index = max(arrayfun(@(x) length(indices{x}), [1:length(indices)]));
max_dim = max(arrayfun(@(x) sum(unique_rows(:,1)==unique_orientations(x)), [1:numunique(1)]));

peak_threshold = 3;
conseq_frames = 3;
alpha = 0.4;
%binary matrix specifying whether spine was active or not during trial
%(peak above 3 std)
active_spine_trial_smooth = zeros(numROIs_fit, length(unique_orientations), max_dim, max_index);
active_spine_trial = zeros(numROIs_fit, length(unique_orientations), max_dim, max_index);

for l = 1:numROIs_fit
    n = 0;
    stim_trace = ROIdata_Z(:,l);
    smoothed_stim_trace = exp_weight_filter(ROIdata_Z(:,l),alpha);
    [pks, lks] = findpeaks(smoothed_stim_trace, 'MinPeakHeight',peak_threshold,'MinPeakDistance',conseq_frames);
    pk_times_smooth = timestamps(lks);
    [pks, lks] = findpeaks(stim_trace, 'MinPeakHeight',peak_threshold,'MinPeakDistance',conseq_frames);
    pk_times = timestamps(lks);
    for ii = 1:numunique(1) % number of unique orietantions or first dimensions
         for jj = 1:sum(unique_rows==unique_orientations(ii)) % number of unique stims for that orientation/first dimension
            n=n+1;  
            for kk = 1:length(indices{n})
                start_onset = stim_pairs(indices{n}(kk),end);
                start_offset = start_onset+stim_dur+0.5;
                inds = find(pk_times_smooth>=start_onset& pk_times_smooth<=start_offset);
                if ~isempty(inds)
                    active_spine_trial_smooth(l,ii,jj,kk) = 1;
                end

                inds = find( pk_times>=start_onset& pk_times<=start_offset-0.5);
                if ~isempty(inds)
                    active_spine_trial(l,ii,jj,kk) = 1;
                end
             end
        end
    end
end


%% Calculate mean response and p values


%reshape artifact_trials into trial x orientation matrix to preserve
%correct order
artifact_trial = reshape(artifact_trial,[max_index,length(unique_orientations)]);
%reshape bAP_trials into trial x orientation matrix to preserve
%correct order
bAP_trial = reshape(bAP_trial,[max_index,length(unique_orientations)]);
%mean of 0 to 1.5 post stim
stim_mean = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);
%mean of pre period -1.5 to 0 seconds
pre_mean = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);
%preallocate matrix for each individual pre trial period
pre_peak=zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);
%Amplitude of trials
trial_amp = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);
%preallocate matrix for trial amplitude, which is trial_data - pre_data;
trial_peak = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);
%Preallocate matrix for peak response to stimuli
peak_data = zeros(numROIs_fit,length(unique_orientations),max_dim);
%Z-test across trials
z_test = nan(numROIs_fit,length(unique_orientations),max_dim);
%T test across trial
t_test = nan(numROIs_fit,length(unique_orientations),max_dim);
%T test across trial without bAPs removed
t_test_bAPs_removed = nan(numROIs_fit,length(unique_orientations),max_dim);
%pre means above 2 STD
large_pre_mean = zeros(numROIs_fit,length(unique_orientations),max_dim);
%Take median of trial amplitudes
median_amplitude =  nan(numROIs_fit,length(unique_orientations),max_dim);
%Take mean of trial amplitudes
mean_amplitude =  nan(numROIs_fit,length(unique_orientations),max_dim);
%Take mean of trial amplitudes with bAPs removed
mean_amplitude_bAPs_removed =  nan(numROIs_fit,length(unique_orientations),max_dim);
%Take std of trial amplitudes
std_amplitude =  nan(numROIs_fit,length(unique_orientations),max_dim);
%Take std of trial amplitudes with bAPs removed
std_amplitude_bAPs_removed =  nan(numROIs_fit,length(unique_orientations),max_dim);
%Fano factor
fano_factor = nan(numROIs_fit,length(unique_orientations),max_dim);
%Reliability index
reliability_index = nan(numROIs_fit,length(unique_orientations),max_dim);
%preallocate matrix with trials that are included
included_trials = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);
%preallocate matrix with trials that are included counting bAPs
included_trials_bAP = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index);

%traces of individual trials
trial_time = (round(-isi_dur*10):1:round(stim_dur*10+isi_dur*10))/10;
trial_traces = zeros(numROIs_fit,length(unique_orientations),max_dim,max_index,length(trial_time));



%Loop through the ROIs
for l=1:numROIs_fit
    n = 0;
    if l ==15
        ts = 1;
    end
    for ii = 1:length(unique_orientations)
        for jj = 1:1:sum(unique_rows==unique_orientations(ii)) 
            n = n+1;
            %create empty matrix to contain interpolated stim traces
            for kk = 1:length(indices{n})
                %interpolate each trial on a 0.1 S scale  
                trial_traces(l,ii,jj,kk,:)=interp1(Datacell{ii,jj,kk}(:,1),Datacell{ii,jj,kk}(:,l+1),(round(-isi_dur*10):1:round(stim_dur*10+isi_dur*10))/10);
                stim_trace = Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=0&Datacell{ii,jj,kk}(:,1)<=1.5,l+1);
                pre_stim_trace = Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=-1.5 & Datacell{ii,jj,kk}(:,1)<0,l+1);                
                %Take mean of 0 to 1.5 seconds post stim
                stim_mean(l,ii,jj,kk) = mean(stim_trace);
                %take mean of pre-stim periods -1.5 to 0
                pre_mean(l,ii,jj,kk) = mean(pre_stim_trace);
                %take peak of pre-stim periods 
                pre_peak(l,ii,jj,kk)=max(pre_stim_trace);
                %take amplitude of individual trials over pre periods
                trial_amp(l,ii,jj,kk) = stim_mean(l,ii,jj,kk) - pre_mean(l,ii,jj,kk);
                %find peak of individual traces
                trial_peak(l,ii,jj,kk)=max(stim_trace);

            end
            %logical array for trials with pre_mean above 2 std
            large_pre_mean(l,ii,jj) = sum(permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) >= peak_threshold);
          
            %Calculate excluding pre periods above 2 STD and artifact
            %trials
            include_trials = permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) < peak_threshold & artifact_trial(:,ii) == 0;
            include_trials_bAP = permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) < peak_threshold & artifact_trial(:,ii) == 0 & bAP_trial(:,ii)==0;

            %Take Z-test of included trials
            [~,z_test(l,ii,jj)] = ztest(stim_mean(l,ii,jj,include_trials),mean(pre_mean(l,ii,jj,include_trials)),std(pre_mean(l,ii,jj,include_trials)),'Tail','right');
            %Take T-test of included trials
            [~,t_test(l,ii,jj)] = ttest(stim_mean(l,ii,jj,include_trials),pre_mean(l,ii,jj,include_trials),'Tail','right');
            %Take T-test of included trials with bAP removed
            [~,t_test_bAPs_removed(l,ii,jj)] = ttest(stim_mean(l,ii,jj,include_trials_bAP),pre_mean(l,ii,jj,include_trials_bAP),'Tail','right');

            %Take median of included trials
            median_amplitude(l,ii,jj) = median(trial_amp(l,ii,jj,include_trials));
            %Take mean of included trials
            mean_amplitude(l,ii,jj) = mean(trial_amp(l,ii,jj,include_trials));
            %Take mean of included trials with bAPs removed
            mean_amplitude_bAPs_removed(l,ii,jj) = mean(trial_amp(l,ii,jj,include_trials_bAP));

            %Take std of included trials
            std_amplitude(l,ii,jj) = std(trial_amp(l,ii,jj,include_trials));
            %Take std of included trials with bAPs removed
            std_amplitude_bAPs_removed(l,ii,jj) = std(trial_amp(l,ii,jj,include_trials_bAP));
          
            %take fano factor of included trials
            fano_factor(l,ii,jj) = (std_amplitude(l,ii,jj))/mean_amplitude(l,ii,jj);
            %Take reliability index of included trials 0 to 1.5 seconds
            corr_trials = corr(permute(trial_traces(l,ii,jj,include_trials,trial_time>=0 & trial_time<=1.5),[5,4,3,2,1])); %take correlation of interpolated trials
            idx = logical(triu(ones(size(corr_trials)),1)); %Make logical index of upper triangle of corr matrix 
            reliability_index(l,ii,jj) = (2/(sum(include_trials)^2 - sum(include_trials))) * sum(corr_trials(idx)); %take average of upper diagonal and divide by n to get mean   
            
            included_trials(l,ii,jj,:) =  include_trials';
            included_trials_bAP(l,ii,jj,:) = include_trials_bAP';
            
            clear include_trials corr_trials idx
        end
    end
end
clear ii jj kk


if graph==1

time_axis = -round(isi_dur):.1:round(stim_dur+isi_dur);
if isfolder(fullfile(files(1).folder, 'figs'))
   delete(fullfile(files(1).folder, 'figs','*.fig'))
end
%Loop through the ROIs
for l=1:numROIs_fit
    responsive = 'unresponsive';
    index_isi = find(time_axis < 0 & time_axis >= -isi_dur/3);
    index_stim = find(time_axis >=0 & time_axis <= 1);
    %time_length = numel(time_axis);
    n=0; 
    %create new figure with ROI name and position. Adjust for your monitor
    %needs. 
    figure('name',sprintf('Plot of ROI %d',l),'numbertitle','off','position',[414,391,472,257]);
    %loop through stims
    time_period = time_axis + isi_dur;
   
    for ii = 1:length(unique_orientations)
        for jj = 1:1:sum(unique_rows(:,1)==unique_orientations(ii)) 
            
            
            n = n+1;  
            %set colors to use for plots
            if contains(filename, 'contra')
                c = 'r';
            elseif contains(filename, 'ipsi')
                c = 'b';
            else
                c = 'y';
            end
            graphsize=ceil(length(unique_orientations)/2);    
           
            %create empty matrix to contain interpolated stim traces
            stim_traces=zeros(round(20*isi_dur + 10*stim_dur),length(indices{n}));
            for kk = 1:length(indices{n})
                %interpolate each trial on a 0.1 S scale  
                stim_traces(:,kk)=interp1(Datacell{ii,jj,kk}(:,1),Datacell{ii,jj,kk}(:,l+1),(round(-isi_dur*10):1:floor(stim_dur*10+isi_dur*10))/10);
                
                %plot stim trace
                plot(time_period([index_isi,index_stim]),stim_traces([index_isi,index_stim],kk),'Color',[0.5,0.5,0.5],'LineWidth',0.25);
                
                if kk ==1
                    hold on          
                    xlabel('Time (s)')
                    ylabel('Z score')
                end
            end
            
            include_trials = permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) < peak_threshold & artifact_trial(:,ii) == 0;
            %Take mean of traces 
            mean_trace=nanmean(stim_traces(:,include_trials),2);
            
            %plot mean trace
            plot(time_period([index_isi,index_stim]),mean_trace([index_isi,index_stim]),'Color',c,'LineWidth',1);
           
                       
            clear mean_trace TF stim_traces trial_traces
        end
        %label subplots with sig responses
        if min(t_test(l,ii,:))<=0.05 && max(mean_amplitude(l,ii,:))>=0.5 && sum(include_trials) >= 5
            patch([time_period(index_stim(1)),time_period(index_stim(end)),time_period(index_stim(end)),time_period(index_stim(1))], [-5, -5,5, 5],c, 'FaceAlpha', 0.3, 'LineStyle', 'none')    
            responsive = 'responsive';
        else
            patch([time_period(index_stim(1)),time_period(index_stim(end)),time_period(index_stim(end)),time_period(index_stim(1))], [-5, -5,5, 5],[0.5,0.5,0.5], 'FaceAlpha', 0.3,'LineStyle', 'none')          
        end
        %hold off
        time_period = time_period + round(stim_dur*10+(isi_dur/2)*10)/10;
    end
    
    

    %pause then advance. Allows user time to look at the plots. 
    pause on
    pause(1);
    
    ylim([-5,5])
    pause off
    
    if ~isfolder(fullfile(files(1).folder, 'figs'))
        mkdir(fullfile(files(1).folder, 'figs'))
    end
     
    savefig(fullfile(files(1).folder, 'figs', ['Plot of ', responsive, ' ROI ', num2str(l), '.fig'])) 
    close all
end
%% concatenate each successive tiff to tiff_stack
ori_num = 1;
trial_num = 1;
count = 1;
for l = 1: size(ROIdata_Z,2)
    roidata = ROIdata_Z(:,l)
    f = figure
    set(gca, 'FontSize', 20)
    start_ind = 1;
    
    plot(timestamps,roidata,'Color','k', 'LineWidth',3); hold on;
           

     for ii = 1:length(startevent)
        patch([startevent(ii),blankframes(ii),blankframes(ii),startevent(ii)], [-5, -5,12, 12],[0.5,0.5,0.5], 'FaceAlpha', 0.3,'LineStyle', 'none') 
        %patch([startevent(ii),blankframes(ii),blankframes(ii),startevent(ii)], [0, 0,2000, 2000],[0.5,0.5,0.5], 'FaceAlpha', 0.3,'LineStyle', 'none') 
    
     
       
        %axis([0 30000 0 20])
       
       
     end
     xlabel('Time(s)');
     ylabel('F (z scored)')
     savefig(fullfile(files(1).folder, 'figs', ['Full trace of  ROI ', num2str(l), '.fig'])) 
     close(f);
end
end
clear ii jj kk  
%Save Data in the current directory
save(fullfile(files(1).folder, 'normalized_data_by_stim'));

end
end
end
%% for spine calcium events
function smooth_data = exp_weight_filter(data,alpha)
    smooth_data = zeros(size(data));
    smooth_data(1) = data(1);
    for t = 2:length(data)
        smooth_data(t) = alpha*data(t) + (1-alpha)*smooth_data(t-1);
    end
end

