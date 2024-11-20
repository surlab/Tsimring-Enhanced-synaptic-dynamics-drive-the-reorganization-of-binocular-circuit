function  normalize_stimevents_soma_longstim(timestamps,ROIdata_scan,graph,files,framelength)
%% Select Electrophys data
%Electrophys file used to calculate timing of stimulation events
%open file
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
   
%% find F for each ROI 

%perform background subtraction of ROIs
%every odd ROI is real signal, even ROI is background
ROIdata_scan = ROIdata_scan(:,1:2:end,:)-0.7*ROIdata_scan(:,2:2:end);

%calculate number of ROIs
[TimeLength,numROIs]=size(ROIdata_scan);
%preallocate matrix for blank frame data (period b/t stims, ISI)
F_data=zeros(1,(numROIs));
%loop to collect ISI 
for ii=1:length(blankframes)
    Chunk = (ROIdata_scan(timestamps(:,1) >=(blankframes(ii)+isi_dur/2)&timestamps(:,1) <=(blankframes(ii)+isi_dur),:));
    F_data=[F_data;Chunk];
end
%exclude first empty row
F_data=F_data(2:end,:);
clear ii 

%% Calculate Z score of F for each ROI
%calculate mean and std of each ROI's ISI
meanF=mean(F_data(1:end,:),1);
stdevF=std(F_data(1:end,:),[],1);

%Delta F/F
ROIdata_Z=zeros(TimeLength,numROIs);
for ii=1:numROIs
    ROIdata_Z(:,ii)=(ROIdata_scan(:,ii)-meanF(1,ii))./stdevF(1,ii);
end
clear ii 

%% process by event
%.mat file that has the stimulus order saved.

load(fullfile(files(1).folder, 'stim', 'stim.mat')); %Your stim file from disp_grating or disp_movie
if ~isempty(stim)
[~,~,num_dims] = size(stim);

%remove redundant dimensions with only one variable
if num_dims>1
    numunique = arrayfun(@(x) length(unique(stim(:,:,x))), [1:num_dims]); 
    numunique_catch = ones(length(numunique),1);
    numunique_catch(numunique>1) = 0;
    stim(:,:,numunique_catch == 1) = [];
    numunique(numunique_catch == 1) = [];
    clear uniquemm numunique_catch
    else
    numunique = length(unique(stim));
end
[~,~,num_dims] = size(stim);
  
%find # stimulations total
stimnumber=numel(stim(:,:,1));

%Resphape into matrix into 2D ordered from 1st row to last row by time.
stim_sets = arrayfun(@(x) reshape(stim(:,:,x)',1,stimnumber), [1:num_dims], 'UniformOutput', false);
stim_setsmat = zeros(stimnumber,num_dims);
for mm = 1:num_dims
    stim_setsmat(:,mm) = stim_sets{mm};
end
clear mm

%append startevent as final column
if length(startevent)<length(stim_setsmat) %If you have fewer voltage events than stim
    for ii = length(startevent):length(stim_setsmat)
        startevent(ii) = startevent(ii-1) + isi_dur + stim_dur;
    end
    disp('more stims reported than show up in electrophys recording. This may be because the electrophys recording cut out early');
end
stim_pairs = [stim_setsmat,startevent];
clear stim_sets stim_setsmat

%Sort by values of 1st dimension
stim_pairs = sortrows(stim_pairs);

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
for ii = 1:numunique(1) %# of unique orietantions or first dimensions
    for jj = 1:sum(unique_rows==unique_orientations(ii)) %# of unique stims for that orientation/first dimension
        n=n+1;  
        for kk = 1:length(indices{n})
            Datacell{ii,jj,kk}=timestamps(timestamps >=(stim_pairs(indices{n}(kk),end)-isi_dur)&timestamps <=(stim_pairs(indices{n}(kk),end)+stim_dur+isi_dur));
            % ROI data for all ROIs for that time range around the stim
            Datacell{ii,jj,kk}(:,2:numROIs+1)=ROIdata_Z(timestamps >=(stim_pairs(indices{n}(kk),end)-isi_dur)&timestamps <=(stim_pairs(indices{n}(kk),end)+stim_dur+isi_dur),:);
            %normalize times of each event to be on a scale with 0 as the event
            %by subtracting the stim time
            Datacell{ii,jj,kk}(:,1)=Datacell{ii,jj,kk}(:,1)-stim_pairs(indices{n}(kk),end);   
        end
    end
end
clear ii jj kk n
 
%% Calculate mean response and p values
pre_peak_thresh = 3;
max_index = max(arrayfun(@(x) length(indices{x}), [1:length(indices)]));
max_dim = max(arrayfun(@(x) sum(unique_rows(:,1)==unique_orientations(x)), [1:numunique(1)]));
time_axis = -round(isi_dur):.1:round(stim_dur+isi_dur);
trial_time = (round(-isi_dur*10):1:round(stim_dur*10+isi_dur*10))/10;
%Preallocate matrix for peak response to stimuli
peak_data = zeros(numROIs,length(unique_orientations),max_dim);
%Preallocate matrix for mean response to stimuli
mean_data=zeros(numROIs,length(unique_orientations),max_dim);
%preallocate matrix for the mean pre-stim period
mean_pre=zeros(numROIs,length(unique_orientations),max_dim);
%preallocate matrix for each individual pre trial period
pre_mean=zeros(numROIs,length(unique_orientations),max_dim,max_index);
%preallocate matrix for each individual pre trial period
pre_peak=zeros(numROIs,length(unique_orientations),max_dim,max_index);

%preallocate matrix for each individual trial period
trial_data = zeros(numROIs,length(unique_orientations),max_dim,max_index);
%preallocate matrix for trial amplitude, which is trial_data - pre_data;
trial_amp = zeros(numROIs,length(unique_orientations),max_dim,max_index);
%preallocate matrix for trial amplitude, which is trial_data - pre_data;
trial_peak = zeros(numROIs,length(unique_orientations),max_dim,max_index);
%preallocate matrix for mean trial amplitude, which is trial_data - pre_data;
mean_amplitude = zeros(numROIs,length(unique_orientations),max_dim);
%preallocate matrix for mean trial amplitude, which is trial_data - pre_data;
std_amplitude = zeros(numROIs,length(unique_orientations),max_dim);
%preallocate matrix for median trial amplitude, which is trial_data - pre_data;
median_amplitude = zeros(numROIs,length(unique_orientations),max_dim);
%preallocate trial matrix 
trial_traces = zeros(numROIs,length(unique_orientations),max_dim,max_index,length(time_axis));
%preallocate matrix with trials that are included
included_trials = zeros(numROIs,length(unique_orientations),max_dim,max_index);
%preallocate reliability matrix 
reliability_index = zeros(numROIs,length(unique_orientations),max_dim);
%preallocate reliability matrix (gray screen)
pre_reliability_index = zeros(numROIs,length(unique_orientations),max_dim);
%Fano factor
fano_factor = nan(numROIs,length(unique_orientations),max_dim);
%Z-test across trials
z_test = nan(numROIs,length(unique_orientations),max_dim);
%T test across trial
t_test = nan(numROIs,length(unique_orientations),max_dim);
%pre means above 2 STD
large_pre_mean = zeros(numROIs,length(unique_orientations),max_dim);


%Loop through the ROIs
for l=1:numROIs
    n = 0;
    for ii = 1:length(unique_orientations)
        for jj = 1:1:sum(unique_rows==unique_orientations(ii)) 
            n = n+1;
            %create empty matrix to contain interpolated stim traces
            stim_traces=zeros(round(20*isi_dur + 10*stim_dur),length(indices{n}));
            for kk = 1:length(indices{n})
                %interpolate each trial on a 0.1 S scale  
                stim_traces(:,kk)=interp1(Datacell{ii,jj,kk}(:,1),Datacell{ii,jj,kk}(:,l+1),(round(-isi_dur*10):1:floor(stim_dur*10+isi_dur*10))/10);
                %take mean of individual trials from 0 to 1.5 seconds (due
                %to decay in stim
                trial_data(l,ii,jj,kk)=mean(Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=0&Datacell{ii,jj,kk}(:,1)<=1.5,l+1));
                %take peak of pre-stim periods 
                 pre_peak(l,ii,jj,kk)=max(Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=-0.5&Datacell{ii,jj,kk}(:,1)<0,l+1));
                %take mean of pre-stim periods 
                pre_mean(l,ii,jj,kk)=mean(Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=-1.5&Datacell{ii,jj,kk}(:,1)<0,l+1));
                %subtract trial_data - pre_data to get amplitude
                trial_amp(l,ii,jj,kk) = trial_data(l,ii,jj,kk) - pre_mean(l,ii,jj,kk);
                %find peak of individual traces
                trial_peak(l,ii,jj,kk)=max(Datacell{ii,jj,kk}(Datacell{ii,jj,kk}(:,1)>=0&Datacell{ii,jj,kk}(:,1)<=stim_dur,l+1));

            end
             %logical array for trials with pre_mean above 2 std
            large_pre_mean(l,ii,jj) = sum(permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) >= pre_peak_thresh);


            include_trials = permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) < pre_peak_thresh;
            
            mean_trace = nanmean(stim_traces(:,include_trials),2);
           
            %take mean of trial amplitudes
            median_amplitude(l,ii,jj) = median(trial_amp(l,ii,jj,include_trials));

            %take mean of trial amplitudes
            mean_amplitude(l,ii,jj) = mean(trial_amp(l,ii,jj,include_trials));
            
            %Take std of trial amplitudes
            std_amplitude(l,ii,jj) = std(trial_amp(l,ii,jj,include_trials));

            %take fano factor of included trials
            fano_factor(l,ii,jj) = (std_amplitude(l,ii,jj))/mean_amplitude(l,ii,jj);
            
            %Take reliability index of included trials 0 to 1.5 seconds
            corr_trials = corr(permute(trial_traces(l,ii,jj,include_trials,trial_time>=0 & trial_time<=1.5),[5,4,3,2,1])); %take correlation of interpolated trials
            idx = logical(triu(ones(size(corr_trials)),1)); %Make logical index of upper triangle of corr matrix 
            reliability_index(l,ii,jj) = (2/(sum(include_trials)^2 - sum(include_trials))) * sum(corr_trials(idx)); %take average of upper diagonal and divide by n to get mean   
           
            %calculate peak of mean trace 
            peak_data(l,ii,jj) = max(mean_trace(round(10*isi_dur)+1:round(10*isi_dur)+1+round(10*stim_dur)+10,:));
            
            %take the mean of the mean trace during stim period. Adding 10 to add 1
            %second to analyze part of the post stim period when still decaying.
            mean_data(l,ii,jj)= mean(mean_trace(round(10*isi_dur)+1:round(10*isi_dur)+1+round(10*stim_dur)+10,:));
           
            %take the mean of the mean trace during pre period
            mean_pre(l,ii,jj)= mean(mean_trace(round(10*isi_dur/2):round(10*isi_dur),:));
            
            %Take Z-test of included trials
            [~,z_test(l,ii,jj)] = ztest(trial_data(l,ii,jj,include_trials),mean(pre_mean(l,ii,jj,include_trials)),std(pre_mean(l,ii,jj,include_trials)),'Tail','right');
            %Take T-test of included trials
            [~,t_test(l,ii,jj)] = ttest(trial_data(l,ii,jj,include_trials),pre_mean(l,ii,jj,include_trials),'Tail','right');
            
            included_trials(l,ii,jj,:) = include_trials';
            
            %Take reliability index of included trials 0 to 1 seconds
            corr_trials = corr(permute(trial_traces(l,ii,jj,include_trials,time_axis>=0 & time_axis<=1),[5,4,3,2,1])); %take correlation of interpolated trials
            idx = logical(triu(ones(size(corr_trials)),1)); %Make logical index of upper triangle of corr matrix 
            reliability_index(l,ii,jj) = (2/(sum(include_trials)^2 - sum(include_trials))) * sum(corr_trials(idx)); %take average of upper diagonal and divide by n to get mean   
            
            %Take reliability index of pre data
            corr_trials = corr(permute(trial_traces(l,ii,jj,include_trials,time_axis>=-1.5 & time_axis<0),[5,4,3,2,1])); %take correlation of interpolated trials
            idx = logical(triu(ones(size(corr_trials)),1)); %Make logical index of upper triangle of corr matrix 
            pre_reliability_index(l,ii,jj) = (2/(sum(include_trials)^2 - sum(include_trials))) * sum(corr_trials(idx)); %take average of upper diagonal and divide by n to get mean   

            clear mean_trace stim_traces TF include_trials corr_trials idx
        end
    end
end

clear ii jj kk
%%
%Graph and calculate the data
if graph==1


%Loop through the ROIs
for l=1:numROIs
    responsive = 'unresponsive';
    index_isi = find(time_axis < 0 & time_axis >= -isi_dur/3);
    index_stim = find(time_axis >=0 & time_axis <= 1);
    %time_length = numel(time_axis);
    n=0; 
    %create new figure with ROI name and position. Adjust for your monitor
    %needs. 
    %figure('name',sprintf('Plot of ROI %d',l),'numbertitle','off','position',[414,391,472,257]);
    %loop through stims
    time_period = time_axis + isi_dur;
    figure('name',sprintf('Plot of ROI %d',l),'numbertitle','off','position',[414,391,472,257]);
    
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
            
            include_trials = permute(pre_peak(l,ii,jj,:),[4,5,3,2,1]) < pre_peak_thresh;

            %Take mean of traces 
            mean_trace=nanmean(stim_traces(:,include_trials),2);
                        
            %plot mean trace
            plot(time_period([index_isi,index_stim]),mean_trace([index_isi,index_stim]),'Color',c,'LineWidth',1);          
                       
            clear mean_trace TF stim_traces trial_traces
        end
        %label subplots with sig responses
        if min(t_test(l,ii,:))<=0.05 && max(mean_amplitude(l,ii,:))>=0.5 && sum(include_trials)>=5
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
for l = 1: numROIs
    roidata = ROIdata_Z(:,l)
    f = figure
    set(gca, 'FontSize', 20)
    start_ind = 1;
    
     plot(timestamps,roidata,'Color','k', 'LineWidth',3); hold on;
           

     for ii = 1:length(startevent)
        patch([startevent(ii),blankframes(ii),blankframes(ii),startevent(ii)], [-5, -5,12, 12],[0.5,0.5,0.5], 'FaceAlpha', 0.3,'LineStyle', 'none') 
        %patch([startevent(ii),blankframes(ii),blankframes(ii),startevent(ii)], [0, 0,200, 200],[0.5,0.5,0.5], 'FaceAlpha', 0.3,'LineStyle', 'none') 
    
     
       
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
%Report number of responsive cells in the command window
Responsive_mean=t_test<=0.05 & mean_amplitude>=0.5;
Responsive_cells = max(Responsive_mean,[],3);
Responsive_cells=sum(Responsive_cells,2);
Responsive_cells(Responsive_cells>=1)=1;
Responsive_cells=sum(Responsive_cells);
Report=sprintf('There were %s responsive cells',num2str(Responsive_cells));
disp(Report)
clear Report Responsive_cells
end
end
end
