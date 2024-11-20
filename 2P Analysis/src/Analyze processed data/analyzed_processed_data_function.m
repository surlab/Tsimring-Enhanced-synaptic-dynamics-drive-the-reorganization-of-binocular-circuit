%% concatenate all spines from each cell from each mouse and from each day into data matrix
%01-25-23 used this for code for my data!!
function analyzed_processed_data_function(inputpath, savepath,filename, mouse_files)
all_stims = [];
stimfiles_name = {'binoc', 'contra','ipsi'};

%stimfiles_name = {'binoc'};
for s = 1:length(stimfiles_name)
    disp(['Doing ', stimfiles_name{s}])
    all_active_spine_trial = [];
    all_active_spine_trial_smooth = [];
    all_included_trials = [];

    all_z_scored_trace = [];
    all_roi_inds = [];
    all_trial_data = [];
    all_trial_amp = [];
    all_large_pre_mean = [];
    all_reliability_index = [];
    all_fano_factor = [];
    all_mean_amp = [];

    all_ttest = [];
    all_artifact_trial = [];
    all_shaft_corr = [];

    all_day = [];
    all_mouse = [];
    all_cells = [];
    all_fovs = [];
    %all_spine = [];

    all_Ori_pref_vector = [];
    all_OSI_vector = [];
    all_Dir_pref_vector = [];
    all_DSI_vector = [];

    all_peak = [];
    all_ISI_data = [];
  
%mouse_files = {'SOMA22'};
    for m = 1:length(mouse_files)
        disp(['Running mouse: ', mouse_files{m}])
       folders = dir(fullfile(inputpath,mouse_files{m}));
       Days = {folders.name};
       Days=Days(~contains({folders.name}, '.'));
       for i = 1:length(Days)
           
            Day = Days{i};
            disp(['Running day: ',Day])
            d = str2double(Day(2:end));
            if d >= 50
                index = 1:2:16; %get 8 directions not 16 
            else
                index = 1:8;
            end
            Cells = dir(fullfile(inputpath,mouse_files{m},Day));
            Cells=Cells(~contains({Cells.name},'.'));
            for cc = 1:length(Cells)
                disp(['Running cell: ',Cells(cc).name])
                stims = dir(fullfile(inputpath,mouse_files{m},Day,Cells(cc).name));
                stims=stims(~contains({stims.name}, '.'));
                for iii = 1:length(stims)
                    disp(['Running segment: ',stims(iii).name])
                    sessions = dir(fullfile(inputpath,mouse_files{m},Day,Cells(cc).name, stims(iii).name));
                    sessions=sessions(~contains({sessions.name}, '.'));
                    if length(sessions)==3
                    sessions=sessions(contains({sessions.name}, stimfiles_name{s},'IgnoreCase',true ));
                   
                    for ss = 1:length(sessions)
                          if contains(sessions(ss).name, 'Cell') || contains(sessions(ss).name, 'Soma')
                              load(fullfile(sessions(ss).folder, sessions(ss).name, 'mean_amp_ori_analysis.mat'), 'OSI', 'DSI', 'Ori_pref', 'Dir_pref')
                              load(fullfile(sessions(ss).folder, sessions(ss).name, 'normalized_data_by_stim.mat'), ...
                              't_test', 'median_amplitude', 'mean_amplitude','timestamps', ...
                              'fano_factor','startevent','reliability_index', 'included_trials',...
                              'large_pre_mean','blankframes','ROIdata_Z', 'peak_data', 'trial_amp');
                              
                              if(size(mean_amplitude,3) > 1)
                                    ind = find_preferred_SF(mean_amplitude,t_test);
                                    mean_amplitude = mean_amplitude(:,:,ind);
                                    median_amplitude = median_amplitude(:,:,ind);
                                    t_test = t_test(:,:,ind);
                                    reliability_index = reliability_index(:,:,ind);
                                    fano_factor = fano_factor(:,:,ind);
                                    large_pre_mean =  large_pre_mean(:,:,ind);
                                    peak_data = peak_data(:,ind);
                                    included_trials = included_trials(:,:,ind,:);
                                    trial_amp = trial_amp(:,:,ind,:);
                              end
                               all_artifact_trial = [all_artifact_trial;repmat({zeros(1,8)},size(max(mean_amplitude,[],2)))];
                               all_shaft_corr = [all_shaft_corr; zeros(size(max(mean_amplitude,[],2)))];

                               all_active_spine_trial = [ all_active_spine_trial; repmat({NaN(1,8,1,10)},size(max(mean_amplitude,[],2)))];
                               all_active_spine_trial_smooth = [ all_active_spine_trial_smooth; repmat({NaN(1,8,1,10)},size(max(mean_amplitude,[],2)))];

                                

                          else
                              if exist(fullfile(sessions(ss).folder, sessions(ss).name,'normalized_data_by_stim.mat'), "file")
                                  load(fullfile(sessions(ss).folder, sessions(ss).name, 'mean_amp_ori_analysis.mat'), 'OSI', 'DSI', 'Ori_pref', 'Dir_pref')
                                  load(fullfile(sessions(ss).folder, sessions(ss).name,'normalized_data_by_stim.mat'), ...
                                  't_test','t_test_bAPs_removed','mean_amplitude_bAPs_removed', ...
                                  'mean_amplitude','artifact_trial','shaft_corr', ...
                                  'fano_factor','startevent','timestamps','reliability_index', ...
                                  'large_pre_mean','blankframes','ROIdata_Z','active_spine_trial_smooth','active_spine_trial', ...
                                  'peak_data','bAP_trial', 'included_trials', 'included_trials_bAP', 'trial_amp', 'median_amplitude');
                                  
                                 
                                   total_artifacts = sum(artifact_trial,1);
                                   %total_artifacts(total_artifacts>1)=1;
                                   all_shaft_corr = [all_shaft_corr; shaft_corr];
                                   

                                   for rois = 1:size(t_test,1)
                                         all_artifact_trial = [all_artifact_trial; {total_artifacts}];
                                         all_active_spine_trial = [ all_active_spine_trial; {active_spine_trial(rois,:,:,:)}];
                                         all_active_spine_trial_smooth = [ all_active_spine_trial_smooth; {active_spine_trial_smooth(rois,:,:,:)}];
                                   end
                              else
                                   all_artifact_trial = [all_artifact_trial; NaN];
                                   all_included_trials = [all_included_trials; NaN];

                                   all_shaft_corr = [all_shaft_corr; NaN];

                                
                                   all_Ori_pref_vector = [all_Ori_pref_vector; NaN];
                                   all_OSI_vector = [all_OSI_vector; NaN];
                                   all_Dir_pref_vector = [all_Dir_pref_vector; NaN];
                                   all_DSI_vector = [all_DSI_vector; NaN];
                                  
                                   all_peak = [all_peak; NaN];
                                   all_roi_inds = [all_roi_inds;NaN];
                                   all_cells = [all_cells; repmat({Cells(cc).name},1)];
                                   all_fovs = [all_fovs; repmat({stims(iii).name},1)];
                                   all_day = [all_day; repmat({Day}, 1)];
                                   all_active_spine_trial = [ all_active_spine_trial; NaN];
                                   all_active_spine_trial_smooth = [ all_active_spine_trial_smooth; NaN];

                                   all_mouse = [all_mouse; repmat({mouse_files{m}}, 1)];
                                   all_reliability_index = [all_reliability_index;NaN];
                        
                                   all_fano_factor = [all_fano_factor; NaN];
                                   all_large_pre_mean = [all_large_pre_mean; NaN];
                                   all_trial_data = [all_trial_data;NaN];
                                   all_ISI_data = [all_ISI_data;NaN];
                                   all_mean_amp = [all_mean_amp;NaN];
                                   all_ttest = [all_ttest;NaN];
                                   
                                   all_trial_amp = [all_trial_amp;NaN];
                                   all_z_scored_trace = [all_z_scored_trace; NaN];

                                   break;
                              end

                        end
                        
                           
                        all_cells = [all_cells; repmat({Cells(cc).name},size(max(mean_amplitude,[],2)))];
                        all_fovs = [all_fovs; repmat({stims(iii).name},size(max(mean_amplitude,[],2)))];
                        all_day = [all_day; repmat({Day}, size(max(mean_amplitude,[],2)))];
                        all_mouse = [all_mouse; repmat({mouse_files{m}}, size(max(mean_amplitude,[],2)))];
                        
                       
                        %all_large_pre_mean = [all_large_pre_mean; max(large_pre_mean,[],2)];
                        
                        %extract flourescence during stim (for spine-spine
                        %correlations) and mean amp for each spine
                        
                        for rois = 1:size(t_test,1)
                            trial_data = [];
                            isi_data = [];
                            for time = 1:length(startevent)
                                if time < length(startevent)
                                    isi_data = [isi_data; ...
                                        ROIdata_Z(timestamps>=blankframes(time)&...
                                        timestamps<startevent(time+1),rois)];
                                end
                                trial_data = [trial_data;...
                                    ROIdata_Z(timestamps>=startevent(time)&...
                                    timestamps<blankframes(time),rois)];
                            end
                            trial_temp = trial_amp(rois,:,:,:);
                            all_trial_amp = [all_trial_amp; {trial_temp}];
                            all_ISI_data = [all_ISI_data; {isi_data}];
                            all_trial_data = [all_trial_data;{trial_data}];
                            all_ttest = [all_ttest;{t_test(rois,:)}];
                            all_mean_amp = [all_mean_amp;{mean_amplitude(rois,:)}];
                      
                            
                            all_large_pre_mean = [all_large_pre_mean; {large_pre_mean(rois,:)}];
                           
                            
                            all_reliability_index = [all_reliability_index;{reliability_index(rois,:)}];
                            all_fano_factor = [all_fano_factor;{fano_factor(rois,:)}];
                            all_included_trials = [all_included_trials; {included_trials(rois,:,:,:)}];
                            all_z_scored_trace = [all_z_scored_trace; {ROIdata_Z(:, rois)}];

                        end
                         

                       
                       all_Ori_pref_vector = [all_Ori_pref_vector; Ori_pref];
                       all_OSI_vector = [all_OSI_vector; OSI];
                       all_Dir_pref_vector = [all_Dir_pref_vector; Dir_pref];
                       all_DSI_vector = [all_DSI_vector; DSI];
                       all_peak = [all_peak; max(peak_data,[],2)];
                       all_roi_inds = [all_roi_inds;[1:size(t_test,1)]'];
                       
                      

                        

                    end
                    end
                end
            end
       end
    end
    stim.all_active_spine_trials = all_active_spine_trial;
    stim.all_active_spine_trials_smooth = all_active_spine_trial_smooth;
    stim.all_z_scored_trace = all_z_scored_trace;
    stim.all_ISI_data = all_ISI_data;
    stim.all_roi_inds = all_roi_inds;
    stim.all_artifact_trial = all_artifact_trial;
    stim.all_included_trial = all_included_trials;
    stim.all_shaft_corr = all_shaft_corr;
    stim.all_large_pre_mean = all_large_pre_mean;

    stim.all_ttest = all_ttest;
    stim.all_fano_factor = all_fano_factor;
    stim.all_trial_data = all_trial_data;
    stim.all_mean_amp = all_mean_amp;
    stim.all_OSI_vector = all_OSI_vector;
    stim.all_DSI_vector = all_DSI_vector;
    stim.all_Dir_pref_vector = all_Dir_pref_vector;
    stim.all_Ori_pref_vector = all_Ori_pref_vector;
    stim.all_trial_amp = all_trial_amp;
    stim.all_cells = all_cells;
    stim.all_fovs = all_fovs;
    stim.all_mouse = all_mouse;
    stim.all_day = all_day;
    stim.all_peak = all_peak;
    stim.name = stimfiles_name{s};
    all_stims = [all_stims, stim];
end

if ~isdir(savepath)
    mkdir(savepath)
end

save(fullfile(savepath, filename),'all_stims');
end

