%% concatenate all spines from each cell from each mouse and from each day into data matrix
%01-25-23 used this for code for my data!!
function analyzed_processed_data_function(path, savefile, mouse_files)
output_path1 = fullfile(path, 'Processed Data');
output_path2 = fullfile(path,'Coords Data');
savepath = fullfile(path, 'Analyzed Data');
savefile = fullfile(savepath,savefile);
all_stims = [];
stimfiles_name = {'binoc', 'contra','ipsi'};

%stimfiles_name = {'binoc'};
for s = 1:length(stimfiles_name)
    disp(['Doing ', stimfiles_name{s}])
    all_active_spine_trial = [];
    all_active_spine_trial_smooth = [];
    all_included_trials = [];
    all_included_trials_bAP = [];
    all_z_scored_trace = [];
    all_bAP_trial = [];
    all_roi_inds = [];
    all_x_coord = [];
    all_y_coord = [];
    all_trial_data = [];
    all_trial_amp = [];
    all_large_pre_mean = [];
    all_reliability_index = [];
    all_fano_factor = [];
    all_max_mean_amp = [];
    all_mean_amp = [];
    all_mean_amp_bAPs_removed = [];
    all_median_amp = [];
    %all_std_amp = [];
    all_ttest = [];
    all_ttest_bAPs_removed = [];
    all_artifact_trial = [];
    all_shaft_corr = [];

    all_day = [];
    all_mouse = [];
    all_cells = [];
    all_fovs = [];
    %all_spine = [];
   
    all_curve_fit_vonMises = [];
    all_Ori_pref_vonMises = [];
    all_OSI_vonMises = [];
    all_Dir_pref_vonMises = [];
    all_DSI_vonMises = [];
    all_GOF = [];
    all_coeffs_dir = [];
   

    all_Ori_pref_vector = [];
    all_OSI_vector = [];
    all_Dir_pref_vector = [];
    all_DSI_vector = [];

    all_Ori_pref_vector_bAPs_removed = [];
    all_OSI_vector_bAPs_removed = [];
    all_Dir_pref_vector_bAPs_removed = [];
    all_DSI_vector_bAPs_removed = [];

    all_peak = [];
    all_ISI_data = [];
  
%mouse_files = {'SOMA22'};
    for m = 1:length(mouse_files)
        disp(['Running mouse: ', mouse_files{m}])
       folders = dir(fullfile(output_path1,mouse_files{m}));
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
            Cells = dir(fullfile(output_path1,mouse_files{m},Day));
            Cells=Cells(~contains({Cells.name},'.'));
            for cc = 1:length(Cells)
                disp(['Running cell: ',Cells(cc).name])
                stims = dir(fullfile(output_path1,mouse_files{m},Day,Cells(cc).name));
                stims=stims(~contains({stims.name}, '.'));
                for iii = 1:length(stims)
                    disp(['Running segment: ',stims(iii).name])
                    sessions = dir(fullfile(output_path1,mouse_files{m},Day,Cells(cc).name, stims(iii).name));
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
                               if exist(fullfile(sessions(ss).folder, sessions(ss).name, 'mean_amp_ori_analysis_vonMises.mat'))
                                      load(fullfile(sessions(ss).folder, sessions(ss).name, 'mean_amp_ori_analysis_vonMises.mat'), 'OSI_vonMises', 'DSI_vonMises',...
                                      'Ori_pref_vonMises', 'Dir_pref_vonMises', 'Curve_Dir', 'coeffs_dir', 'GOF_vonMises_Dir');
                                       all_Ori_pref_vonMises = [all_Ori_pref_vonMises; Ori_pref_vonMises];
                                       all_OSI_vonMises = [all_OSI_vonMises; OSI_vonMises];
                                       all_Dir_pref_vonMises = [all_Dir_pref_vonMises; Dir_pref_vonMises];
                                       all_DSI_vonMises = [all_DSI_vonMises; DSI_vonMises];
                                      for rois = 1:size(t_test,1)
                                            all_curve_fit_vonMises = [all_curve_fit_vonMises; {Curve_Dir(rois,:)}];
                                            all_coeffs_dir = [all_coeffs_dir; {coeffs_dir(rois,:)}];
                                       end
                                       all_GOF = [all_GOF; GOF_vonMises_Dir];
                              else
                                    all_Ori_pref_vonMises = [all_Ori_pref_vonMises; NaN(size(max(mean_amplitude,[],2)))];
                                       all_OSI_vonMises = [all_OSI_vonMises; NaN(size(max(mean_amplitude,[],2)))];
                                       all_Dir_pref_vonMises = [all_Dir_pref_vonMises; NaN(size(max(mean_amplitude,[],2)))];
                                       all_DSI_vonMises = [all_DSI_vonMises;NaN(size(max(mean_amplitude,[],2)))];
                                       all_curve_fit_vonMises = [all_curve_fit_vonMises; repmat({NaN(1,360)},size(max(mean_amplitude,[],2)))];
                                       all_GOF = [all_GOF; NaN(size(max(mean_amplitude,[],2)))];
                                       all_coeffs_dir = [all_coeffs_dir; repmat({NaN(1,5)},size(max(mean_amplitude,[],2)))];
                              end 
                               all_artifact_trial = [all_artifact_trial;repmat({zeros(1,8)},size(max(mean_amplitude,[],2)))];
                               all_shaft_corr = [all_shaft_corr; zeros(size(max(mean_amplitude,[],2)))];
                               all_x_coord = [all_x_coord;NaN(size(max(mean_amplitude,[],2)))];
                               all_y_coord = [all_y_coord;NaN(size(max(mean_amplitude,[],2)))];

                                all_bAP_trial = [all_bAP_trial; repmat({NaN(1,8)},size(max(mean_amplitude,[],2)))];
                                all_included_trials_bAP = [all_included_trials_bAP;repmat({NaN(8,10)},size(max(mean_amplitude,[],2)))];
                                all_mean_amp_bAPs_removed = [all_mean_amp_bAPs_removed;repmat({NaN(1,8)},size(max(mean_amplitude,[],2)))];
                                all_ttest_bAPs_removed = [all_ttest_bAPs_removed; repmat({NaN(1,8)},size(max(mean_amplitude,[],2)))];  
                               
                               all_Ori_pref_vector_bAPs_removed = [all_Ori_pref_vector_bAPs_removed; zeros(size(max(mean_amplitude,[],2)))];
                               all_Dir_pref_vector_bAPs_removed = [all_Dir_pref_vector_bAPs_removed;  zeros(size(max(mean_amplitude,[],2)))];
                               all_OSI_vector_bAPs_removed = [all_OSI_vector_bAPs_removed;  zeros(size(max(mean_amplitude,[],2)))];
                               all_DSI_vector_bAPs_removed = [all_DSI_vector_bAPs_removed;  zeros(size(max(mean_amplitude,[],2)))];
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
                                  if exist(fullfile(sessions(ss).folder, sessions(ss).name, 'mean_amp_ori_analysis_bAPs_removed.mat.mat'))
                                        bap = load(fullfile(sessions(ss).folder, sessions(ss).name, 'mean_amp_ori_analysis_bAPs_removed.mat.mat'), 'OSI', 'DSI', 'Ori_pref', 'Dir_pref');
                                        all_Ori_pref_vector_bAPs_removed = [all_Ori_pref_vector_bAPs_removed; bap.Ori_pref];
                                        all_Dir_pref_vector_bAPs_removed = [all_Dir_pref_vector_bAPs_removed; bap.Dir_pref];
                                        all_OSI_vector_bAPs_removed = [all_OSI_vector_bAPs_removed; bap.OSI];
                                        all_DSI_vector_bAPs_removed = [all_DSI_vector_bAPs_removed; bap.DSI];
                                  else
                                       all_Ori_pref_vector_bAPs_removed = [all_Ori_pref_vector_bAPs_removed; zeros(size(max(mean_amplitude,[],2)))];
                                       all_Dir_pref_vector_bAPs_removed = [all_Dir_pref_vector_bAPs_removed;  zeros(size(max(mean_amplitude,[],2)))];
                                       all_OSI_vector_bAPs_removed = [all_OSI_vector_bAPs_removed;  zeros(size(max(mean_amplitude,[],2)))];
                                       all_DSI_vector_bAPs_removed = [all_DSI_vector_bAPs_removed;  zeros(size(max(mean_amplitude,[],2)))];
                                       all_bAP_trial = [all_bAP_trial; repmat({NaN(1,8)},size(max(mean_amplitude,[],2)))];
                                       all_included_trials_bAP = [all_included_trials_bAP;repmat({NaN(8,10)},size(max(mean_amplitude,[],2)))];
                                       all_mean_amp_bAPs_removed = [all_mean_amp_bAPs_removed;repmat({NaN(1,8)},size(max(mean_amplitude,[],2)))];
                                       all_ttest_bAPs_removed = [all_ttest_bAPs_removed; repmat({NaN(1,8)},size(max(mean_amplitude,[],2)))];  

                                  end
                                  if exist(fullfile(sessions(ss).folder, sessions(ss).name, 'mean_amp_ori_analysis_vonMises.mat'))
                                      load(fullfile(sessions(ss).folder, sessions(ss).name, 'mean_amp_ori_analysis_vonMises.mat'), 'OSI_vonMises', 'DSI_vonMises',...
                                      'Ori_pref_vonMises', 'Dir_pref_vonMises', 'Curve_Dir', 'coeffs_dir', 'GOF_vonMises_Dir');
                                       all_Ori_pref_vonMises = [all_Ori_pref_vonMises; Ori_pref_vonMises];
                                       all_OSI_vonMises = [all_OSI_vonMises; OSI_vonMises];
                                       all_Dir_pref_vonMises = [all_Dir_pref_vonMises; Dir_pref_vonMises];
                                       all_DSI_vonMises = [all_DSI_vonMises; DSI_vonMises];
                                      
                                       all_GOF = [all_GOF; GOF_vonMises_Dir];
                                       for rois = 1:size(t_test,1)
                                            all_curve_fit_vonMises = [all_curve_fit_vonMises; {Curve_Dir(rois,:)}];
                                            all_coeffs_dir = [all_coeffs_dir; {coeffs_dir(rois,:)}];
                                       end
                                  else
                                       all_Ori_pref_vonMises = [all_Ori_pref_vonMises; NaN(size(max(mean_amplitude,[],2)))];
                                       all_OSI_vonMises = [all_OSI_vonMises; NaN(size(max(mean_amplitude,[],2)))];
                                       all_Dir_pref_vonMises = [all_Dir_pref_vonMises; NaN(size(max(mean_amplitude,[],2)))];
                                       all_DSI_vonMises = [all_DSI_vonMises;NaN(size(max(mean_amplitude,[],2)))];
                                       all_curve_fit_vonMises = [all_curve_fit_vonMises; repmat({NaN(1,360)},size(max(mean_amplitude,[],2)))];
                                       all_GOF = [all_GOF; NaN(size(max(mean_amplitude,[],2)))];
                                       all_coeffs_dir = [all_coeffs_dir; repmat({NaN(1,5)},size(max(mean_amplitude,[],2)))];                                  end
                                  
                        
                                  if exist(fullfile(output_path2,mouse_files{m},Day,Cells(cc).name,[stims(iii).name,'.mat'],"file"))
                                        load(fullfile(output_path2,mouse_files{m},Day,Cells(cc).name,[stims(iii).name,'.mat']));
                                        all_x_coord = [all_x_coord; x_coords'];
                                        all_y_coord = [all_y_coord;y_coords'];
                                  else
                                        all_x_coord = [all_x_coord;NaN(size(max(mean_amplitude,[],2)))];
                                        all_y_coord = [all_y_coord;NaN(size(max(mean_amplitude,[],2)))];
    
                                  end
                                   total_artifacts = sum(artifact_trial,1);
                                   %total_artifacts(total_artifacts>1)=1;
                                   all_shaft_corr = [all_shaft_corr; shaft_corr];
                                   

                                   for rois = 1:size(t_test,1)
                                         all_artifact_trial = [all_artifact_trial; {total_artifacts}];
                                         all_active_spine_trial = [ all_active_spine_trial; {active_spine_trial(rois,:,:,:)}];
                                         all_active_spine_trial_smooth = [ all_active_spine_trial_smooth; {active_spine_trial_smooth(rois,:,:,:)}];
                                        
                                         all_mean_amp_bAPs_removed = [all_mean_amp_bAPs_removed; {mean_amplitude_bAPs_removed(rois,:)}];
                                         all_ttest_bAPs_removed = [all_ttest_bAPs_removed; {t_test_bAPs_removed(rois,:)}];  
                                         all_bAP_trial = [all_bAP_trial; {bAP_trial}];
                                         all_included_trials_bAP = [all_included_trials_bAP;{squeeze(included_trials_bAP(rois, :, :, :))}];                                     
                                  

                                   end
                              else
                                   all_artifact_trial = [all_artifact_trial; NaN];
                                   all_bAP_trial = [all_bAP_trial; NaN];
                                   all_included_trials_bAP = [all_included_trials_bAP;NaN];
                                   all_included_trials = [all_included_trials; NaN];

                                   all_shaft_corr = [all_shaft_corr; NaN];

                                   all_x_coord = [all_x_coord;NaN];
                                   all_y_coord = [all_y_coord;NaN];
                                   
                                   all_Ori_pref_vonMises = [all_Ori_pref_vonMises; NaN];
                                   all_OSI_vonMises = [all_OSI_vonMises; NaN];
                                   all_Dir_pref_vonMises = [all_Dir_pref_vonMises; NaN];
                                   all_DSI_vonMises = [all_DSI_vonMises; NaN];
                                   all_curve_fit_vonMises = [all_curve_fit_vonMises; NaN];
                                   all_GOF = [all_GOF; NaN];
                                   all_coeffs_dir = [all_coeffs_dir; NaN];

                                   all_Ori_pref_vector = [all_Ori_pref_vector; NaN];
                                   all_OSI_vector = [all_OSI_vector; NaN];
                                   all_Dir_pref_vector = [all_Dir_pref_vector; NaN];
                                   all_DSI_vector = [all_DSI_vector; NaN];

                                   all_Ori_pref_vector_bAPs_removed = [all_Ori_pref_vector_bAPs_removed; NaN];
                                   all_Dir_pref_vector_bAPs_removed = [all_Dir_pref_vector_bAPs_removed;  NaN];
                                   all_OSI_vector_bAPs_removed = [all_OSI_vector_bAPs_removed;  NaN];
                                   all_DSI_vector_bAPs_removed = [all_DSI_vector_bAPs_removed;  NaN];

                                  
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
                                   all_mean_amp_bAPs_removed = [all_mean_amp_bAPs_removed; NaN];
                                   all_ttest_bAPs_removed = [all_ttest_bAPs_removed; NaN];
                                   all_median_amp = [all_median_amp;NaN];
                                   all_ttest = [all_ttest;NaN];
                                   all_max_mean_amp = [all_max_mean_amp;NaN];
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
                      
                            
                            all_median_amp = [all_median_amp;{median_amplitude(rois,:)}];
                            all_large_pre_mean = [all_large_pre_mean; {large_pre_mean(rois,:)}];
                           
                            
                            all_reliability_index = [all_reliability_index;{reliability_index(rois,:)}];
                            all_fano_factor = [all_fano_factor;{fano_factor(rois,:)}];
                            all_included_trials = [all_included_trials; {included_trials(rois,:,:,:)}];
                            all_z_scored_trace = [all_z_scored_trace; {ROIdata_Z(:, rois)}];

                        end
                         

                       all_max_mean_amp = [all_max_mean_amp;max(mean_amplitude,[],2)];
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
    stim.all_bAP_trial = all_bAP_trial;
    stim.all_included_trial = all_included_trials;
    stim.all_included_trials_bAP = all_included_trials_bAP;
    stim.all_shaft_corr = all_shaft_corr;
    stim.all_large_pre_mean = all_large_pre_mean;

    stim.all_ttest = all_ttest;
    stim.all_ttest_bAPs_removed = all_ttest_bAPs_removed;
    stim.all_fano_factor = all_fano_factor;
    stim.all_trial_data = all_trial_data;
    stim.all_median_amp = all_median_amp;
    stim.all_mean_amp = all_mean_amp;
    stim.all_mean_amp_bAPs_removed = all_mean_amp_bAPs_removed;
    stim.all_max_reliability_index_day = all_reliability_index;
    stim.all_OSI_vector = all_OSI_vector;
    stim.all_DSI_vector = all_DSI_vector;
    stim.all_Dir_pref_vector = all_Dir_pref_vector;
    stim.all_Ori_pref_vector = all_Ori_pref_vector;
    stim.all_OSI_vector_bAPs_removed = all_OSI_vector_bAPs_removed;
    stim.all_DSI_vector_bAPs_removed = all_DSI_vector_bAPs_removed;
    stim.all_Dir_pref_vector_bAPs_removed = all_Dir_pref_vector_bAPs_removed;
    stim.all_Ori_pref_vector_bAPs_removed = all_Ori_pref_vector_bAPs_removed;
    stim.all_trial_amp = all_trial_amp;
    stim.all_OSI_vonMises = all_OSI_vonMises;
    stim.all_DSI_vonMises = all_DSI_vonMises;
    stim.all_Dir_vonMises = all_Dir_pref_vonMises;
    stim.all_Ori_vonMises = all_Ori_pref_vonMises;
    stim.all_curve_fit = all_curve_fit_vonMises;
    stim.all_GOF = all_GOF;
    stim.all_coeffs_dir = all_coeffs_dir;
    stim.all_mean_Z_day = all_max_mean_amp;
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

save(savefile, 'all_stims');
end

