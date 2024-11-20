%% Creates a table for every spine that was tracked between two timepoints
function create_table_tracked_spines_function(distance_path, path, spinefile, savefile)
chronic_path = fullfile(distance_path, 'FOV_alignment');

processed_path = fullfile(path, 'Analyzed Data');
load(fullfile(processed_path, spinefile), 'all_stim_table');

%%
extract_columns = {'all_roi_inds','nonresp','resp','soma_resp',...
                        'all_Dir_pref_vector','all_DSI_vector','all_Ori_pref_vector','all_OSI_vector',...
                        'all_Dir_vonMises','all_DSI_vonMises','all_OSI_vonMises',...
                        'soma_dir_pref_vector', 'soma_DSI_vector', 'soma_ori_pref_vector','soma_OSI_vector', ...
                        'soma_dir_pref_vonMises', 'soma_OSI_vonMises', 'soma_DSI_vonMises',...
                        'all_ori_alignment_vector','all_dir_corr','all_max_reliability_index_day', 'all_fano_factor', ...
                        'co_activity_trial', 'neighbor_spine_distance', 'neighbor_spine_resp', 'all_shaft_corr', 'mean_co_activity_trial',...
                        'co_activity_ISI', 'mean_co_activity_ISI', 'mean_co_activity_mean_amp', 'mean_co_activity_trial_amp','spine_area',...
                         'all_trial_data', 'all_z_scored_trace', 'all_active_spine_trials_smooth', ...
                        'spine_mean_resp', 'soma_mean_resp','co_activity_mean_amp', 'co_activity_trial_amp', 'mean_ori_offset_spines', 'ori_offset_spines',...
                        'all_trial_amp', 'all_included_trial', 'all_day'};

mice = dir(chronic_path);
mice = mice(~ismember({mice.name},{'.','..', '.DS_Store', 'ReadMe.txt', 'annotated_images'}));
D1_D5_table = table;
D5_D10_table = table;

sessions = {'contra', 'ipsi', 'binoc'};
for m = 1:length(mice)
    if ismember(mice(m).name, 'BM016')
        g = 1;
    end
    is_blinded = 0;
    disp(mice(m).name);
    if contains(mice(m).name,'BM016')
        k = 1;
    end
    cells = dir(fullfile(mice(m).folder, mice(m).name));
    cells = cells(~ismember({cells.name},{'.','..', '.DS_Store', 'ReadMe.txt'}));
    T = table;
    blind_folder = fullfile(mice(m).folder, mice(m).name, 'key_to_day_folders.csv');
    if exist(blind_folder,"file")
        is_blinded = 1;
        tbl = readtable(blind_folder);
       
    end
    for c = 1:length(cells)
        disp(cells(c).name)
        comp_days = dir(fullfile(cells(c).folder, cells(c).name));
        comp_days = comp_days(contains({comp_days.name},'Compare'));
       
        [~,inds] = sort({comp_days.name});
        comp_days = comp_days(inds);
        mouse_cell = [mice(m).name, '-', cells(c).name];
        for d = 1:length(comp_days)
            comp_day = comp_days(d).name;
            disp(comp_day)
            comp_day = split(comp_day, '-');
            
            d1 = comp_day{1};
            d2 = comp_day{2};
            if is_blinded
                 tbl_filt = tbl(ismember(tbl.Var1_3, cells(c).name),{'Var1_5';'Var1_4'});
                 real_date = tbl_filt.Var1_4;
                 blind_date = tbl_filt.Var1_5;
                 d1 = erase(d1, 'Compare: ');
                 ind1 = find(contains(blind_date, d1));
                 ind2 = find(contains(blind_date, d2)); 
                 d1 = real_date(ind1);
                 d2 = real_date(ind2);
            end
            T1 = which_day(d1);
            
            T2 = which_day(d2);
            dends = dir(fullfile(comp_days(d).folder,comp_days(d).name));
            dends = dends(~ismember({dends.name},{'.','..', '.DS_Store', 'ReadMe.txt'}));
            for dd = 1:length(dends)
                disp(dends(dd).name)
                filter_table_T1 = all_stim_table(contains...
                (all_stim_table.mouse_cell, mouse_cell)...
                &ismember(all_stim_table.days, T1) & ...
                ~isnan(all_stim_table.all_peak) &contains...
                (all_stim_table.all_fovs,dends(dd).name),:);
                
                filter_table_T2 = all_stim_table(contains...
                (all_stim_table.mouse_cell, mouse_cell)...
                &ismember(all_stim_table.days, T2) & ...
                ~isnan(all_stim_table.all_peak) &contains...
                (all_stim_table.all_fovs,dends(dd).name),:);
                if ~isempty(filter_table_T1) & ~isempty(filter_table_T2)
                load(fullfile(dends(dd).folder, dends(dd).name, 'alignment_based_on_fiducial.mat'),'modify_mat');
                filter_modify_mat = modify_mat;
                f1 = modify_mat(:,1)==0;
                f2 = modify_mat(:,2)==0;
                filter_modify_mat = filter_modify_mat(~(f1&f2),:);
                f1 = filter_modify_mat(:,1)>0;
                f2 = filter_modify_mat(:,2)>0;
                f3 = filter_modify_mat(:,2)==0;
                f4 = filter_modify_mat(:,1)==0;
                [uniqueA i j] = unique(filter_modify_mat(f1,1),'first');
                indexToDupes1 = find(not(ismember(1:numel(filter_modify_mat(f1,1)),i)))
                [uniqueA i j] = unique(filter_modify_mat(f2,2),'first');
                indexToDupes2 = find(not(ismember(1:numel(filter_modify_mat(f2,2)),i)))
                
                if ~isempty(indexToDupes1) || ~isempty(indexToDupes2) 
                    k = 1;
                end

                diff_f1 = diff(sort(filter_modify_mat(f2,2)));
                diff_f2 = diff(sort(filter_modify_mat(f1,1)));
                if ~isempty(find(diff_f1>1)) || ~isempty(find(diff_f2>1)) 
                    k = 1;
                end
                for s = 1:length(sessions)
                    
                    filter_session1 = filter_table_T1.session;
                    filter_session2 = filter_table_T2.session;
                    session_table_T1 = filter_table_T1(...
                        contains([filter_session1{:}], sessions(s)),:);
                    session_table_T2 = filter_table_T2(...
                        contains([filter_session2{:}], sessions(s)),:);
                    if ~isempty(session_table_T1) & ~isempty(session_table_T2)
                    eliminated_ind = find(f1&f3);
                    formed_ind = find(f4&f2);
                    retained_ind = find(f1&f2);

                    retained_spines1 = session_table_T1(filter_modify_mat(retained_ind,1),extract_columns);
                    retained_spines2 = session_table_T2(filter_modify_mat(retained_ind,2),extract_columns);
                    eliminated_spines = session_table_T1(filter_modify_mat(eliminated_ind,1),extract_columns);
                    formed_spines = session_table_T2(filter_modify_mat(formed_ind,2),extract_columns);
                    
                    %extract retained spines
                    retained_spines = table(retained_spines1,retained_spines2);
                    retained_spines = renamevars(retained_spines,["retained_spines1","retained_spines2"],[string(T1), string(T2)]);
                    %retained_spines = splitvars(retained_spines);
                    retained = repmat("retained", height(retained_spines),1);
  
                    A1 = array2table(zeros(height(eliminated_spines),length(extract_columns)), 'VariableNames',extract_columns);
                    A1.soma_mean_resp = repmat({nan}, height(eliminated_spines),1);
                    A1.spine_mean_resp = repmat({nan}, height(eliminated_spines),1);
                    A1.all_max_reliability_index_day = repmat({nan}, height(eliminated_spines),1);
                    A1.all_fano_factor = repmat({nan}, height(eliminated_spines),1);
                    A1.all_trial_data =  repmat({nan}, height(eliminated_spines),1);
                    A1.all_z_scored_trace = repmat({nan}, height(eliminated_spines),1);
                    A1.all_active_spine_trials_smooth = repmat({nan}, height(eliminated_spines),1);
                    A1.all_trial_amp = repmat({nan}, height(eliminated_spines),1);
                    A1.all_included_trial = repmat({nan}, height(eliminated_spines),1);
                    A1.all_day = repmat({nan}, height(eliminated_spines),1);
                    
                    A2 = array2table(zeros(height(formed_spines),length(extract_columns)), 'VariableNames',extract_columns);
                    A2.soma_mean_resp = repmat({nan}, height(formed_spines),1);
                    A2.spine_mean_resp = repmat({nan}, height(formed_spines),1);
                    A2.all_max_reliability_index_day = repmat({nan}, height(formed_spines),1);
                    A2.all_fano_factor = repmat({nan}, height(formed_spines),1);
                    A2.all_trial_data =  repmat({nan}, height(formed_spines),1);
                    A2.all_z_scored_trace =  repmat({nan}, height(formed_spines),1);
                    A2.all_active_spine_trials_smooth = repmat({nan}, height(formed_spines),1);
                    A2.all_trial_amp = repmat({nan}, height(formed_spines),1);
                    A2.all_included_trial = repmat({nan}, height(formed_spines),1);
                    A2.all_day = repmat({nan}, height(formed_spines),1);

                    %extract eliminated spines
                    eliminated_spines = table(eliminated_spines, A1);
                    eliminated_spines = renamevars(eliminated_spines,["eliminated_spines","A1"],[string(T1), string(T2)]);
                    %eliminated_spines = splitvars(eliminated_spines);
                    lost = repmat("lost", height(eliminated_spines),1);
                    
                    %extract formed spines
                    formed_spines = table(A2, formed_spines);
                    formed_spines = renamevars(formed_spines,["A2","formed_spines"],[string(T1), string(T2)]);
                    %formed_spines = splitvars(formed_spines);
                    formed = repmat("formed", height(formed_spines),1);
                    
                    temp_table = [retained_spines; eliminated_spines;formed_spines];
                    temp_table.all_mice_cells = repmat(string(mouse_cell), height(temp_table),1);
                    temp_table.all_fovs = repmat(string(dends(dd).name), height(temp_table),1);
                    temp_table.session = repmat(sessions(s), height(temp_table),1);
                    temp_table.structure_type = [retained; lost; formed];
                   
                    temp_table.dendritic_type = repmat(filter_table_T1(1,:).dend_type, height(temp_table),1);
                    temp_table.dendritic_depth = repmat(filter_table_T1(1,:).dend_depth, height(temp_table),1);
                    temp_table.soma_depth = repmat(filter_table_T1(1,:).soma_depth, height(temp_table),1);

                    if ismember(T1, 'D1')
                        D1_D5_table = [D1_D5_table; temp_table];

                    elseif ismember(T1, 'D5')
                        D5_D10_table = [D5_D10_table; temp_table];
                       
                    end
                    end
                end
                end
            end
                


        end
    end
end
save(fullfile(distance_path, savefile), 'D1_D5_table', 'D5_D10_table');
%%
function day = which_day(day)
    if contains(day, 'p22') || contains(day, 'p23') || contains(day, 'p24')
        day = 'D1';
    elseif contains(day, 'p27') ||contains(day, 'p28')|| contains(day, 'p29') 
        day = 'D5';
    else
        day = 'D10';
    end
end
end
