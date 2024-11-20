%% Create table for all spines tuning, all sessions, all days, alignment to soma
% update 9/28/2023 to account for resizing distances, and corrected soma to
% update 02/03/2024 to account for spine size
% spine alignment
function create_table_of_all3sessions_function(savepath, input_path, all_stims, all_vis_stims, savefile)
%% separate 3 sessions into table format
stim_binoc = all_stims(strcmp({all_stims.name},'binoc')==1);
stim_contra = all_stims(strcmp({all_stims.name}, 'contra')==1);
stim_ipsi = all_stims(strcmp({all_stims.name}, 'ipsi')==1);
stim_binoc = rmfield(stim_binoc,'name');
stim_ipsi = rmfield(stim_ipsi,'name');
stim_contra = rmfield(stim_contra, 'name');

vis_stim_binoc = all_vis_stims(strcmp({all_vis_stims.session},'binoc')==1);
vis_stim_contra = all_vis_stims(strcmp({all_vis_stims.session}, 'contra')==1);
vis_stim_ipsi = all_vis_stims(strcmp({all_vis_stims.session}, 'ipsi')==1);
vis_stim_binoc = rmfield(vis_stim_binoc,'session');
vis_stim_ipsi = rmfield(vis_stim_ipsi,'session');
vis_stim_contra = rmfield(vis_stim_contra, 'session');

stim_table_binoc = struct2table(stim_binoc);
stim_table_contra = struct2table(stim_contra);
stim_table_ipsi = struct2table(stim_ipsi);

vis_stim_table_binoc = struct2table(vis_stim_binoc);
vis_stim_table_contra = struct2table(vis_stim_contra);
vis_stim_table_ipsi = struct2table(vis_stim_ipsi);
%% Combine vis responsive info to table
combine_stim_binoc = [vis_stim_table_binoc,stim_table_binoc];
combine_stim_contra = [vis_stim_table_contra, stim_table_contra];
combine_stim_ipsi = [vis_stim_table_ipsi, stim_table_ipsi];


%% separate dendrites
%binoc
stim_table_dend_binoc = combine_stim_binoc(~contains(combine_stim_binoc.all_fovs,'Soma_0'),:);
stim_table_dend_binoc.mouse_cell = strcat(stim_table_dend_binoc.all_mouse, '-', stim_table_dend_binoc.all_cells);
stim_table_dend_binoc = removevars(stim_table_dend_binoc, 'all_mouse');
stim_table_dend_binoc = removevars(stim_table_dend_binoc, 'all_cells');

% contra
stim_table_dend_contra = combine_stim_contra(~contains(combine_stim_contra.all_fovs,'Soma_0'),:);
stim_table_dend_contra.mouse_cell = strcat(stim_table_dend_contra.all_mouse, '-', stim_table_dend_contra.all_cells);
stim_table_dend_contra = removevars(stim_table_dend_contra, 'all_mouse');
stim_table_dend_contra = removevars(stim_table_dend_contra, 'all_cells');

%ipsi
stim_table_dend_ipsi = combine_stim_ipsi(~contains(combine_stim_ipsi.all_fovs,'Soma_0'),:);
stim_table_dend_ipsi.mouse_cell = strcat(stim_table_dend_ipsi.all_mouse, '-', stim_table_dend_ipsi.all_cells);
stim_table_dend_ipsi = removevars(stim_table_dend_ipsi, 'all_mouse');
stim_table_dend_ipsi = removevars(stim_table_dend_ipsi, 'all_cells');

%% separate somas 
% binoc
soma_binoc = combine_stim_binoc(contains(combine_stim_binoc.all_fovs,'Soma'),:);
soma_binoc.mouse_cell = strcat(soma_binoc.all_mouse, '-', soma_binoc.all_cells);
soma_binoc = removevars(soma_binoc, 'all_mouse');
soma_binoc = removevars(soma_binoc, 'all_cells');
%retrieve only soma indices 
soma_binoc = soma_binoc(soma_binoc.all_roi_inds==1,:);

% contra
soma_contra = combine_stim_contra(contains(combine_stim_contra.all_fovs,'Soma'),:);
soma_contra.mouse_cell = strcat(soma_contra.all_mouse, '-', soma_contra.all_cells);
soma_contra = removevars(soma_contra, 'all_mouse');
soma_contra = removevars(soma_contra, 'all_cells');
soma_contra = soma_contra(soma_contra.all_roi_inds==1,:);

% ipsi
soma_ipsi = combine_stim_ipsi(contains(combine_stim_ipsi.all_fovs,'Soma'),:);
soma_ipsi.mouse_cell = strcat(soma_ipsi.all_mouse, '-', soma_ipsi.all_cells);
soma_ipsi = removevars(soma_ipsi, 'all_mouse');
soma_ipsi = removevars(soma_ipsi, 'all_cells');
soma_ipsi = soma_ipsi(soma_ipsi.all_roi_inds==1,:);


%% Add which session table belongs to
soma_ipsi.session = repmat({"ipsi"}, height(soma_ipsi),1);
soma_binoc.session = repmat({"binoc"}, height(soma_binoc),1);
soma_contra.session = repmat({"contra"}, height(soma_contra),1);

%% Retrieve aligmnent for binoc, ipsi, contra session
soma_spine_binoc = retrieve_soma_spine_alignment(stim_table_dend_binoc,soma_binoc);
soma_spine_contra = retrieve_soma_spine_alignment(stim_table_dend_contra,soma_contra);
soma_spine_ipsi = retrieve_soma_spine_alignment(stim_table_dend_ipsi,soma_ipsi);
% 
%% Add which session table belongs to
soma_spine_binoc.session = repmat({"binoc"}, height(soma_spine_binoc),1);
soma_spine_contra.session = repmat({"contra"}, height(soma_spine_contra),1);
soma_spine_ipsi.session = repmat({"ipsi"}, height(soma_spine_ipsi),1);

%% Retrieve alignment to neighboring spines 
soma_spine_contra = retrieve_neighbor_spine_alignment(input_path,soma_spine_contra);
soma_spine_ipsi = retrieve_neighbor_spine_alignment(input_path,soma_spine_ipsi);
soma_spine_binoc = retrieve_neighbor_spine_alignment(input_path,soma_spine_binoc);

%% Retrieve mean alignment to neighboring spines within 10 um
soma_spine_contra = retrieve_mean_neighbor_spine_alignment(input_path,soma_spine_contra,5);
soma_spine_ipsi = retrieve_mean_neighbor_spine_alignment(input_path,soma_spine_ipsi,5);
soma_spine_binoc = retrieve_mean_neighbor_spine_alignment(input_path,soma_spine_binoc,5);

%% Orientation preference is based on direction of movement --> switch to orientation
all_stim_table = [soma_spine_contra;soma_spine_ipsi;soma_spine_binoc];
all_day = all_stim_table.all_day;
all_stim_table.days = all_day;

%vector
Oris = all_stim_table.all_Ori_pref_vector;
Oris = Oris-90;
Oris(Oris<0) = Oris(Oris<0)+180; 
all_stim_table.all_Ori_pref_vector = Oris;

spine_save_mat = ['spine_',savefile];
save(fullfile(savepath, spine_save_mat), "all_stim_table");

%% Orientation preference is based on direction of movement --> switch to orientation
all_soma_stim_table = [soma_contra;soma_ipsi;soma_binoc];
all_day = all_soma_stim_table.all_day;
all_soma_stim_table.days = all_day;

%vector
Oris = all_stim_table.all_Ori_pref_vector;
Oris = Oris-90;
Oris(Oris<0) = Oris(Oris<0)+180; 
all_stim_table.all_Ori_pref_vector = Oris;


soma_save_mat = ['soma_',savefile];
save(fullfile(savepath, soma_save_mat), 'all_soma_stim_table');


%%
function day = which_day(day)
    if contains(day, '22') || contains(day, '23') || contains(day, '24')
        day = 'D1';
    elseif contains(day, '27')||contains(day, '28') || contains(day, '29') 
        day = 'D5';
    else
        day = 'D10';
    end
end
%% retrieve soma-spine alignment 
function table_dend = retrieve_soma_spine_alignment(table_dend,table_soma)

    soma_ori_pref_vector = [];
    soma_dir_pref_vector = [];
    soma_DSI_vector = [];
    soma_OSI_vector = [];
    spine_mean_resp = [];
    soma_mean_resp = [];


    ori_alignment_vector = [];
    dir_alignment_vector = [];
    dir_corr = [];
    soma_vis = [];
    for i = 1:height(table_dend)
        spine = table_dend(i, :);
        mouse_cell = spine.mouse_cell;
        day = spine.all_day;
        ind = find(contains(table_soma.mouse_cell, mouse_cell) & ...
            contains(table_soma.all_day, day));
        soma = table_soma(ind,:);

        % Ori diff vector
        ori_diff1 = abs(spine.all_Ori_pref_vector - soma.all_Ori_pref_vector);
        if ori_diff1 > 90
            ori_diff1 = abs(180 - ori_diff1);
        end

        % Dir diff vector
        dir_diff1 = abs(spine.all_Dir_pref_vector - soma.all_Dir_pref_vector);
        if dir_diff1 > 180
            dir_diff1 = abs(360 - dir_diff1);
        end

        % correlation b/w mean amps to 8 directions
        if isnan(spine.all_Ori_pref_vector)||isempty(soma.all_mean_amp)
            ori_alignment_vector = [ori_alignment_vector; nan];
            dir_alignment_vector = [dir_alignment_vector; nan];
                        dir_corr = [dir_corr; nan];
            soma_vis = [soma_vis; nan];
            soma_ori_pref_vector = [soma_ori_pref_vector; nan];
            soma_dir_pref_vector = [soma_dir_pref_vector; nan];
            soma_DSI_vector = [soma_DSI_vector; nan];
            soma_OSI_vector = [soma_OSI_vector; nan];
            spine_mean_resp = [spine_mean_resp; {nan}];
            soma_mean_resp = [soma_mean_resp; {nan}];

        else
            dir = corr(spine.all_mean_amp{:}',soma.all_mean_amp{:}');
            dir_alignment_vector = [dir_alignment_vector; dir_diff1];
            ori_alignment_vector = [ori_alignment_vector; ori_diff1];
            dir_corr = [dir_corr; dir];
            soma_vis = [soma_vis; soma.resp];
            soma_ori_vector = soma.all_Ori_pref_vector - 90;
            if soma_ori_vector < 0
                soma_ori_vector = soma_ori_vector + 180;
            end
            
           
            spine_mean_resp = [spine_mean_resp; spine.all_mean_amp];
            soma_mean_resp = [soma_mean_resp; soma.all_mean_amp];
            soma_ori_pref_vector = [soma_ori_pref_vector; soma_ori_vector];
            soma_dir_pref_vector = [soma_dir_pref_vector; soma.all_Dir_pref_vector];
            soma_DSI_vector = [soma_DSI_vector; soma.all_DSI_vector];
            soma_OSI_vector = [soma_OSI_vector; soma.all_OSI_vector];

                       
        end
        

    end
    table_dend = removevars(table_dend, 'all_mean_amp');
    table_dend.spine_mean_resp = spine_mean_resp;
    table_dend.soma_mean_resp = soma_mean_resp;

    table_dend.soma_DSI_vector = soma_DSI_vector;
    table_dend.soma_OSI_vector = soma_OSI_vector;
    table_dend.soma_dir_pref_vector = soma_dir_pref_vector;
    table_dend.soma_ori_pref_vector = soma_ori_pref_vector;

    table_dend.all_ori_alignment_vector = ori_alignment_vector;
    table_dend.all_dir_alignment_vector = dir_alignment_vector;
    
    table_dend.all_dir_corr = dir_corr;
    table_dend.soma_resp = soma_vis;
end
%% retrieve spine-spine alignment for signal and noise corr to neighbor
%distances resized by 0.606: added 9/28/2023
function table_dend = retrieve_neighbor_spine_alignment(distance_path,table_dend)
    co_activity_trial = [];
    co_activity_ISI = [];
    co_activity_mean_amp = [];
    co_activity_trial_amp = [];
    all_ori_offset_spine = [];
    neighbor_spine_distance = [];
    neighbor_spine_resp = [];
    
    for i = 1:height(table_dend)
        disp(i)
        
        is_blinded = 0;
         
        spine = table_dend(i, :);
        
        %get data of spines
        roi = spine.all_roi_inds;
        day = spine.all_day{:};
        mouse_cell = spine.mouse_cell;
        

        split_mouse_cell = split(mouse_cell, '-');
        mouse = split_mouse_cell{1};
        cells = split_mouse_cell{2};
        dend_segment = spine.all_fovs;
        
        %get all spines from that mouse,cell, fov
        filter = contains(table_dend.mouse_cell,mouse_cell)...
            &contains(table_dend.all_fovs, dend_segment)&...
            contains(table_dend.all_day, day);
        all_spines = table_dend(filter,:);
        split_dend = split(dend_segment, '_');
        dend = split_dend{1};

        %check if folder is blinded
        blind_folder = fullfile(distance_path, mouse, 'key_to_day_folders.csv');
        if exist(blind_folder,"file")
             is_blinded = 1;
             tbl = readtable(blind_folder);
             tbl_filt = tbl(ismember(tbl.Var1_3, cells),{'Var1_5';'Var1_4'});
             real_date = tbl_filt.Var1_4;
             blind_date = tbl_filt.Var1_5;
                 
             ind1 = find(contains(real_date, day));
             day = blind_date{ind1};
               
         end

        stats_file = fullfile(distance_path, mouse, day,cells,  ...
            dend, 'spine_stats.csv');
        distance_file = fullfile(distance_path, mouse, day,cells,  ...
            dend, 'dendritic_distance.csv');
        if exist(stats_file,"file") & ~isnan(roi)
            stats = readtable(stats_file);
            distances = table2array(readtable(distance_file)).*0.606; %resize distances
            distances_by_seg = stats.source_file;
            seg = ['RoiSet1.zip'];
            inds = find(contains(distances_by_seg,seg));
            filter_distances = distances(inds,inds);
            temp = [1:size(filter_distances,1)];
            
            
            min_distance = min(filter_distances(roi,temp~=roi));
            if ~isempty(min_distance)
                roi_inds = all_spines.all_roi_inds(temp~=roi);
                neighbor_spine_ind = roi_inds(filter_distances(roi,temp~=roi)==min_distance);
                neighbor_spine = all_spines(neighbor_spine_ind,:);
                if length(neighbor_spine_ind)>1
                    resp_neighbor = find(neighbor_spine.resp);
                    if length(resp_neighbor)>1
                        corr_spines_trial = [];
                        ori_offset_spines = [];
                        corr_spines_ISI = [];
                        corr_spines_mean_amp = [];
                        corr_spines_trial_amp = [];
                        for ii = 1:length(neighbor_spine_ind)
                            corr_spine_trial = corr(spine.all_trial_data{:}, neighbor_spine(ii,:).all_trial_data{:});
                            corr_spines_trial = [corr_spines_trial, corr_spine_trial];
                           
                            corr_spine_ISI = corr(spine.all_ISI_data{:}, neighbor_spine(ii,:).all_ISI_data{:});
                            corr_spines_ISI = [corr_spines_ISI, corr_spine_ISI];

                            corr_spine_mean_amp = corr(spine.spine_mean_resp{:}', neighbor_spine(ii,:).spine_mean_resp{:}');
                            corr_spines_mean_amp = [corr_spines_mean_amp, corr_spine_mean_amp];
                            
                            spine_trial_amp = spine.all_trial_amp{:};
                            neighbor_spine_trial_amp = neighbor_spine(ii,:).all_trial_amp{:};
                            corr_spine_trial_amp = corr(reshape(spine_trial_amp, [size(spine_trial_amp,2)*size(spine_trial_amp,3)*size(spine_trial_amp,4),1])...
                                , reshape(neighbor_spine_trial_amp, [size(neighbor_spine_trial_amp,2)*size(neighbor_spine_trial_amp,3)*size(neighbor_spine_trial_amp,4),1]));

                            corr_spines_trial_amp = [ corr_spines_trial_amp, corr_spine_trial_amp];
                           
                            diff1 = abs(spine.all_Ori_pref_vector - neighbor_spine(ii,:).all_Ori_pref_vector);
                            diff1(diff1>90) = 180 - diff1(diff1>90);
                            ori_offset_spines = [ori_offset_spines, diff1];
                        end
                        corr_spine_trial_amp = mean(corr_spines_trial_amp);
                        corr_spine_mean_amp = mean(corr_spines_mean_amp);
                        corr_spine_trial = mean(corr_spines_trial);
                        corr_spine_ISI = mean(corr_spines_ISI);
                        ori_offset_spine = mean(ori_offset_spines);

                        resp = 1;
                    elseif isempty(resp_neighbor)
                        corr_spines_trial = [];
                        ori_offset_spines = [];
                        corr_spines_ISI = [];
                        corr_spines_mean_amp = [];
                        corr_spines_trial_amp = [];
                        for ii = 1:length(neighbor_spine_ind)
                            corr_spine_trial = corr(spine.all_trial_data{:}, neighbor_spine(ii,:).all_trial_data{:});
                            corr_spines_trial = [corr_spines_trial, corr_spine_trial];
                            
                            corr_spine_ISI = corr(spine.all_ISI_data{:}, neighbor_spine(ii,:).all_ISI_data{:});
                            corr_spines_ISI = [corr_spines_ISI, corr_spine_ISI];

                            corr_spine_mean_amp = corr(spine.spine_mean_resp{:}', neighbor_spine(ii,:).spine_mean_resp{:}');
                            corr_spines_mean_amp = [corr_spines_mean_amp, corr_spine_mean_amp];
                            
                            spine_trial_amp = spine.all_trial_amp{:};
                            neighbor_spine_trial_amp = neighbor_spine(ii,:).all_trial_amp{:};
                            corr_spine_trial_amp = corr(reshape(spine_trial_amp, [size(spine_trial_amp,2)*size(spine_trial_amp,3)*size(spine_trial_amp,4),1])...
                                , reshape(neighbor_spine_trial_amp, [size(neighbor_spine_trial_amp,2)*size(neighbor_spine_trial_amp,3)*size(neighbor_spine_trial_amp,4),1]));

                            corr_spines_trial_amp = [ corr_spines_trial_amp, corr_spine_trial_amp];

                            diff1 = abs(spine.all_Ori_pref_vector - neighbor_spine(ii,:).all_Ori_pref_vector);
                            diff1(diff1>90) = 180 - diff1(diff1>90);
                            ori_offset_spines = [ori_offset_spines, diff1];

                        end
                        corr_spine_trial_amp = mean(corr_spines_trial_amp);
                        corr_spine_mean_amp = mean(corr_spines_mean_amp);
                        corr_spine_trial = mean(corr_spines_trial);
                        corr_spine_ISI = mean(corr_spines_ISI);
                        ori_offset_spine = mean(ori_offset_spines);
                        resp = 0;
                    else
                        corr_spine_trial = corr(spine.all_trial_data{:}, neighbor_spine(resp_neighbor,:).all_trial_data{:});
                        corr_spine_ISI = corr(spine.all_ISI_data{:}, neighbor_spine(resp_neighbor,:).all_ISI_data{:});
                        corr_spine_mean_amp = corr(spine.spine_mean_resp{:}', neighbor_spine(resp_neighbor,:).spine_mean_resp{:}');
                       
                        spine_trial_amp = spine.all_trial_amp{:};
                        neighbor_spine_trial_amp = neighbor_spine(resp_neighbor,:).all_trial_amp{:};
                        corr_spine_trial_amp = corr(reshape(spine_trial_amp, [size(spine_trial_amp,2)*size(spine_trial_amp,3)*size(spine_trial_amp,4),1])...
                        , reshape(neighbor_spine_trial_amp, [size(neighbor_spine_trial_amp,2)*size(neighbor_spine_trial_amp,3)*size(neighbor_spine_trial_amp,4),1]));

                        diff1 = abs(spine.all_Ori_pref_vector - neighbor_spine(resp_neighbor,:).all_Ori_pref_vector);
                        diff1(diff1>90) = 180 - diff1(diff1>90);
                        ori_offset_spine = diff1;

                        resp = 1;
                    end
                else
                    corr_spine_trial = corr(spine.all_trial_data{:}, neighbor_spine.all_trial_data{:});
                    corr_spine_ISI = corr(spine.all_ISI_data{:}, neighbor_spine.all_ISI_data{:});
                    corr_spine_mean_amp = corr(spine.spine_mean_resp{:}', neighbor_spine.spine_mean_resp{:}');
                    spine_trial_amp = spine.all_trial_amp{:};
                    neighbor_spine_trial_amp = neighbor_spine.all_trial_amp{:};
                    corr_spine_trial_amp = corr(reshape(spine_trial_amp, [size(spine_trial_amp,2)*size(spine_trial_amp,3)*size(spine_trial_amp,4),1])...
                        , reshape(neighbor_spine_trial_amp, [size(neighbor_spine_trial_amp,2)*size(neighbor_spine_trial_amp,3)*size(neighbor_spine_trial_amp,4),1]));

                    diff1 = abs(spine.all_Ori_pref_vector - neighbor_spine.all_Ori_pref_vector);
                    diff1(diff1>90) = 180 - diff1(diff1>90);
                    ori_offset_spine = diff1;
                    
                    resp = neighbor_spine.resp;
                end
                co_activity_trial = [co_activity_trial;corr_spine_trial];
                co_activity_ISI = [co_activity_ISI;corr_spine_ISI];
                co_activity_mean_amp = [co_activity_mean_amp; corr_spine_mean_amp];
                co_activity_trial_amp = [co_activity_trial_amp; corr_spine_trial_amp];
                all_ori_offset_spine = [all_ori_offset_spine; ori_offset_spine];
                
                neighbor_spine_distance = [neighbor_spine_distance; min_distance];
                neighbor_spine_resp = [neighbor_spine_resp; resp];
            else
                co_activity_trial = [co_activity_trial;NaN];
                co_activity_ISI = [co_activity_ISI;NaN];
                co_activity_mean_amp = [co_activity_mean_amp; NaN];
                co_activity_trial_amp = [co_activity_trial_amp; NaN];
                all_ori_offset_spine = [all_ori_offset_spine;NaN];
                neighbor_spine_distance = [neighbor_spine_distance; NaN];
                neighbor_spine_resp = [neighbor_spine_resp; NaN];

            end
        else
            co_activity_trial = [co_activity_trial;NaN];
            co_activity_ISI = [co_activity_ISI;NaN];
            co_activity_mean_amp = [co_activity_mean_amp; NaN];
            co_activity_trial_amp = [co_activity_trial_amp; NaN];
            all_ori_offset_spine = [all_ori_offset_spine;NaN];
            neighbor_spine_distance = [neighbor_spine_distance; NaN];
            neighbor_spine_resp = [neighbor_spine_resp; NaN];
        end


    
    %end
    
   
    end
    table_dend.ori_offset_spines = all_ori_offset_spine;
    table_dend.co_activity_mean_amp = co_activity_mean_amp;
    table_dend.co_activity_trial_amp = co_activity_trial_amp;
    table_dend.co_activity_ISI = co_activity_ISI;
    table_dend.co_activity_trial = co_activity_trial;
    table_dend.neighbor_spine_distance = neighbor_spine_distance;
    table_dend.neighbor_spine_resp = neighbor_spine_resp;
    end



%% retrieve average spine-spine alignment 
function table_dend = retrieve_mean_neighbor_spine_alignment(distance_path,table_dend, min_dist)
    co_activity_trial = [];
    co_activity_ISI = [];
    co_activity_mean_amp = [];
    co_activity_trial_amp = [];
    all_ori_offset_spine = [];

    
    for i = 1:height(table_dend)
        disp(i)
         is_blinded = 0;
         
        spine = table_dend(i, :);
        
        %get data of spines
        roi = spine.all_roi_inds;
        day = spine.all_day{:};
        mouse_cell = spine.mouse_cell;

        split_mouse_cell = split(mouse_cell, '-');
        mouse = split_mouse_cell{1};
        cells = split_mouse_cell{2};
        dend_segment = spine.all_fovs;
        
        %get all spines from that mouse,cell, fov
        filter = contains(table_dend.mouse_cell,mouse_cell)...
            &contains(table_dend.all_fovs, dend_segment)&...
            contains(table_dend.all_day, day);
        all_spines = table_dend(filter,:);
        split_dend = split(dend_segment, '_');
        dend = split_dend{1};

        %check if folder is blinded
        blind_folder = fullfile(distance_path, mouse, 'key_to_day_folders.csv');
        if exist(blind_folder,"file")
             is_blinded = 1;
             tbl = readtable(blind_folder);
             tbl_filt = tbl(ismember(tbl.Var1_3, cells),{'Var1_5';'Var1_4'});
             real_date = tbl_filt.Var1_4;
             blind_date = tbl_filt.Var1_5;
                 
             ind1 = find(contains(real_date, day));
             day = blind_date{ind1};
               
         end

        stats_file = fullfile(distance_path, mouse, day,cells,  ...
            dend, 'spine_stats.csv');
        distance_file = fullfile(distance_path, mouse, day,cells,  ...
            dend, 'dendritic_distance.csv');
        if exist(stats_file,"file") & ~isnan(roi)
            stats = readtable(stats_file);
            distances = table2array(readtable(distance_file))*0.606 % kyle's conversion factor;
            distances_by_seg = stats.source_file;
            seg = ['RoiSet1.zip'];
            inds = find(contains(distances_by_seg,seg));
            filter_distances = distances(inds,inds);
            temp = [1:size(filter_distances,1)];
           
            roi_inds = all_spines.all_roi_inds(temp~=roi);
            neighbor_spine_ind = roi_inds(filter_distances(roi,temp~=roi)<=min_dist);
            neighbor_spine = all_spines(neighbor_spine_ind,:);

            corr_spines_trial = [];
            corr_spines_ISI = [];
            corr_spines_mean_amp = [];
            corr_spines_trial_amp = [];
            ori_offset_spines = [];
            if ~isempty(neighbor_spine_ind)
                for ii = 1:length(neighbor_spine_ind)
                    corr_spine_trial = corr(spine.all_trial_data{:}, neighbor_spine(ii,:).all_trial_data{:});
                    corr_spines_trial = [corr_spines_trial, corr_spine_trial];
                    
                    corr_spine_ISI = corr(spine.all_ISI_data{:}, neighbor_spine(ii,:).all_ISI_data{:});
                    corr_spines_ISI = [corr_spines_ISI, corr_spine_ISI];

                    corr_spine_mean_amp = corr(spine.spine_mean_resp{:}', neighbor_spine(ii,:).spine_mean_resp{:}');
                    corr_spines_mean_amp = [corr_spines_mean_amp, corr_spine_mean_amp];
                    
                    
                    spine_trial_amp = spine.all_trial_amp{:};
                    neighbor_spine_trial_amp = neighbor_spine(ii,:).all_trial_amp{:};
                    corr_spine_trial_amp = corr(reshape(spine_trial_amp, [size(spine_trial_amp,2)*size(spine_trial_amp,3)*size(spine_trial_amp,4),1])...
                        , reshape(neighbor_spine_trial_amp, [size(neighbor_spine_trial_amp,2)*size(neighbor_spine_trial_amp,3)*size(neighbor_spine_trial_amp,4),1]));
                    corr_spines_trial_amp = [corr_spines_trial_amp,corr_spine_trial_amp];
                    
                    diff1 = abs(spine.all_Ori_pref_vector - neighbor_spine(ii,:).all_Ori_pref_vector);
                    diff1(diff1>90) = 180 - diff1(diff1>90);
                    ori_offset_spines = [ori_offset_spines, diff1];
    
                end
                co_activity_trial = [co_activity_trial;nanmean(corr_spines_trial)];
                co_activity_ISI = [co_activity_ISI;nanmean(corr_spines_ISI)];
                co_activity_mean_amp = [co_activity_mean_amp; nanmean(corr_spines_mean_amp)];
                co_activity_trial_amp = [co_activity_trial_amp; nanmean(corr_spines_trial_amp)];
                all_ori_offset_spine = [all_ori_offset_spine; nanmean(ori_offset_spines)];
                                                    
            else
                co_activity_trial = [co_activity_trial;NaN];
                co_activity_ISI = [co_activity_ISI;NaN];
                co_activity_mean_amp = [co_activity_mean_amp; NaN];
                co_activity_trial_amp = [co_activity_trial_amp; NaN];
                all_ori_offset_spine = [all_ori_offset_spine;NaN];
            end
        else
            co_activity_trial = [co_activity_trial;NaN];
            co_activity_ISI = [co_activity_ISI;NaN];
            co_activity_mean_amp = [co_activity_mean_amp; NaN];
            co_activity_trial_amp = [co_activity_trial_amp; NaN];
            all_ori_offset_spine = [all_ori_offset_spine;NaN];

        end


    
    %end
    
   
    end
    table_dend.mean_ori_offset_spines = all_ori_offset_spine;
    table_dend.mean_co_activity_mean_amp = co_activity_mean_amp;
    table_dend.mean_co_activity_trial_amp = co_activity_trial_amp;
    table_dend.mean_co_activity_ISI = co_activity_ISI;
    table_dend.mean_co_activity_trial = co_activity_trial;

end




end


