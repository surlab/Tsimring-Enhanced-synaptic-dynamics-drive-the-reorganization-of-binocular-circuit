function tbl = find_nearest_vis_resp_neighbor(tracked_temp,all_spine_temp,d, distance_path)
structure_type = [];
all_distances = [];
for i = 1:height(tracked_temp)
        mouse_cell = tracked_temp(i,:).all_mice_cells;
        split_mouse_cell = split(mouse_cell, '-');
        mouse = split_mouse_cell{1};
        cell = split_mouse_cell{2};
        fov = tracked_temp(i,:).all_fovs;
        day = tracked_temp(i,:).all_day;
        roi = tracked_temp(i,:).all_roi_inds;
        table_spines = all_spine_temp(strcmp(all_spine_temp.mouse_cell, mouse_cell)&...
            strcmp(all_spine_temp.all_fovs, fov)&strcmp(all_spine_temp.all_day, day),:);
        spine_ind = table_spines.resp>0;
        is_blinded = 0;
        blind_folder = fullfile(distance_path,mouse, 'key_to_day_folders.csv');
        if exist(blind_folder,"file")
            is_blinded = 1;
            tbl = readtable(blind_folder);
           
        end
        split_dend = split(fov, '_');
        dend = split_dend{1};
        day = day{:};

        if is_blinded
             tbl_filt = tbl(ismember(tbl.Var1_3, cell),{'Var1_5';'Var1_4'});
             real_date = tbl_filt.Var1_4;
             blind_date = tbl_filt.Var1_5;
             ind1 = find(contains(real_date, day));
             day = blind_date{ind1};
       
        end
        stats_file = fullfile(distance_path, mouse, cell, day, ...
            dend, 'dendritic_distance', 'spine_stats.csv');
        distance_file = fullfile(distance_path, mouse, cell, day, ...
            dend, 'dendritic_distance', 'dendritic_distance.csv');
        
        stats = readtable(stats_file);
        distances = table2array(readtable(distance_file))*0.606;
        distances_by_seg = stats.source_file;
        seg = ['RoiSet', num2str(str2num(split_dend{2})+1), '.zip'];
        inds = find(contains(distances_by_seg,seg));
        filter_distances = distances(inds,inds);
        idx=ismember(1:length(filter_distances),roi);
        roi_distances = filter_distances(roi,~idx);
        
        resp_filter_distances = roi_distances(spine_ind(~idx));
        min_distance = NaN;
        if ~isempty(resp_filter_distances)>0
            min_distance = min(resp_filter_distances);
        end
        
        all_distances = [all_distances; min_distance];
end
tbl = tracked_temp;
tbl.distances = all_distances;
end