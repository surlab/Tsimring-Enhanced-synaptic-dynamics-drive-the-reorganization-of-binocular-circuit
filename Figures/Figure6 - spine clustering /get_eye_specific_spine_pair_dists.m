%% Obtain spine-spine correlations or offset for spine pairs in dendritic segment
function [all_resp_dists] =  get_eye_specific_spine_pair_dists(all_stim_table, unique_day, distance_path, combin)

    all_resp_dists = [];
    for d = 1:length(unique_day)
        resp_distances = [];
        table_mouse = all_stim_table(ismember(all_stim_table.days, unique_day(d)),:);
        unique_mouse_cell = unique(table_mouse.mouse_cell);
            
        for m = 1:length(unique_mouse_cell)
            split_mouse_cell = split(unique_mouse_cell(m), '-');
            mouse = split_mouse_cell{1};
            cell = split_mouse_cell{2};
            inds = ismember(table_mouse.mouse_cell, unique_mouse_cell(m));
            table_dend = table_mouse(inds,:);
            unique_dend = unique(table_dend.all_fovs);
            is_blinded = 0;
            blind_folder = fullfile(distance_path,mouse, 'key_to_day_folders.csv');
            if exist(blind_folder,"file")
                is_blinded = 1;
                tbl = readtable(blind_folder);
               
            end
            for dd = 1:length(unique_dend)
                dend_segment =  unique_dend(dd);
                inds = ismember(table_dend.all_fovs, unique_dend(dd));
                table_spines = table_dend(inds,:);
                split_dend = split(dend_segment, '_');
                dend = split_dend{1};
                day = table2array(table_spines(1,"all_day"));
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
                if exist(stats_file,"file") & ~isnan(table2array(table_spines(1,"all_roi_inds")))
                    stats = readtable(stats_file);
                    distances = table2array(readtable(distance_file))*0.606;
                    distances_by_seg = stats.source_file;
                    seg = ['RoiSet', num2str(str2num(split_dend{2})+1), '.zip'];
                    inds = find(contains(distances_by_seg,seg));
                    filter_distances = distances(inds,inds);
                    
                    for c = 1:length(combin)
                        pair = split(combin{c}, '_');

                        spine_ind1 = find(strcmp(table_spines.ci_spine, pair(1)));
                        spine_ind2 = find(strcmp(table_spines.ci_spine, pair(2)));
                        resp_filter_distances = filter_distances(spine_ind1, spine_ind2);
                        % if pair is the same type
                        if pair{1} == pair{2}
                            mask = triu(true(size(resp_filter_distances)),1);
                            combin_pair_distances.(['pair',pair{1}, 'to',pair{2}]) = resp_filter_distances(mask);

                        else
                           combin_pair_distances.(['pair',pair{1}, 'to',pair{2}]) = resp_filter_distances(:);

                        end
                  
                    end
                    all_pairs_mask = triu(true(size(filter_distances)),1);
                    combin_pair_distances.all_dists = filter_distances(all_pairs_mask);


                    
                end
                resp_distances = [resp_distances, {combin_pair_distances}];
            end
        end
        all_resp_dists = [all_resp_dists, {resp_distances}];
    end
end
