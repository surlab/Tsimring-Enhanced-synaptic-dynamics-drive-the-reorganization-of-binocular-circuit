function [all_spine_density_table,all_mouse_table] = plot_spine_density(distance_path, temp, unique_day, get_resp, type_resp, resp)

unique_mouse_cell = unique(temp.mouse_cell);
all_spine_density_table = [];
all_mouse_table = [];
for i = 1:length(unique_mouse_cell)
    spine_density_table = NaN(1,3);
    mouse_cell = unique_mouse_cell{i};
    split_mouse_cell = split(mouse_cell, '-');
    mouse = split_mouse_cell{1};
    cells = split_mouse_cell{2};
    vals = temp(ismember(temp.mouse_cell,mouse_cell),:);
    for d = 1:length(unique_day)
        rois = vals(ismember(vals.days,unique_day{d}),:); 
        
        if ~isempty(rois)
            
            dends = table2array(rois(:,"all_fovs"));
            temp_dends = cellfun(@(x) strsplit(x, '_'), dends, 'UniformOutput', false);
            branch = cellfun(@(x) x{:,1}, temp_dends, 'UniformOutput', false); 
            branch_seg = cellfun(@(x) x{:,2}, temp_dends, 'UniformOutput', false); 
            temp_dends = string(branch);
                
            unique_dends = unique(branch);
            day = unique(rois.all_day);
            for dd = 1:length(unique_dends)
              
              blind_folder = fullfile(distance_path, mouse, 'key_to_day_folders.csv');
              if exist(blind_folder,"file")
                 is_blinded = 1;
                 tbl = readtable(blind_folder);
                 tbl_filt = tbl(ismember(tbl.Var1_3, cells),{'Var1_5';'Var1_4'});
                 real_date = tbl_filt.Var1_4;
                 blind_date = tbl_filt.Var1_5;
                     
                 ind1 = find(contains(real_date, day));
                 day = blind_date(ind1);
              end
              stats_file = fullfile(distance_path, mouse, cells, day{:}, unique_dends{dd}, 'dendritic_distance', 'spine_stats.csv');
              distance_file = fullfile(distance_path, mouse, cells, day{:}, unique_dends{dd}, 'dendritic_distance', 'dendritic_distance.csv');
              
              if exist(stats_file,"file")
                branch = readtable(stats_file);
                distances_by_seg = branch.source_file;
                num_seg = unique(distances_by_seg);
                distances = table2array(readtable(distance_file))*0.606; % kyle's conversion factor;
                total_spines = 0;
                total_dist = 0;
                for n = 1:length(num_seg)
                    ind = find(strcmp(distances_by_seg, num_seg{n}));
                    max_dist = max(max(distances(ind, ind)));
                    
                    if(get_resp)
                        branch_type = num_seg{n};
                        branch_type = str2num(branch_type(end-4));
                        spines = rois(strcmp(rois.all_fovs,[unique_dends{dd},'_', num2str(branch_type-1)]),:);
                        num_spines = 0;
                        if isfloat(type_resp)
                            num_spines = sum(spines.(resp)==type_resp);
                        else
                            for ii = 1:length(type_resp)
                                num_spines = num_spines+sum(strcmp(spines.(resp),type_resp(ii)));
                            end
                        end
                    else
                        num_spines = length(ind);
                    end
                    total_dist = total_dist + max_dist;
                    total_spines = num_spines+total_spines;
                   
               
                end
              end
            end
            spine_density_table(d) = total_spines/total_dist;
            mouse_table = {mouse_cell};
        end
end
all_spine_density_table = [all_spine_density_table; spine_density_table];
all_mouse_table = [all_mouse_table; mouse_table];
end

end
