%% Obtain spine-spine correlations or offset for spine pairs in dendritic segment
function [all_mouse_cell_by_session] =  get_spine_pair_correlations_by_dist_matrix_version(feature,soma_feature,get_soma_resp, all_stim_table, unique_sessions, unique_day, distance_path)
temp = all_stim_table;
if get_soma_resp == 1
    temp = temp(temp.soma_resp>0,:);
end
spine_criteria_func = @(x,y) x>0; %& y<=0.5;
temp.session = [temp.session{:}]';

for s = 1:length(unique_sessions)

    all_mouse_cell_by_day = [];
    for d = 1:length(unique_day)
        
        
        all_mouse_cell = [];
        
        inds = contains(temp.session, unique_sessions(s)) ...
            & ismember(temp.days, unique_day(d));
        table_mouse = temp(inds,:);
        unique_mouse_cell = unique(table_mouse.mouse_cell);
            
        for m = 1:length(unique_mouse_cell)
            corr_activity = [];
            mean_soma_resp_corrs = [];
            corr_resp_distances = [];
            all_distances = [];
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
                spine_ind = find(spine_criteria_func(table_spines.resp,table_spines.all_shaft_corr));
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
                    resp_filter_distances = filter_distances(spine_ind, spine_ind);
                    spines_activity = table2array(table_spines(spine_ind,feature));
                    spine_soma = table2array(table_spines(spine_ind,soma_feature));
                    mean_spine_soma = []; 
                    corr_spine = [];
                    if length(spines_activity)>2
                        if contains(feature, 'all_Ori')
                            corr_spine = find_ori_offset_matrix(spines_activity);
                        else
                            if ismember(feature,'spine_mean_resp')
                                 spines_activity = cellfun(@(x) x', spines_activity, 'UniformOutput', false);
                            elseif ismember(feature, 'all_trial_amp')|ismember(feature,'all_active_spine_trials_smooth')
                                 spines_activity = cellfun(@(x) reshape(x,[size(x,2)*size(x,3)*size(x,4),1]), spines_activity, 'UniformOutput', false);
                            end
                            
                            corr_spine = corr([spines_activity{:}], [spines_activity{:}]);
                            mean_spine_soma = (spine_soma + spine_soma')/2; %pairwise mean for each spine corr to soma
                        end
                    end
                    all_mask = triu(true(size(filter_distances)),1);
  
                    corr_activity = [corr_activity; {corr_spine}];
                    corr_resp_distances = [corr_resp_distances; {resp_filter_distances}];
                    mean_soma_resp_corrs = [mean_soma_resp_corrs; {mean_spine_soma}];
                    all_distances = [all_distances;{filter_distances}];
                end
            end
            inds = cellfun(@(x) ~isempty(x), corr_activity);
            corr_activity = corr_activity(inds);
            corr_resp_distances = corr_resp_distances(inds);
            all_distances = all_distances(inds);
            if ~isempty(find(inds))
                mouse_cell.all_resp_dists = corr_resp_distances;
                mouse_cell.all_dists = all_distances;
                
                mouse_cell.all_corrs = corr_activity;
            end
            mouse_cell.num_dends = length(unique_dend);
            mouse_cell.soma_resp = table_spines(1,:).soma_resp;
            mouse_cell.mouse_cell = unique_mouse_cell{m};
            all_mouse_cell = [all_mouse_cell, {mouse_cell}];
            clear mouse_cell inds corr_activity corr_resp_distances all_distances

        end
        all_mouse_cell_by_day = [all_mouse_cell_by_day, {all_mouse_cell}];

    end
    all_mouse_cell_by_session.(unique_sessions{s}) = all_mouse_cell_by_day;
end
end

%% find orientation offset 
function corr_spine = find_ori_offset_matrix(corr_activity)
    corr_spine = NaN(length(corr_activity),length(corr_activity));
    for i = 1:length(corr_spine)
        for ii = 1:length(corr_spine)
            diff1 = abs(corr_activity(i)-corr_activity(ii));
            if diff1>90
                diff1 = 180 - diff1;
            end
            corr_spine(i,ii) = diff1;
        end
    end
end
