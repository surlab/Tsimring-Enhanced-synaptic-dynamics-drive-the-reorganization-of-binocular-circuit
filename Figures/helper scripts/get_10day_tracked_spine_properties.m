function [p,stats, tbl] = get_10day_tracked_spine_properties(feature, all_temp, unique_sessions, xtick_labels)
for s = 1:length(unique_sessions)
    data_feature = [];
    for i = 1:length(all_temp)
        temp  = all_temp{i};
        temp = temp(strcmp(temp.session_D1_D5,unique_sessions{s}), :).D5_D1_D5;
        if ismember(feature, 'all_dir_corr')|contains(feature, 'all_ori_alignment')
            temp = temp(temp.soma_resp>0&temp.resp>0,:);
        elseif contains(feature, 'mean')
            temp = temp(temp.resp>0,:);
        elseif contains(feature, 'co_activity') | contains(feature, 'ori_offset')
            temp = temp(temp.neighbor_spine_distance<5&temp.resp>0,:);
        elseif contains(feature, 'OSI')
            temp = temp(temp.resp>1,:);
        end
        data_feature = [data_feature, {temp.(feature)}];
    end
   
   session.(unique_sessions{s}) = data_feature;
   
end

figure
set(gcf, 'Position', [440,242,750,555])
for i = 1:length(unique_sessions)
    subplot(length(unique_sessions),1,i)
    data_feature = session.(unique_sessions{i});
    mean_data= cellfun(@(x) nanmean(x), data_feature);
    err_data = cellfun(@(x) nanstd(x)/sqrt(length(x)), data_feature);


    bar(mean_data); hold on; 
    errorbar([1:length(mean_data)], mean_data, err_data,'.k');
    for ii = 1:length(data_feature)
        scatter(ii*ones(size(data_feature{ii})), data_feature{ii},5, 'k', 'filled', 'jitter', 'on', 'jitterAmount', 0.05);
    end
    if i ==length(unique_sessions)
    xticklabels(xtick_labels)
    end
    title(unique_sessions{i})
    [p,tbl, stats] = perform_kruskallis_wallis(session)

end
%% Kruskallis wallis test 
function [p, tbl, stats] = perform_kruskallis_wallis(session)
fields = fieldnames(session);
structure = {"retained D1-D10","retained D5, lost D10", "added D5, lost D10", "added D5 retained D10"};
for i = 1:length(fields)
    data = session.(fields{i});
    temp = vertcat(data{:});

    concat_structure = cellfun(@(x,y)repmat(x, length(y),1),structure, data, 'UniformOutput', false);
    concat_structure = vertcat(concat_structure{:});
    [p,tbl, stats] = kruskalwallis(temp, concat_structure)
    [results,~,~,gnames] = multcompare(stats);
    tbl = array2table(results,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl.("Group A")=gnames(tbl.("Group A"));
    tbl.("Group B")=gnames(tbl.("Group B"))
end