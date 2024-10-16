function compare_corrs_between_days(all_corrs_by_session,all_dists_by_session, bins,unique_days)
sessions = fieldnames(all_corrs_by_session);
for s = 1:length(sessions)
    corrs_by_days = all_corrs_by_session.(sessions{s});
    dists_by_days = all_dists_by_session.(sessions{s});
    all_days = [];
    all_bins = [];
    all_corr_grids = [];
    days = [1,3];
    for d = days
        corrs_by_branch = corrs_by_days{d};
        dists_by_branch = dists_by_days{d};
        
        
        [~, ~, corr_grids] = bin_corrs_output_mean_std(bins...
            , cell2mat(corrs_by_branch), cell2mat(dists_by_branch) );
        all_days = [all_days; repmat(unique_days(d),size(cell2mat(corr_grids')))];
        for i = 1:length(corr_grids)
            all_bins = [all_bins; repmat(bins(i), size(corr_grids{i}))];
        end
        all_corr_grids = [all_corr_grids;cell2mat(corr_grids')];
    end
    [p,stats,tbl] = anovan(all_corr_grids, {all_days, all_bins}, 'model', 'interaction');
    [results,~,~,gnames]=multcompare(tbl, 'Dimension', [1,2])
    tbl2 = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl2.("Group A")=gnames(tbl2.("Group A"));
    tbl2.("Group B")=gnames(tbl2.("Group B"))
end