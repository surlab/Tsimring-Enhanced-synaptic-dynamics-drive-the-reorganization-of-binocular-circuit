function all_mean_grids = compare_model_data_delta_corr(corrs_by_days,dists_by_days, bins)

    all_mean_grids = [];
    for d = 1:length(corrs_by_days)
        corrs_by_branch = corrs_by_days{d};
        dists_by_branch = dists_by_days{d};
        if iscell(corrs_by_branch)
            [mean_grids, ~, ~] = bin_corrs_output_mean_std(bins...
                , cell2mat(corrs_by_branch), cell2mat(dists_by_branch) );
        else
             [mean_grids, ~, ~] = bin_corrs_output_mean_std(bins...
                , corrs_by_branch, dists_by_branch);
        end
        all_mean_grids = [all_mean_grids;mean_grids];

    end
end
