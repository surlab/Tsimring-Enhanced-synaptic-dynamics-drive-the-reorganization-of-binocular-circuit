function [mean_grids] = shuffle_corrs(bins,corrs_by_branch, dists)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
num_dends = length(corrs_by_branch);
all_corrs = [];
for ii = 1:num_dends
    corrs = corrs_by_branch{ii};
    if ~isempty(corrs)
        temp_corrs = corrs(randperm(length(corrs)));
        all_corrs = [all_corrs; temp_corrs];
    end
end


[mean_grids, ~] = bin_corrs_output_mean_std(bins, all_corrs, dists);
end
