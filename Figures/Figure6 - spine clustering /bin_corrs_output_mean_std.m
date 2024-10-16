function [mean_grid, std_grid, corr_grids] = bin_corrs_output_mean_std(bins, corrs, dists )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
mean_grid = [];
std_grid = [];
corr_grids = [];
for i = 1:length(bins)-1
    inds = dists>=bins(i)&dists<bins(i+1);
    n_corrs = corrs(inds);
   n_corrs = n_corrs(~isnan(n_corrs));
   corr_grids = [corr_grids, {n_corrs}];
   mean_grid = [mean_grid,mean(n_corrs)];
   std_grid = [std_grid, std(n_corrs)/sqrt(length(n_corrs))];
 
end

end

