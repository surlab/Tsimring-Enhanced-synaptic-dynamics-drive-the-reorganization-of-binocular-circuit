function [mean_corrs] = get_mean_corrs_matrix_exclude_diagonal(corrs)
mean_corrs = [];

 for i = 1:length(corrs)
    inds = ones(size(corrs,2),1);
    inds(i) =0;
    
    mean_corr = mean(corrs(i,find(inds)));
    mean_corrs = [mean_corrs;mean_corr];
 end
end