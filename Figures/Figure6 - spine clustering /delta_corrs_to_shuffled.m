function [mean_corrs, stds_corrs] = delta_corrs_to_shuffled(bins, corrs, dists, shuffled_means)
all_corrs = zeros(length(dists),2);

mean_corrs = [];
stds_corrs = [];
for i = 1:length(bins)-1
    inds = dists>=bins(i)&dists<bins(i+1);
    n_corrs = corrs(inds);
    n_corrs = n_corrs(~isnan(n_corrs));
    mean_corrs = [mean_corrs;mean(n_corrs-mean(shuffled_means(i,:)))];
    stds_corrs = [stds_corrs;std(n_corrs-mean(shuffled_means(i,:)))/sqrt(length(n_corrs))];

end
end
