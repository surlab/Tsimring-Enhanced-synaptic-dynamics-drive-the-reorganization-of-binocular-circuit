function [ps] = measure_p_val_rank_sum(bins, corrs, dists, shuffled_means)

ps = [];  
for i = 1:length(bins)-1
    ind = dists>=bins(i)&dists<bins(i+1);
    compare_corr = corrs(ind);
    mean_corr = mean(compare_corr);
    diff_mean = shuffled_means(i,:)-mean_corr;
    ps = [ps min(sum(diff_mean > 0), sum(diff_mean < 0))/length(shuffled_means)];
end
end

