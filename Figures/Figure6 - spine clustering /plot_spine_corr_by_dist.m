function plot_spine_corr_by_dist(f,g,h,all_corrs_by_session,all_dists_by_session,CI_by_session,shuffled_means_by_sessions, bins,days, feature)
count = 1;
sessions = fieldnames(all_corrs_by_session);
all_mean_corrs = [];
all_std_corrs = [];
for s = 1:length(sessions)
    corrs_by_days = all_corrs_by_session.(sessions{s});
    dists_by_days = all_dists_by_session.(sessions{s});
    CI_by_days = CI_by_session.(sessions{s});
    shuffled_mean_by_days = shuffled_means_by_sessions.(sessions{s});
    all_mean_grids = [];
    for d = 1:length(days)
        corrs_by_branch = corrs_by_days{d};
        dists_by_branch = dists_by_days{d};
        shuffled_means = shuffled_mean_by_days{d};
        CIs = [CI_by_days{d,:}];
        
        [mean_grids, std_grids, ~] = bin_corrs_output_mean_std(bins...
            , cell2mat(corrs_by_branch), cell2mat(dists_by_branch) );
        ps = measure_p_val_rank_sum(bins, cell2mat(corrs_by_branch), ...
            cell2mat(dists_by_branch), shuffled_means);

        figure(f)
        subplot(length(days),length(sessions), count);
        plot(bins(2:end), CIs(:,1), 'Color',[0.93,0.69,0.13]); hold on;
        plot(bins(2:end), CIs(:,2), 'Color',[0.93,0.69,0.13]);
        plot(bins(2:end), mean_grids); hold on;
        errorbar(bins(2:end), mean_grids, std_grids);
        all_mean_grids = [all_mean_grids; mean_grids]
        
        plot_p = 0.45;
        for p = 1:length(ps)
             if ps(p)<0.05
                scatter(bins(p+1),plot_p, '*', 'k');
             end
        end
        xlabel('spine pair distance')
        ylabel(feature)
        title([sessions{s}, ', ',   days{d}]);
       
        ylim([0,0.5])

        xlim([4,21])
       
        figure(g)
        subplot(length(sessions),length(days), count);
        scatter(cell2mat(dists_by_branch), cell2mat(corrs_by_branch));
        xlabel('spine pair distance')
        ylabel(feature)
        title([sessions{s}, ', ',   days{d}]);
  

         count = count +1;
        
    end
end
