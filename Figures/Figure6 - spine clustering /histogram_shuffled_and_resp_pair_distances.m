function histogram_shuffled_and_resp_pair_distances(bins,all_resps_dists_by_session,shuffled_dists_by_session, days)

sessions = fieldnames(all_resps_dists_by_session);
for s = 1:length(sessions)
    figure("Name",sessions{s})
    resp_dists_by_days = all_resps_dists_by_session.(sessions{s});
    dists_by_days = shuffled_dists_by_session.(sessions{s});
    for d = 1:length(days)
        subplot(length(days),1,d)
        resp_dists_by_branch = resp_dists_by_days{d};
        dists_by_branch = dists_by_days{d};
       
        all_shuffle = cell2mat(dists_by_branch);
        median_resp_dists = median(cell2mat(resp_dists_by_branch));
        median_shuffle_dists = median(all_shuffle,1);
        p = min(sum(median_resp_dists>median_shuffle_dists),sum(median_resp_dists<median_shuffle_dists))/size(all_shuffle,2)
        histogram(median_shuffle_dists, bins, 'Normalization','probability'); hold on;
        plot([median_resp_dists,median_resp_dists],[0,0.2], '-k', 'LineWidth',1);
      
        plot_sig(mean(bins),0.25,p);
        xlabel('Spine Pair Distance (um)')
        ylabel('Fraction')
        title(days{d})
        ylim([0,0.25])
     end
end
       
