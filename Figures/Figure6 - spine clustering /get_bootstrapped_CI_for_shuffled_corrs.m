function [CI_by_session, shuffled_means_by_sessions, shuffled_dists_by_sessions] = get_bootstrapped_CI_for_shuffled_corrs(bins, shuffle_times, all_corrs_by_session, all_resp_dists_by_session, all_dists_by_session)

sessions = fieldnames(all_corrs_by_session);
for s = 1:length(sessions)
    CI_by_day = [];
    shuffled_means_by_day = [];
    shuffled_dists_by_day = [];
    corrs_by_days = all_corrs_by_session.(sessions{s});
    resp_dists_by_days = all_resp_dists_by_session.(sessions{s});
    dists_by_days = all_dists_by_session.(sessions{s});
    for d = 1:length(corrs_by_days)
        corrs_by_branch = corrs_by_days{d};
        resp_dists_by_branch = resp_dists_by_days{d};
        dists_by_branch = dists_by_days{d};
        shuffled_means = zeros(length(bins)-1, shuffle_times);
        shuffled_dists = [];
        for t = 1:shuffle_times
            shuffled_dists = [shuffled_dists, get_shuffle_dists(dists_by_branch,corrs_by_branch)];
            shuffle_mean_grids = shuffle_corrs(bins, corrs_by_branch, cell2mat(resp_dists_by_branch));
            shuffled_means(:,t) = shuffle_mean_grids;
        end
        CIs = sort(shuffled_means,2);
        CI_97_5 = CIs(:,0.975*10000);
        CI_2_5 = CIs(:,0.025*10000);
        CI_by_day = [CI_by_day; {CI_2_5, CI_97_5}];
        shuffled_means_by_day = [shuffled_means_by_day; {shuffled_means}];
        shuffled_dists_by_day = [shuffled_dists_by_day; {shuffled_dists}];


    end
    CI_by_session.(sessions{s}) = CI_by_day;
    shuffled_means_by_sessions.(sessions{s}) =  shuffled_means_by_day;
    shuffled_dists_by_sessions.(sessions{s}) = shuffled_dists_by_day;

end
