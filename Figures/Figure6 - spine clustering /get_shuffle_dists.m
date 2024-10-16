function shuffled_dists = get_shuffle_dists(dists_by_branch,corrs_by_branch)
num_dends = length(corrs_by_branch);
shuffled_dists = [];
for ii = 1:num_dends
    num_resp_pairs = length(corrs_by_branch{ii});
    distances = dists_by_branch{ii};
    all_pairs = length(distances);
    shuffled_dists = [shuffled_dists; {distances(randperm(all_pairs,num_resp_pairs))}];
end

end