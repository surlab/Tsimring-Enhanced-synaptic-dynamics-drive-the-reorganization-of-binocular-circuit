function plot_pair_dists(all_combo_dists,combinations, type_inds,unique_day)
xticks_count = [];
xticks_label = []
for i = 1:length(unique_day)
    count = i;
    combo_dists_by_day = all_combo_dists{i};
    for c = 1:length(combinations(type_inds))
        pair = split(combinations(type_inds(c)), '_');
        pair1 = pair{1};
        pair2 = pair{2};
        string_pair = ['pair', pair1,'to',pair2];
        all_pair_dists = cellfun(@(x) x.(string_pair),combo_dists_by_day, 'UniformOutput',false);
        all_pair_dists = cell2mat(all_pair_dists');
        %swarmchart(count*ones(size(all_pair_dists)),all_pair_dists); hold on;
        bar(count,mean(all_pair_dists));hold on;
        errorbar(count,mean(all_pair_dists),std(all_pair_dists)/sqrt(length(all_pair_dists)), '.k'); 
        scatter(count*ones(size(all_pair_dists)), all_pair_dists, '.k', 'jitter', 'on','jitterAmount', 0.1);
        if i ==1
            text(count+(i/2),10,string_pair)
        end
        xticks_count = [xticks_count, count];
        xticks_label = [xticks_label, unique_day(i)];
        count = count + 3;
    end
end
[~, inds] = sort(xticks_count);

xticks(xticks_count(inds));
xticklabels(xticks_label(inds));
ylim([0,35])
ylabel('Spine Pair Distance')