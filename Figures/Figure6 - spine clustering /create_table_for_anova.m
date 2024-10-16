function [anova_tbl] = create_table_for_anova(all_combo_dists,combinations, unique_day)
all_data = [];
all_eye_combo = []
all_day = [];
for i = 1:length(unique_day)
    combo_dists_by_day = all_combo_dists{i};
    for c = 1:length(combinations)
        pair = split(combinations(c), '_');
        pair1 = pair{1};
        pair2 = pair{2};
        string_pair = ['pair', pair1,'to',pair2];
        all_pair_dists = cellfun(@(x) x.(string_pair),combo_dists_by_day, 'UniformOutput',false);
        all_pair_dists = cell2mat(all_pair_dists');

        all_data = [all_data; all_pair_dists];
        all_day = [all_day; repmat(unique_day(i), size(all_pair_dists))];
        all_eye_combo = [all_eye_combo; repmat({string_pair}, size(all_pair_dists))];


    end
  
end
anova_tbl = table(all_data,all_day,all_eye_combo);
end