function ind = find_preferred_SF(mean_amplitude, t_test)

criteria = t_test <= 0.05 & mean_amplitude >= 0.5;
if size(criteria,1) == 1
    find_sum_criteria = find(sum(criteria));
else
    find_sum_criteria = find(sum(sum(criteria)));
end
if isempty(find_sum_criteria)
    [max_am,ind_am] = max(mean_amplitude(1,:,:),[],3);
    [min_t,ind_t] = min(t_test(1,:,:),[],3);
    [~,ind2_am] = max(max_am,[],2);
    [~,ind2_t] = min(min_t,[],2);
    if ind2_t == ind2_am
        ind = ind_am(ind2_am);
    else
        ind = ind_t(ind2_t);
    end    
elseif length(find_sum_criteria) == 1
    ind = find_sum_criteria;                          
else
    soma_mean_amp = mean_amplitude(1,:,:);
    [max_am, ind1] = max(soma_mean_amp(:,:, find_sum_criteria), [], 3);
    [~, ind2] = max(max_am, [], 2);
    ind = find_sum_criteria(ind1(ind2));
    if criteria(1,ind2, ind) == 0
        k = 1;
    end
end

