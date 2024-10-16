function [ind_corrs] = get_thresh_corrs(corrs, median_resp_corrs)
ind_corrs = [];
for i = 1:length(corrs)
    inds = corrs(i,:)>median_resp_corrs(i);
    ind_corrs = [ind_corrs;inds];
end
end