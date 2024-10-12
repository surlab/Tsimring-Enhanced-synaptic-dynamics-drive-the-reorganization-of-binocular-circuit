function x = plot_errorbar_scatter(y, err,data)
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i),'.k');   
    for ii = 1:length(x)
        %scatter(x(ii)*ones(size(data{ii,i})), data{ii,i}, '.k', 'jitter', 'on', 'jitterAmount', 0.05);
    end
end
hold off
end
