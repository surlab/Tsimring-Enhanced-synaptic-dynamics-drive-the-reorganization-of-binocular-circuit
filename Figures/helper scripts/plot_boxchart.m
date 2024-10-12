function x = plot_boxchart(all_data, x)
    ngroups = size(all_data, 1);
    nbars = size(all_data, 2);
    all_data_temp = [];
    all_x = [];
    for i = 1:ngroups
       for ii = 1:nbars
            all_data_temp = [all_data_temp; all_data{i,ii}];
            all_x = [all_x; repmat(x(i,ii), size(all_data{i,ii}))];
       end
    end
    boxchart(all_x,all_data_temp, 'JitterOutliers','on', 'MarkerStyle','.');
    hold off
end
