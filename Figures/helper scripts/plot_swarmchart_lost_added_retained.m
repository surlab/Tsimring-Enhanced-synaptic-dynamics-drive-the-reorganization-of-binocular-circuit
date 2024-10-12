function plot_swarmchart_lost_added_retained(data_feature_not_retained,data_feature_retained, ylim_axis,y_label)
    mean_data_not_retained = cellfun(@(x) nanmean(x), data_feature_not_retained);
    err_data_not_retained = cellfun(@(x) nanstd(x)/sqrt(length(x)), data_feature_not_retained);
    
    mean_data_retained = cellfun(@(x) nanmean(x), data_feature_retained);
    err_data_retained = cellfun(@(x) nanstd(x)/sqrt(length(x)), data_feature_retained);
    
    all_data = [data_feature_not_retained; data_feature_retained];
    mean_data = [mean_data_not_retained; mean_data_retained];
    error_data = [err_data_not_retained; err_data_retained];
    
    %D1 to D5
    subplot(1,2,1)
    %bar(mean_data(:,[1,2])'); hold on; 
    plot_errorbar_swarm(mean_data(:,[1,2])',error_data(:,[1,2])', all_data([1,2],[1,2])', [1,2;4,5])
    [p,a]= ranksum(all_data{1,1}, all_data{2,1})
    plot_sig(1.5,max(ylim_axis),p);
    [p,a]= ranksum(all_data{1,2}, all_data{2,2})
    plot_sig(4.5,max(ylim_axis),p);
    set(gca, 'FontSize',10)
    ylim(ylim_axis)
    title(['D1 to D5'])
    legend({'not retained', 'retained'})
    xticklabels({'D1 lost', 'D5 formed'})  
    ylabel(y_label)
    
    %D5 to D10
    subplot(1,2,2)
    %bar(mean_data(:,[3,4])'); hold on;
    plot_errorbar_swarm(mean_data(:,[3,4])',error_data(:,[3,4])',all_data([1,2],[3,4]), [1,2;4,5])
    [p,a]= ranksum(all_data{1,3}, all_data{2,3})
    plot_sig(1.5,max(ylim_axis),p);
    [p,a]= ranksum(all_data{1,4}, all_data{2,4})
    plot_sig(4.5,max(ylim_axis),p);
    set(gca, 'FontSize',10)
    ylim(ylim_axis)
    title(['D5 to D10'])
    legend({'not retained', 'retained'})
    xticklabels({'D5 lost', 'D10 formed'})
end

