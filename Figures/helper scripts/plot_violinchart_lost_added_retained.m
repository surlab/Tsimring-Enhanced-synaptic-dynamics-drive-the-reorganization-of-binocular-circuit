function plot_violinchart_lost_added_retained(data_feature_not_retained,data_feature_retained,y_label,ylims)
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
    %y = 100;
    colors = [1.00,0.73,0.59;0.8,0.8,0.8;1.00,0.59,0.80;0.8,0.8,0.8]
    plot_violinplot(all_data([1,2],[1,2])', [1,2;3,4], colors)
    [p,a]= ranksum(all_data{1,1}, all_data{2,1})
    plot_sig(1.5,ylims(2),p);
    [p,a]= ranksum(all_data{1,2}, all_data{2,2})
    plot_sig(3.5,ylims(2),p);
    set(gca, 'FontSize',10)
    ylim(ylims)
    title(['D1 to D5'])
    xticks([1.5,3.5])
    xticklabels({'D1','D5'})  
    ylabel(y_label)
    
    %D5 to D10
    subplot(1,2,2)
   % y =100;
    %bar(mean_data(:,[3,4])'); hold on;
    plot_violinplot(all_data([1,2],[3,4])', [1,2;3,4], colors)
    [p,a]= ranksum(all_data{1,3}, all_data{2,3})
    plot_sig(1.5,ylims(2),p);
    [p,a]= ranksum(all_data{1,4}, all_data{2,4})
    plot_sig(3.5,ylims(2),p);
    set(gca, 'FontSize',10)
    ylim(ylims)
    title(['D5 to D10'])
    xticks([1.5,3.5])
    xticklabels({'D5','D10'})  
end

