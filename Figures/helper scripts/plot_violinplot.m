function x = plot_violinplot(all_data, x, colors)
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
    v = violinplot(all_data_temp,all_x);
    
    for i = 1:length(v)
        v(i).ViolinColor = colors(i,:);
        v(i).ViolinAlpha = 0.8;
        v(i).MedianPlot.SizeData = 20;
        v(i).BoxWidth = 0.08;
        v(i).ScatterPlot.Marker = 'o';
        v(i).ScatterPlot.SizeData = 0.2;
        v(i).ScatterPlot.MarkerFaceColor = 'k';
        v(i).ScatterPlot.MarkerFaceAlpha = 0.3;
    end
    hold off
end
