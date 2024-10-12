function x = plot_errorbar_swarm(y, err,all_data, x)
    ngroups = size(y, 1);
    nbars = size(y, 2);

    for i = 1:ngroups
       for ii = 1:nbars
            swarmchart(x(i,ii)*ones(size(all_data{i,ii})), all_data{i,ii}); hold on;
            errorbar(x(i,ii), y(i,ii), err(i,ii),'Marker','_', 'Color','k', 'LineStyle','none');   
       end
    end
    hold off
end
