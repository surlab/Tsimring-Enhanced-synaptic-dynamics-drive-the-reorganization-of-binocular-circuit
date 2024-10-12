%% For Panel D to F, plot fraction of responsive spines by neuron
% color code by whether soma was responsive or not
close all
unique_mouse = unique(contra_ipsi_binoc.mouse_cell);
session = { '','_contra', '_ipsi'};%empty string is the binoc session
unique_days = {'D1', 'D5', 'D10'};
colors = {'g','r', 'b' };
all_soma_resp_by_session = [];
all_session = [];
for i = 1:length(session)
    all_resp = 0;
    resp_type = ['resp', session{i}];
    
    %Plot fraction of responsive spines on somas pooled by day
    [fraction_spine_resp_by_neuron_by_day,soma_resp_by_day] = get_eye_resp_neuron_day(unique_mouse,unique_days, contra_ipsi_binoc, all_resp,resp_type);
    figure("Position", [283,132,150,344])
    subplot(2,1,1)
    plot_resp_by_mouse(fraction_spine_resp_by_neuron_by_day,soma_resp_by_day,colors{i})
    xticks([1:3])
    xticklabels(unique_days)
    xlim([0,4])
    ylabel('% responsive spines')
    ylim([0,100])
    
    subplot(2,1,2)  
    %Plot fraction of responsive spines on unresponsive vs responsive somas pooled across days (inset)
    disp(session{i})
    soma_resp_by_session = compare_resp_unresp_soma_spine_fractions(fraction_spine_resp_by_neuron_by_day, soma_resp_by_day);
    all_soma_resp_by_session = [all_soma_resp_by_session; soma_resp_by_session];
    
    all_session = [all_session; repmat(colors(i),(size(soma_resp_by_session)))];


    ylabel('% responsive spines')
    ylim([0,100])
end

%%
function plot_resp_by_mouse(fraction_spine_resp_by_neuron_by_day,soma_resp_by_day,c)
mean_fract = squeeze(nanmean(fraction_spine_resp_by_neuron_by_day,1));
num_cells_by_day = sum(~isnan(fraction_spine_resp_by_neuron_by_day));
sem_fract = squeeze(nanstd(fraction_spine_resp_by_neuron_by_day,1))./sqrt(num_cells_by_day);
for i = 1:size(fraction_spine_resp_by_neuron_by_day,1)
    plot([1:3], fraction_spine_resp_by_neuron_by_day(i,:),'Color', [0.5,0.5,0.5],'LineWidth',0.5); hold on;

    for ii = 1:size(fraction_spine_resp_by_neuron_by_day,2)
         if soma_resp_by_day(i,ii)==1
            scatter(ii, fraction_spine_resp_by_neuron_by_day(i,ii), 100,'.',c, 'MarkerFaceAlpha',0.6); hold on;
         else
            scatter(ii, fraction_spine_resp_by_neuron_by_day(i,ii), 100,'.','k','MarkerFaceAlpha',0.6); hold on;
         end
    end
    %errorbar([1:size(mean_fract,2)], mean_fract, sem_fract, '-o', 'MarkerFaceColor', 'auto', 'MarkerSize',5, 'LineWidth',2); hold on;
end

end

%%
function [fraction_resp_by_mouse_by_day,soma_resp_by_mouse_by_day] = get_eye_resp_neuron_day(unique_mouse,unique_days, contra_ipsi_binoc, all_resp,resp_type)
    fraction_resp_by_mouse_by_day = NaN(length(unique_mouse), length(unique_days));
    soma_resp_by_mouse_by_day = NaN(length(unique_mouse), length(unique_days));
    for i = 1:length(unique_mouse)
         for iii = 1:length(unique_days)
             temp = contra_ipsi_binoc(strcmp(contra_ipsi_binoc.mouse_cell, unique_mouse{i})&...
                 strcmp(contra_ipsi_binoc.days, unique_days{iii}),:);
             if ~isempty(temp)
                 if all_resp
                     resp = height(temp(~strcmp(temp.cib_spine, resp_type),:));
                     soma_resp = max(~strcmp(temp.cib_soma, resp_type));
                 else
                     resp = height(temp(temp.(resp_type)==1,:));
                     soma_resp = max(temp.(['soma_',resp_type]));
                 end
                 fraction_resp_by_mouse_by_day(i,iii) = resp*100/height(temp);
                 soma_resp_by_mouse_by_day(i,iii) = soma_resp;
             end
         end
    end
end

%%
function [soma_resp_fract_spines] =  compare_resp_unresp_soma_spine_fractions(fraction_resp_by_mouse_by_day, soma_resp_by_mouse_by_day)
    soma_resp_fract_spines = fraction_resp_by_mouse_by_day(soma_resp_by_mouse_by_day==1);
    soma_unresp_fract_spines = fraction_resp_by_mouse_by_day(soma_resp_by_mouse_by_day==0);
    mean_soma_resp_fract_spines = mean(soma_resp_fract_spines);
    mean_soma_unresp_fract_spines = mean(soma_unresp_fract_spines);
    sem_soma_resp_fract_spines = std(soma_resp_fract_spines)/sqrt(length(soma_resp_fract_spines));
    sem_soma_unresp_fract_spines = std(soma_resp_fract_spines)/sqrt(length(soma_unresp_fract_spines));
    
    bar([mean_soma_unresp_fract_spines,mean_soma_resp_fract_spines]); hold on;
    errorbar([mean_soma_unresp_fract_spines,mean_soma_resp_fract_spines],[sem_soma_unresp_fract_spines, sem_soma_resp_fract_spines], '.k');
    scatter(ones(size(soma_unresp_fract_spines)), soma_unresp_fract_spines, '.k','jitter', 'on', 'jitterAmount', 0.05);
    scatter(2*ones(size(soma_resp_fract_spines)), soma_resp_fract_spines,'.k', 'jitter', 'on', 'jitterAmount', 0.05);
    [p,a] = ranksum(soma_resp_fract_spines, soma_unresp_fract_spines);
    plot_sig(1.5,100,p)
    disp(['mean unresp', num2str(mean_soma_unresp_fract_spines)]);
    disp(['sem unresp', num2str(sem_soma_unresp_fract_spines)]);
    disp(['mean resp', num2str(mean_soma_resp_fract_spines)]);
    disp(['sem resp', num2str(sem_soma_resp_fract_spines)]);


end


     

