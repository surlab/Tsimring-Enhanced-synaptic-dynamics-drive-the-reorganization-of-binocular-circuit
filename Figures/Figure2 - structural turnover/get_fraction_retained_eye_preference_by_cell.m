function get_fraction_retained_eye_preference_by_cell(all_temps)

for i = 1:length(all_temps)
    temp = all_temps{i};
    fraction_contra_retained = [];
    fraction_ipsi_retained = [];
    fraction_contra_ipsi_retained = [];
    fraction_unresp_retained = [];
    fraction_binoc_retained = [];
    if i == 1
        D1 = 'D1';
        D5 = 'D5';
    else
        D1 = 'D5';
        D5 = 'D10';
    end

    lost = splitvars(temp(contains(temp.structure_type, 'lost'), {D1, 'all_mice_cells', 'all_fovs', 'session', 'structure_type'}));
    retained = splitvars(temp(contains(temp.structure_type, 'retained'),{D1, 'all_mice_cells', 'all_fovs', 'session','structure_type'}));
    
    temp = [lost;retained];
    temp.roi_fovs_mouse_cell = strcat(num2str(temp.all_roi_inds),'_', temp.all_fovs, '_', temp.all_mice_cells, '_', temp.structure_type);
    

    features = {'resp','roi_fovs_mouse_cell', 'structure_type', 'all_mice_cells'};
    join_key = 'roi_fovs_mouse_cell';
 
    contra_ipsi_binoc = concate_contra_ipsi_binoc(temp, features, join_key);
    
    contra_ipsi_binoc.ci = arrayfun(@(x,y) [num2str(x);num2str(y)], contra_ipsi_binoc.('resp_contra'),...
        contra_ipsi_binoc.('resp_ipsi'), 'UniformOutput',false );
    contra_ipsi_binoc.cib = arrayfun(@(x,y,z) [num2str(x);num2str(y);num2str(z)], contra_ipsi_binoc.('resp_contra'),...
        contra_ipsi_binoc.('resp_ipsi'),contra_ipsi_binoc.('resp'), 'UniformOutput',false );

    unique_mouse_cell = unique(contra_ipsi_binoc.all_mice_cells);


    for ii = 1:length(unique_mouse_cell)
        spines = contra_ipsi_binoc(ismember(contra_ipsi_binoc.all_mice_cells,unique_mouse_cell{ii}),:);
        
        contra_spines = spines(strcmp(spines.ci, "10"),:);
        ipsi_spines = spines(strcmp(spines.ci, "01"),:);
        contra_ipsi_spines = spines(strcmp(spines.ci, "11"),:);
        unresp_spines = spines(strcmp(spines.cib, "000"),:);
        binoc_spines = spines(spines.resp==1,:);
        
       
        fraction_contra_retained = [fraction_contra_retained, height(contra_spines(strcmp(contra_spines.structure_type, 'retained'),:))/height(contra_spines)];
        fraction_ipsi_retained = [fraction_ipsi_retained, height(ipsi_spines(strcmp(ipsi_spines.structure_type, 'retained'),:))/height(ipsi_spines)];
        fraction_contra_ipsi_retained = [fraction_contra_ipsi_retained, height(contra_ipsi_spines(strcmp(contra_ipsi_spines.structure_type, 'retained'),:))/height(contra_ipsi_spines)];
        fraction_binoc_retained = [fraction_binoc_retained, height(binoc_spines(strcmp(binoc_spines.structure_type, 'retained'),:))/height(binoc_spines)]
        fraction_unresp_retained = [fraction_unresp_retained, height(unresp_spines(strcmp(unresp_spines.structure_type, 'retained'),:))/height(unresp_spines)];  
    end
   
    g = figure("Position",[369,69,413,489]);
    set(gca, 'FontSize',8, 'FontName', 'Arial', 'Position',[268,346,377,250])
    h = figure;
    data = [fraction_binoc_retained;fraction_contra_retained;fraction_ipsi_retained;fraction_contra_ipsi_retained;fraction_unresp_retained];
  
 
    mean_data = nanmean(data, 2);
    err_data = nanstd(data')/sqrt(size(data,2));
    figure(g)
    subplot(2,2,1)
    bar(mean_data); hold on;
    errorbar(mean_data, err_data, '.k');
    for sc = 1:size(data,1)
        scatter(sc*ones(size(data(sc,:))), data(sc,:), 4, 'k', 'jitter', 'on', 'jitterAmount',0.05);
    end
    
    ylabel('Fraction')
    xticklabels({'B','C', 'I', 'C+I', 'U'})
    ylim([0,1])

    temp = data;
    temp = temp(:)';
    num_cells = size(data,2);
    eye_pref = repmat(["binoc","contra", "ipsi", "contra_ipsi", "unresp"], 1,num_cells);
    %eye_pref = repmat(["contra", "ipsi", "contra_ipsi"], 1,num_cells);

    ind = ~isnan(temp);
    temp = temp(ind);
    eye_pref = eye_pref(ind);
    figure(h)
    [p,tbl,stats] = kruskalwallis(temp,eye_pref);
    figure(g)
    plot_sig(2.5,1,p)
    figure(h)
    [results,~,~,gnames] = multcompare(stats);
    tbl = array2table(results,"VariableNames", ...
        ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
    tbl.("Group A")=gnames(tbl.("Group A"));
    tbl.("Group B")=gnames(tbl.("Group B"))
end