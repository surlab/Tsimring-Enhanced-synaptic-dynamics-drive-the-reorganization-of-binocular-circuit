function all_counts = get_10day_retained_spine_soma_identity(retain_D1_D5, retain_D5_D10,mice_cells_D1_D10,resp)
%% get 3 time point data
temp = [];
for i = 1:height(retain_D1_D5)
    if ismember(retain_D1_D5(i,:).fovs_mouse_cell,mice_cells_D1_D10)
        temp = [temp; retain_D1_D5(i,:)];
    end
end
retain_D1_D5 = temp;

temp = [];
for i = 1:height(retain_D5_D10)
    if ismember(retain_D5_D10(i,:).fovs_mouse_cell,mice_cells_D1_D10)
        temp = [temp; retain_D5_D10(i,:)];
    end
end
retain_D5_D10 = temp;
%% get 3 timepoint spine tracking
D1_D10 = outerjoin(retain_D1_D5, retain_D5_D10, "Keys","fovs_mouse_cell_roi", "MergeKeys",true);
D1_D10=D1_D10(~any(ismissing(D1_D10),2),:);

%% retain spines
retain_D1_D5 = D1_D10(:, {'D1', 'D5_retain_D1_D5', 'session_retain_D1_D5', 'all_fovs_retain_D1_D5', 'all_mice_cells_retain_D1_D5'});
retain_D5_D10 = D1_D10(:, {'D10', 'D5_retain_D5_D10', 'session_retain_D5_D10', 'all_fovs_retain_D5_D10', 'all_mice_cells_retain_D5_D10'});
retain_D1_D5 = renamevars(retain_D1_D5, ["session_retain_D1_D5",  "D5_retain_D1_D5", "all_fovs_retain_D1_D5", "all_mice_cells_retain_D1_D5"], ...
    ["session", "D5", "all_fovs", "all_mice_cells"]);
retain_D5_D10 = renamevars(retain_D5_D10, ["session_retain_D5_D10",  "D5_retain_D5_D10", "all_fovs_retain_D5_D10", "all_mice_cells_retain_D5_D10"],...
    ["session", "D5", "all_fovs", "all_mice_cells"]); 
%% get eye-specific identity shifts
all_temps = {retain_D1_D5, retain_D5_D10};
all_counts = [];
days = {'D1', 'D5', 'D10'};
elements = {0:1, 0:1}; %cell array with N vectors to combine
result = get_combin(elements);
combinations = {};
for i = 1:size(result,1)
    combinations = [ combinations; strcat(num2str(result(i,1)),'_',num2str(result(i,2)))];
end
for i = 1:length(all_temps)
    temp = all_temps{i};
    d1 = days{i};
    d5 = days{i+1};
    temp = splitvars(temp);
    d1_var = strcat(d1, '_', 'resp');
    d5_var = strcat(d5, '_', 'resp');
    d1_soma_var = strcat(d1, '_', 'soma_resp');
    d5_soma_var = strcat(d5, '_', 'soma_resp');
    d1_roi = strcat(d1, '_', 'all_roi_inds');
    d5_roi = strcat(d5, '_', 'all_roi_inds');
    temp.roi_fovs_mouse_cell = strcat(num2str(temp.(d1_roi)),'_',num2str(temp.(d5_roi)), '_', temp.all_fovs, '_', temp.all_mice_cells);

    varnames = temp.Properties.VariableNames;
    filter = {'roi_fovs_mouse_cell', 'session'};
    var_filter = setxor(varnames, filter);
    
    %extract each session and then combine
    contra = temp(contains(temp.session, 'contra'),:);
    contra = contra(:, {d1_var, d5_var,d1_soma_var, d5_soma_var,'roi_fovs_mouse_cell'});
    ipsi = temp(contains(temp.session, 'ipsi'),:);
    ipsi = ipsi(:, {d1_var, d5_var,d1_soma_var, d5_soma_var,'roi_fovs_mouse_cell'});
    contra_ipsi = outerjoin(contra, ipsi, "Keys","roi_fovs_mouse_cell", "MergeKeys",true);
    
    %remove rows with NaN values 
    contra_ipsi=contra_ipsi(~any(ismissing(contra_ipsi),2),:);

    d1_contra = strcat(d1, '_resp_contra');
    d1_ipsi = strcat(d1, '_resp_ipsi');
    d1_contra_soma = strcat(d1, '_soma_resp_contra');
    d1_ipsi_soma = strcat(d1, '_soma_resp_ipsi');
    d1_spine = arrayfun(@(x,y) strcat(x,y),num2str(contra_ipsi.(d1_contra)),num2str(contra_ipsi.(d1_ipsi)), 'UniformOutput',false);
    d1_soma = arrayfun(@(x,y) strcat(x,y),num2str(contra_ipsi.(d1_contra_soma)),num2str(contra_ipsi.(d1_ipsi_soma)), 'UniformOutput',false);
    
    d5_contra = strcat(d5, '_resp_contra');
    d5_ipsi = strcat(d5, '_resp_ipsi');
    d5_contra_soma = strcat(d5, '_soma_resp_contra');
    d5_ipsi_soma = strcat(d5, '_soma_resp_ipsi');
    d5_spine = arrayfun(@(x,y) strcat(x,y),num2str(contra_ipsi.(d5_contra)),num2str(contra_ipsi.(d5_ipsi)), 'UniformOutput',false);
    d5_soma = arrayfun(@(x,y) strcat(x,y),num2str(contra_ipsi.(d5_contra_soma)),num2str(contra_ipsi.(d5_ipsi_soma)), 'UniformOutput',false);
       
    if resp
        inds = ~strcmp(d1_soma, '00')&~strcmp(d5_soma, '00');
        contra_ipsi = contra_ipsi(inds,:);
        d1_soma = d1_soma(inds);
        d5_soma = d5_soma(inds);
        d1_spine = d1_spine(inds);
        d5_spine = d5_spine(inds);
    end
    d1_match = cellfun(@(x,y) num2str(strcmp(x,y)), d1_spine, d1_soma);
    d5_match = cellfun(@(x,y) num2str(strcmp(x,y)), d5_spine, d5_soma);
    contra_ipsi.d1_d5_match = strcat(d1_match, '_', d5_match);
     
    %label each spine by eye-specific identity for D1 and D5 and count
    %combinations
    a = contra_ipsi.d1_d5_match;
    a = cellstr(a);
   
    for cc = 1:length(combinations)
        counts = length(find(ismember(a,combinations(cc))));
        combos = split(combinations(cc), '_');
        str_count = [d1,' ',combos{1}, ' [', num2str(counts), '] ', d5, ' ', combos{2}];
        all_counts = [all_counts; {str_count}];
    end
end
end


