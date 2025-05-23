function all_counts = get_10day_soma_identity_contra_ipsi(data, unique_mice_cells_D1_D10)
%% get eye-specific identity shifts
elements = {0:1, 0:1, 0:1, 0:1}; %cell array with N vectors to combine
result = get_combin(elements);
combinations = {};
for i = 1:size(result,1)
    combinations = [ combinations; strcat(num2str(result(i,1)),num2str(result(i,2)), '_', num2str(result(i,3)), num2str(result(i,4)))];

end
all_counts = [];
for i = 1:length(data)-1
    temp = data{i};
    d1 = unique(temp.days);
    
    temp2 = data{i+1};
    d5 = unique(temp2.days);

    contra_ipsi_d1 = get_contra_ipsi_table(temp,unique_mice_cells_D1_D10);
    contra_ipsi_d5 = get_contra_ipsi_table(temp2,unique_mice_cells_D1_D10);
    
    a = strcat(contra_ipsi_d1.ci, '_',contra_ipsi_d5.ci);
    a = cellstr(a);
   
    for cc = 1:length(combinations)
        counts = length(find(ismember(a,combinations(cc))));
        combos = split(combinations(cc), '_');
        str_count = [d1{1},' ',combos{1}, ' [', num2str(counts), '] ', d5{1}, ' ', combos{2}];
        all_counts = [all_counts; {str_count}];
    end
end
end
%% filter data and get the joined tbl by contra and ipsi session
function contra_ipsi_tbl = get_contra_ipsi_table(temp, unique_mice_cells_D1_D10)
    filter_temp = [];
    for ii = 1:height(temp)
        if ismember(temp(ii,:).mouse_cell,unique_mice_cells_D1_D10)
            filter_temp = [filter_temp; temp(ii,:)];
        end
    end
    contra = filter_temp(strcmp([filter_temp.session{:}], "contra"),:);
    ipsi = filter_temp(strcmp([filter_temp.session{:}], "ipsi"),:);
    contra_ipsi_tbl = outerjoin(contra, ipsi, "Keys","mouse_cell", "MergeKeys",true);
    contra_ipsi_tbl=contra_ipsi_tbl(~any(ismissing(contra_ipsi_tbl),2),:);   
    contra_ipsi_tbl.ci = strcat(num2str(contra_ipsi_tbl.resp_contra),num2str(contra_ipsi_tbl.resp_ipsi));
end

