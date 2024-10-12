function tbl = get_10day_spine_fates(D1_D5_table, D5_D10_table)
vars = {'D5','session', 'all_mice_cells', 'all_fovs', 'structure_type', 'dendritic_type'};
D1_D5 = D1_D5_table(:,vars);
D1_D5.fovs_mouse_cell = strcat(D1_D5.all_mice_cells, "_", D1_D5.all_fovs);
D1_D5.D5_roi = num2str(D1_D5.D5.all_roi_inds);
D1_D5.fovs_mouse_cell_roi = strcat(D1_D5.all_mice_cells, "_", D1_D5.all_fovs,"_", D1_D5.D5_roi);
mice_cells_D1_D5 = unique(D1_D5.fovs_mouse_cell);

D5_D10 = D5_D10_table(:,vars);
D5_D10.fovs_mouse_cell = strcat(D5_D10.all_mice_cells, "_", D5_D10.all_fovs);
D5_D10.D5_roi = num2str(D5_D10.D5.all_roi_inds);
D5_D10.fovs_mouse_cell_roi = strcat(D5_D10.all_mice_cells, "_", D5_D10.all_fovs,"_", D5_D10.D5_roi);

mice_cells_D5_D10 = unique(D5_D10.fovs_mouse_cell);
mice_cells_D1_D10 = intersect(mice_cells_D5_D10,mice_cells_D1_D5);

temp = [];
for i = 1:height(D1_D5)
    if ismember(D1_D5(i,:).fovs_mouse_cell,mice_cells_D1_D10)
        temp = [temp; D1_D5(i,:)];
    end
end
D1_D5 = temp(contains(temp.session,"binoc"),{'structure_type', 'fovs_mouse_cell_roi', 'D5_roi', 'dendritic_type'});

temp = [];
for i = 1:height(D5_D10)
    if ismember(D5_D10(i,:).fovs_mouse_cell,mice_cells_D1_D10)
        temp = [temp; D5_D10(i,:)];
    end
end
D5_D10 = temp(contains(temp.session,"binoc"),{'structure_type', 'fovs_mouse_cell_roi','D5_roi', 'dendritic_type'});
%% Find spines that were lost on D5 and spines that were added on D10 and remove from table
num_lost = length(find(contains(D1_D5.structure_type, "lost")));
D1_D5 = D1_D5(~contains(D1_D5.structure_type, "lost"),:);
num_formed = length(find(contains(D5_D10.structure_type, "formed")));
D5_D10 = D5_D10(~contains(D5_D10.structure_type, "formed"),:);

%% count retained, retained then lost, formed then lost, formed then retained
D1_D10 = outerjoin(D1_D5, D5_D10, "Keys","fovs_mouse_cell_roi", "MergeKeys",true);
D1_D10=D1_D10(~any(ismissing(D1_D10),2),:);

temp = D1_D10;
type = '';
temp.structure_type_D1_D10 = strcat(temp.(['structure_type_D1_D5',type]), '_', temp.(['structure_type_D5_D10',type]));
combinations = ["retained_retained"; "retained_lost"; "formed_lost"; "formed_retained"]; 
a = temp.structure_type_D1_D10;
a = cellstr(a);
counts = [];
for cc = 1:length(combinations)
    index = find(ismember(a,combinations(cc)));
    counts = [counts, length(index)];
end

%% Create table for sankey 
string_D1_D5_lost = ['D1 [', num2str(num_lost), '] ', 'D5 lost'];
string_D5_D10_added = ['D5 [', num2str(num_lost), '] ', 'D10 added'];
string_D5_retained_D10_retained = ['D5 retained [', num2str(counts(1)), '] ', 'D10 retained'];
string_D5_retained_D10_lost = ['D5 retained [', num2str(counts(2)), '] ', 'D10 lost'];
string_D5_added_D10_lost = ['D5 added [', num2str(counts(3)), '] ', 'D10 lost'];
string_D5_added_D10_retained = ['D5 added [', num2str(counts(4)), '] ', 'D10 retained'];

tbl = {string_D1_D5_lost;string_D5_D10_added;string_D5_retained_D10_retained;...
    string_D5_retained_D10_lost;string_D5_added_D10_lost;string_D5_added_D10_retained};

end
