%% Create folders for each mouse, cell with folder names of timepoints that are blinded


clear all

%path = '/Users/ktsimring/Dropbox (MIT)/Sur Lab/Projects/Development project/Binocular Matching/Spine imaging/';
path = 'C:/Users/Katya-PC/Dropbox (MIT)/Sur Lab/Projects/Development project/Binocular Matching/Spine imaging/';
input_path = 'D:BinocularMatching\Spines\2P_imaging';

tablepath = fullfile(path, 'Analyzed Data');
load(fullfile(tablepath,'all_analyzed_spines_BM014-BM021.mat')); 

savepath = fullfile(path, 'Chronic Imaging', 'FOV_alignment');

stim_binoc = all_stims(strcmp({all_stims.name},'binoc')==1);
stim_binoc = rmfield(stim_binoc,'name');
stim_table_binoc = struct2table(stim_binoc);
stim_table_binoc.mouse_cell = strcat(stim_table_binoc.all_mouse, '-', stim_table_binoc.all_cells);

%%
unique_mouse_cell = unique(stim_table_binoc.mouse_cell);
all_mice_cell = [];
all_day = [];
all_blind_day = [];
s='A':'J';
blind_key = [];
for i = 1:length(unique_mouse_cell)
    
    mouse_cell = unique_mouse_cell{i};
    split_mouse_cell = split(mouse_cell,'-');
    mouse = split_mouse_cell{1};
    cell = split_mouse_cell{2};
    vals = stim_table_binoc(contains(stim_table_binoc.mouse_cell,mouse_cell),:);
    unique_day = sort(unique(vals.all_day));
    if isdir(fullfile(input_path, mouse))
        if length(unique_day)>=2
    %         day = cellfun(@(x) split(x, '_'), unique_day, 'UniformOutput',false);
    %         day = cellfun(@(x) x{1,:}, day, 'UniformOutput',false);
    %         day = cellfun(@(x) x(2:end), day, 'UniformOutput',false);
    %         day = cellfun(@(x) str2num(x), day, 'UniformOutput',false);
            
            %blind_days = rand(length(unique_day),1);
            str=s(randi(numel(s), length(unique_day)));
            for d = 1:length(unique_day)
                dends = vals(contains(vals.all_day, unique_day{d}), 'all_fovs');
                dends = table2array(dends);
                dends = unique(dends);
                day_name = ['Day_', str(d,:)];
                all_mice_cell = [all_mice_cell; {mouse_cell}];
                all_day = [all_day; {unique_day{d}}];
                all_blind_day = [all_blind_day; {day_name}];
                for dd = 1:length(dends)
                    if contains(dends{dd}, 'Dend')
                       
                       dend = split(dends{dd}, '_');
                       dend = dend{1};
                       savefile = fullfile(savepath,mouse,cell,day_name,dend);
                       mkdir(savefile);
                       input_file  = fullfile(input_path,mouse,unique_day{d},cell,'Data',dend);
                       copyfile(fullfile(input_file, 'REFs', '/*.tif'), savefile)
                       copyfile(fullfile(input_file, 'ROIs',  '/*.zip'), savefile)
                    end
    
    
    
                 end
            end
        end
    end
end
%%
table_blind_key = table([all_mice_cell, all_day, all_blind_day])
       
writetable(table_blind_key,fullfile(savepath, "key_to_day_folders.csv"));

