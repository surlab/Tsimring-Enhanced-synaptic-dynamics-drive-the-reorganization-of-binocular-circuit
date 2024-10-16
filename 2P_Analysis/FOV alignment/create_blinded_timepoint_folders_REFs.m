%% Create folders for each mouse, cell with folder names of timepoints that are blinded


clear all
%path = 'C:/Users/Katya-PC/Dropbox (MIT)/Sur Lab/Projects/Development project/Binocular Matching/Spine imaging/';
path = 'G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';
path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';
input_path = '/Volumes/Katya5/BinocularMatching/Spines/2P_imaging';
input_path = fullfile(input_path, "Cells");
savepath = fullfile(path, 'Chronic Imaging', 'FOV_alignment');


%%


s='A':'J';
blind_key = [];
mice = {'BM017'};
for i = 1:length(mice)
    all_mice_cell = [];
    all_day = [];
    all_blind_day = [];
    mouse = mice{i};
    days = dir(fullfile(input_path, mouse));
    days = days(~contains({days.name}, '.'));
    str=s(randi(numel(s), length(days)));
    for ii = 1:length(days)
        day_name = ['Day_' str(ii,1:3)];
        day = days(ii).name;
        %cells = dir(fullfile(input_path, mouse, day, 'Data'));
        cells = dir(fullfile(input_path, mouse, day));
        cells = cells(~contains({cells.name}, '.'));
        for iii = 1:length(cells)
                cell1 = cells(iii).name;
                all_mice_cell = [all_mice_cell; {mouse, '_', cell1}];
                all_day = [all_day; {day}];
                all_blind_day = [all_blind_day; {day_name}];
                %dends = dir(fullfile(input_path, mouse, day, 'Data', cell1));
                dends = dir(fullfile(input_path, mouse, day, cell1,'Data'));
                dends = dends(~contains({dends.name}, '.'));
                
                for dd = 1:length(dends)
                    
                    dend = dends(dd).name;
                    if contains(dend, 'Dend')
                       
                       
                       savefile = fullfile(savepath,mouse,cell1,day_name,dend);
                      
                       %input_file  = fullfile(input_path,mouse,day,'Data',cell1,dend);
                       input_file  = fullfile(input_path,mouse,day,cell1,'Data',dend);
                       reffile_name = fullfile(input_file, 'REFs', '/*.tif');
                       file = dir(reffile_name);
                       if ~isempty(file)
                            mkdir(savefile);
                            copyfile(reffile_name, savefile)
                            roifile_name = fullfile(input_file, 'ROIs',  '/*.zip');
                            if ~isempty(dir(roifile_name))
                                copyfile(roifile_name, savefile);
                            end
                       end

                    end
    
    
    
                 end
            end
    end
    table_blind_key = table([all_mice_cell, all_day, all_blind_day])
       
    writetable(table_blind_key,fullfile(savepath,mice{i}, "key_to_day_folders.csv"));

end


%% Copy any additional REF files to blinded files
tbl = readtable(fullfile(savepath,mice, "key_to_day_folders.csv"));
arr = table2array(tbl);
num_folders = size(arr,1);
blinded_ind = 5; %blinded folder index in array
real_ind = 4; %real folder index in array
cell_ind = 3; %cell index in array
mice = {'BM023'};
for m = 1:length(mice)
    mouse = mice{m};
    for i = 1:num_folders
        day = arr{i,real_ind};
        blind_day = arr{i,blinded_ind};
        cell1 = arr{i,cell_ind};
        dends = dir(fullfile(input_path, mouse, day, 'Data', cell1));
        dends = dends(~contains({dends.name}, '.'));
        for dd = 1:length(dends)
            dend = dends(dd).name;
            if contains(dend, 'Dend')
                savefile = fullfile(savepath,mouse,cell1,blind_day,dend);
                input_file  = fullfile(input_path,mouse,day,'Data',cell1,dend);
                file = dir(fullfile(input_file, 'REFs', '/*.tif'));
                       if ~isempty(file)
                            mkdir(savefile);
                        copyfile(fullfile(input_file, 'REFs', '/*.tif'), savefile)
                       end
                       %copyfile(fullfile(input_file, 'ROIs',  '/*.zip'), savefile)
            end
    
    
    
         end
    end
end
%% Copy ROI files from blinded folders to hard drive 
savepath = 'G:\My Drive\Sur Lab\Development project\Binocular_Matching\Spine_imaging\Chronic Imaging\FOV_alignment\';
savepath = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/FOV_alignment/';
input_path = 'D:BinocularMatching\Spines\2P_imaging';
input_path = '/Volumes/Katya5/BinocularMatching/Spines/2P_imaging';

input_path = fullfile(input_path, "Cells");



blinded_ind = 5; %blinded folder index in array
real_ind = 4; %real folder index in array
cell_ind = 3; %cell index in array
mice = {'BM016';'BM017'};
for m = 1:length(mice)
    tbl = readtable(fullfile(savepath,mice{m}, "key_to_day_folders.csv"));
    arr = table2array(tbl);
    num_folders = size(arr,1);
    mouse = mice{m};
    for i = 1:num_folders
        day = arr{i,real_ind};
        blind_day = arr{i,blinded_ind};
        cell1 = arr{i,cell_ind};
        dends = dir(fullfile(input_path, mouse, day, cell1,'Data'));
        dends = dends(~contains({dends.name}, '.'));
        for dd = 1:length(dends)
            dend = dends(dd).name;
            if contains(dend, 'Dend')
                file = fullfile(savepath,mouse,cell1,blind_day,dend);
                input_file  = fullfile(input_path,mouse,day,cell1,'Data',dend);
                savefile = fullfile(input_file, 'ROIs');
                       if ~isempty(dir(fullfile(file,  '/*.zip')))
                            mkdir(savefile);
                       
                        
                        copyfile(fullfile(file,  '/*.zip'), savefile)
                       end
                       %copyfile(fullfile(input_file, 'ROIs',  '/*.zip'), savefile)
            end
    
    
    
         end
    end
end