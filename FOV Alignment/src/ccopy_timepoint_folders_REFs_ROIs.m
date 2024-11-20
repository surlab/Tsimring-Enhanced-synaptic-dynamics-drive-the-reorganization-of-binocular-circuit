%% Copy Tif files and ROISet.zip into blinded or non blinded folders 
% Only do this for updated files

clear all

%path = '/Users/ktsimring/Dropbox (MIT)/Sur Lab/Projects/Development project/Binocular Matching/Spine imaging/';
path = 'G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/';
%input_path = 'D:BinocularMatching\Spines\2P_imaging';
input_path = 'F:\Binocular_matching\Spine_imaging\2P_imaging';
%input_path = '/Volumes/Katya6/BinocularMatching/Spines/2P_imaging';

savepath = fullfile(path, 'Chronic Imaging', 'FOV_alignment');


%% Copy from dendritic distance file to hard drive
mice = dir(savepath);
mice = mice(~ismember({mice.name},{'.','..', '.DS_Store', 'ReadMe.txt'}));
mice = mice(~contains({mice.name}, 'annotate'));

for m = 1:length(mice)
    cells = dir(fullfile(mice(m).folder, mice(m).name));
    cells = cells(~ismember({cells.name},{'.','..', '.DS_Store', 'ReadMe.txt'}));
    T = table;
    blind_folder = fullfile(distance_path, mouse, 'key_to_day_folders.csv');
    if exist(blind_folder,"file")
         is_blinded = 1;
         tbl = readtable(blind_folder);
         tbl_filt = tbl(ismember(tbl.Var1_3, cells),{'Var1_5';'Var1_4'});
         real_date = tbl_filt.Var1_4;
         blind_date = tbl_filt.Var1_5;
             
         ind1 = find(contains(real_date, day));
         day = blind_date{ind1};
           
     end
    for c = 1:length(cells)
        days = dir(fullfile(cells(c).folder, cells(c).name));
        days = days(~ismember({days.name},{'.','..', '.DS_Store'}));
        days = days(~contains({days.name}, 'Compare'));

        x_t = [];
        cell_segments = [];
        [~,inds] = sort({days.name});
        days = days(inds);
        for d = 1:length(days)
            dends = dir(fullfile(days(d).folder, days(d).name));
            dends = dends(~ismember({dends.name},{'.','..', '.DS_Store'}));
            for dd = 1:length(dends)
                input_file  = fullfile(input_path,mice(m).name,days(d).name,cells(c).name,'Data',dends(dd).name, 'ROIs');
                copyfile_path = fullfile(dends(dd).folder, dends(dd).name, '*.zip');     
                copyfile(copyfile_path, input_file);
            end
                   
    


        end
    end
end
%% Copy from hard drive to dendritic distance file 
mice = dir(input_path);
mice = mice(~ismember({mice.name},{'.','..', '.DS_Store'}));

for m = 1:length(mice)
    mouse_name = mice(m).name;
    disp(["Running mouse: ", mouse_name]);
       
    days = dir(fullfile(mice(m).folder, mice(m).name));
    days =days(~ismember({days.name},{'.','..', '.DS_Store'}));

    for d = 1:length(days)
        days_name = days(d).name;
        disp(["Running day: ", days_name]);
        
        cells = dir(fullfile(days(d).folder, days(d).name));
        cells = cells(~ismember({cells.name},{'test_vis','.','..', '.DS_Store'}));
        for c = 1:length(cells)
            cell_name = cells(c).name;
            disp(["Running cell: ", cell_name]);
            %check if dendritic distance file is blinded
            blind_folder = fullfile(savepath, mouse_name, 'key_to_day_folders.csv');

            if exist(blind_folder,"file")
                 tbl = readtable(blind_folder);
                 tbl_filt = tbl(ismember(tbl.Var1_3, cell_name),{'Var1_5';'Var1_4'});
                 real_date = tbl_filt.Var1_4;
                 blind_date = tbl_filt.Var1_5;
                     
                 ind1 = find(contains(real_date, days_name));
                 if isempty(ind1)
                     continue;
                 else
                    day = blind_date{ind1};
                 end
               
            else
                day = days_name;
            end
            dends = dir(fullfile(cells(c).folder, cells(c).name, 'Data'));
            dends = dends(~ismember({dends.name},{'Soma','.','..', '.DS_Store'}));

            for dd = 1:length(dends)
                disp(["Running dendrite: ",dends(dd).name]);
                input_file_rois_dir  = dir(fullfile(dends(dd).folder,dends(dd).name, 'ROIs'));
                input_file_rois_dir = input_file_rois_dir(~ismember({input_file_rois_dir.name},{'.','..', '.DS_Store'}));
                if ~isempty(input_file_rois_dir)
                %copy dend.zip path
                input_file_dend_dir = input_file_rois_dir(contains({input_file_rois_dir.name},{'dend'})&~contains({input_file_rois_dir.name},{'._R'}));
                file_data = split(input_file_dend_dir.date, ' ');
                file_data = file_data{1};
                save_file = fullfile(savepath,mouse_name,cell_name,day, dends(dd).name,input_file_dend_dir.name);
                copyfile_path = fullfile(input_file_dend_dir.folder,input_file_dend_dir.name);
                if ~exist(save_file, 'file')
                        mkdir(fullfile(savepath,mouse_name,cell_name,day, dends(dd).name))
                        copyfile(copyfile_path, save_file);
                
                else
                        dendritic_file = dir(save_file);
                        dend_date = split(dendritic_file.date);
                        dend_date = dend_date{1};
                        if ~strcmp(dend_date,file_data)
                            copyfile(copyfile_path, save_file);
                        end
                end
               
                
                %copy roi.zip path
                input_file_rois_dir = input_file_rois_dir(~contains({input_file_rois_dir.name},{'dend'}));
                for r = 1:length(input_file_rois_dir)
                    filename = input_file_rois_dir(r).name;
                    file_data = split(input_file_rois_dir(r).date, ' ');
                    file_data = file_data{1};
                    save_file = fullfile(savepath,mouse_name,cell_name,day, dends(dd).name, filename);
                    copyfile_path = fullfile(input_file_rois_dir(r).folder,input_file_rois_dir(r).name);
                    if ~exist(save_file, 'file')
                        mkdir(fullfile(savepath,mouse_name,cell_name,day, dends(dd).name))
                        copyfile(copyfile_path, save_file);
                
                    else
                        dendritic_file = dir(save_file);
                        dend_date = split(dendritic_file.date);
                        dend_date = dend_date{1};
                        if ~strcmp(dend_date,file_data)
                            copyfile(copyfile_path, save_file);
                        end
                    end
                end

                %copy .tif path
                input_file_ref_dir  = dir(fullfile(dends(dd).folder,dends(dd).name, 'REFs'));
                input_file_ref_dir = input_file_ref_dir(~ismember({input_file_ref_dir.name},{'.','..', '.DS_Store'}));
                
                save_file = fullfile(savepath,mouse_name,cell_name,day, dends(dd).name, [dends(dd).name,'.tif']);
                dendritic_file = dir(save_file);
                file_data = split(input_file_ref_dir.date, ' ');
                file_data = file_data{1};
                copyfile_path = fullfile(input_file_ref_dir.folder,input_file_ref_dir.name);
                if ~exist(save_file, 'file')
                        mkdir(fullfile(savepath,mouse_name,cell_name,day, dends(dd).name))
                        copyfile(copyfile_path, save_file);
                
                else
                        dendritic_file = dir(save_file);
                        dend_date = split(dendritic_file.date);
                        dend_date = dend_date{1};
                        if ~strcmp(dend_date,file_data)
                            copyfile(copyfile_path, save_file);
                        end
                end
                end
            end
                   
    


        end
    end
end
