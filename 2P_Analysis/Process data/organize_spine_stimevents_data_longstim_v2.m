function organize_spine_stimevents_data_longstim_v2(FOV,stimname,mc_path, importpath, videopath_name,stimevent_path,savepath,stimfiles_names)
%import processed data
A = importdata(fullfile(importpath,FOV));
F = A.data';
top=1; %1 or 2 depending on data
if (A.textdata{1}) == ' '
    left= 2;
else
    left = 1;
end
F = F(left:end,top:end);
ROIdata = F(any(diff(F,[],2) ~= 0, 2),:)'; %remove repeating elements

%import order of concatenated stim sessions 
stimorder = importdata(fullfile(mc_path, 'order_of_files.txt'));
stims = [];

stimorder_name = split(stimorder{1},'_');
for kk = 1:length(stimfiles_names)
    ind = find(contains(stimorder_name,stimfiles_names{kk}));
    if ~isempty(ind)
        break;
    end
end




start_num = 1;
%load(fullfile(stimevent_path))
%vars = who();
videofiles = dir(fullfile(videopath_name));
for k = 1:length(stimorder) 
    stimorder_name_temp = split(stimorder{k},'_');
    stimorder_name_temp = stimorder_name_temp{ind};
    if contains(stimorder_name_temp, '-');
        stimorder_name_temp = split(stimorder_name_temp, '-');
        stimorder_name_temp = stimorder_name_temp{1};
    end
    if contains(stimorder_name_temp, '/');
        stimorder_name_temp = split(stimorder_name_temp, '/');
        stimorder_name_temp = stimorder_name_temp{1};
    end
    simplify_name = stimorder_name_temp;
    if exist(fullfile(stimevent_path,[simplify_name,'_stim_gratings.mat']))
        load(fullfile(stimevent_path, [simplify_name, '_stim_gratings.mat']));
    elseif exist(fullfile(stimevent_path,[simplify_name,'_stim.mat']))
        load(fullfile(stimevent_path, [simplify_name, '_stim.mat']));
    else
        stim = [];
    end
    
    videofile = videofiles(contains({videofiles.name},simplify_name));
    info = imfinfo(fullfile(videofile.folder,videofile.name, 'Tiffs','Ch2.tif'));
    numfiles = length(info);
    saveFolder = fullfile(savepath,[stimname, '_',simplify_name]);
   
    if ~isfolder(saveFolder)
        mkdir(saveFolder)
    end
    if (start_num+numfiles-1)>length(ROIdata)
        fileID = fopen(fullfile(saveFolder,'error.txt'),'w');
        error_message = "motion correct file shorter in length";
        fprintf(fileID, "%s", error_message);
        fclose(fileID);
    else

    rois = ROIdata(start_num:start_num+numfiles-1,:);
    start_num = start_num + numfiles;
    saveStimFolder = fullfile(savepath,[stimname,'_', simplify_name], 'stim');
    
    
    if ~isfolder(saveStimFolder)
        mkdir(saveStimFolder)
    end
    
    copyfile(fullfile(videopath_name, videofile.name, '/*.env'), saveFolder)
    copyfile(fullfile(videopath_name, videofile.name,  '/*.csv'), saveFolder)
    copyfile(fullfile(videopath_name,  videofile.name, '/*.xml'), saveFolder)
    save(fullfile(saveStimFolder, 'stim'), 'stim');
    save(fullfile(saveFolder,['ROIdata.mat']),'rois'); 
end
end


