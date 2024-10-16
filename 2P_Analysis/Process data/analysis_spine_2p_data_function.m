function analysis_spine_2p_data_function(input_path, stim_path, path, mouse_files)
% Both Kyle and Katya worked on this script 

% program is designed to import 2P data. Select data, and then select
% timestamps. If anaylzing multiple recordings in a
% cocatenated stack, enter number of frames for each recording in the order
% in which they are in the stack as a vector. 

%scans: vector of frames in each scan in the order they are cocatenated

%% save processed data into Dropbox folder 
output_path = fullfile(path, 'Processed Data');
stimfiles_names = {'binoc';'contra'; 'ipsi' };
% disp('Saving recordings and stim data into folders');
 for m = 1:length(mouse_files)
   disp(['Running mouse: ', mouse_files{m}])
   folders = dir(fullfile(input_path,mouse_files{m}));
   Days = {folders.name};
   Days=Days(~contains({folders.name}, '.'));
   for i = 2:length(Days)
        Day = Days{i};
        disp(['Running day: ',Day])
        Cells = dir(fullfile(input_path,mouse_files{m},Day));
        Cells=Cells(~contains({Cells.name}, '.'));
        Cells=Cells(~contains({Cells.name}, 'X_'));
        for cc = 1:length(Cells)
            disp(['Running cell: ',Cells(cc).name])
            %check if folder contains "test_vis" folder for cases when I
            %test visually responsive neurons
            if contains(Cells(cc).name, "test_vis")
                test_stims = dir(fullfile(Cells(cc).folder, "test_vis"));
                test_stims=test_stims(~contains({test_stims.name}, '.'));
                for t = 1:length(test_stims)
                    if isfolder(fullfile(test_stims(t).folder,test_stims(t).name, "Results"))
                        results_path = fullfile(test_stims(t).folder,test_stims(t).name, 'Results');
                        video_path = fullfile(test_stims(t).folder,test_stims(t).name, 'Videos');
                        mc_path = fullfile(test_stims(t).folder,test_stims(t).name, 'MotionCorrection');
                        %stimevent_path = fullfile(stim_path, mouse_files{m}, Day, Cells(cc).name,[stimname,'.mat']);
                        stimevent_path = fullfile(stim_path, mouse_files{m}, Day, "test_vis",test_stims(t).name);
                        save_path = fullfile(output_path,mouse_files{m},Day,test_stims(t).name, "Soma_0");
                        file_results = dir(results_path);
                        file_results=file_results(contains({file_results.name}, '.txt'));
                        results_name = file_results(1).name;
                        if isdir(save_path)
                            organize_spine_stimevents_data_longstim_v2(results_name,stimname, mc_path, results_path,...
                            video_path,stimevent_path, save_path, stimfiles_names);
                       end
                    end
                end
            else
            stims = dir(fullfile(input_path,mouse_files{m},Day,Cells(cc).name,'Data'));
            stims=stims(~contains({stims.name}, '.'));
            stims=stims(~contains({stims.name}, 'X_'));
            for iii = 1:length(stims)
                stimname = stims(iii).name;
                disp(['Running segment: ',stimname])
                results_path = fullfile(input_path,mouse_files{m},Day,Cells(cc).name,'Data', stimname, 'Results');
                video_path = fullfile(input_path,mouse_files{m},Day,Cells(cc).name, 'Data',stimname, 'Videos');
                mc_path = fullfile(input_path,mouse_files{m},Day,Cells(cc).name, 'Data', stimname, 'MotionCorrection');
                stimevent_path = fullfile(stim_path, mouse_files{m}, Day, Cells(cc).name,[stimname,'.mat']);
                stimevent_path = fullfile(stim_path, mouse_files{m}, Day, Cells(cc).name,stimname);
                segments = dir(results_path);
                segments = segments(~[segments.isdir]);
                if ~isempty(segments)
                    for s = 1:length(segments)
                        results_name = segments(s).name;
                        ind = strfind(results_name,'Set');
                        if ~isempty(ind)
                            num = str2num(results_name(ind+3));
                            segmentname = [stimname,'_',num2str(num-1)];
                        else
                             segmentname = [stimname,'_',num2str(s-1)];
                        end
                        save_path = fullfile(output_path,mouse_files{m},Day,Cells(cc).name, segmentname);
                        if ~isdir(save_path)
                            organize_spine_stimevents_data_longstim_v2(results_name,stimname, mc_path, results_path,...
                                video_path,stimevent_path, save_path, stimfiles_names);
                        end
                    end
                end
            end
            
            end
        end
   end
 end


%% analyze processed data 
disp('Extract z scored data');
graph = 0;
for m = 1:length(mouse_files)
    disp(['Running mouse: ', mouse_files{m}])
   folders = dir(fullfile(output_path,mouse_files{m}));
   Days = {folders.name};
   Days=Days(~contains({folders.name}, '.'));
   for i = 1:length(Days)
        Day = Days{i};
        disp(['Running day: ',Day])
        Cells = dir(fullfile(output_path,mouse_files{m},Day));
        Cells=Cells(~contains({Cells.name}, '.'));
        for cc = 1:length(Cells)
            disp(['Running cell: ',Cells(cc).name])
            stims = dir(fullfile(output_path,mouse_files{m},Day,Cells(cc).name));
            stims=stims(~contains({stims.name}, '.'));
            for iii = 1:length(stims)
                disp(['Running segment: ',stims(iii).name])
                sessions = dir(fullfile(output_path,mouse_files{m},Day,Cells(cc).name, stims(iii).name));
                sessions=sessions(~contains({sessions.name}, '.'));
                for s = 1:length(sessions)
                    if ~exist(fullfile(sessions(s).folder, sessions(s).name, 'normalized_data_by_stim.mat'))
                    load(fullfile(sessions(s).folder, sessions(s).name, 'ROIdata.mat'))
                    ROIdata_scan = rois;
                     filename = sessions(s).name;
                    % Get Timestamps of frames
                    %open and import XML file
                    files = dir(fullfile(sessions(s).folder, sessions(s).name));
                    bruker_xml_file = files(contains({files.name}, '.xml') & ~contains({files.name}, 'VoltageRecording') & ~contains({files.name}, '._'));
                    file = fullfile(bruker_xml_file(1).folder, bruker_xml_file(1).name);
                    data = xml2struct_bruker(file);
        
                    %determine length
                    framelength = length(data.PVScan.Sequence.Frame);
                    %create vector of proper size for timestamps
                    timestamps=zeros(framelength,1);
                    %loop to extract timestamps, either relative or absolute
                    if (data.PVScan.Sequence.Frame{1, 2}.Attributes.relativeTime>0)
                        for i = 1:framelength
                            timestamps(i,1)=str2double(data.PVScan.Sequence.Frame{1, i}.Attributes.relativeTime);
                        end
                    else
                        for i = 1:framelength
                            timestamps(i,1)=str2double(data.PVScan.Sequence.Frame{1, i}.Attributes.absoluteTime);
                        end  
                    end
                    clear i

                    if length(unique(timestamps))>1
                    % Convert to deltaF/F, analyze by unique stim, and save output
                        if  (contains(sessions(s).name, 'Cell')...
                                || contains(sessions(s).name, 'Soma'))
                            
                             normalize_stimevents_soma_longstim(timestamps,ROIdata_scan,graph,files,framelength)
                         else
                             normalize_stimevents_spine_longstim(timestamps,ROIdata_scan,graph,files,framelength)
                        end
                        end
                    end
                end
            end
            
        end
   end
end
%%
Analyze OSI and DSI for each FOV stim 
disp('Analyze vector based tuning properties');
for m = 1:length(mouse_files)
    disp(['Running mouse: ', mouse_files{m}])
   folders = dir(fullfile(output_path,mouse_files{m}));
   Days = {folders.name};
   Days=Days(~contains({folders.name}, '.'));
   for i = 1:length(Days)
       Day = Days{i}; 
       disp(['Running day: ',Day])
        
        Cells = dir(fullfile(output_path,mouse_files{m},Day));
        Cells=Cells(~contains({Cells.name}, '.'));
        for cc = 1:length(Cells)
            disp(['Running cell: ',Cells(cc).name])
            stims = dir(fullfile(output_path,mouse_files{m},Day,Cells(cc).name));
            stims=stims(~contains({stims.name}, '.'));
            for iii = 1:length(stims)
                 disp(['Running segment: ',stims(iii).name])
                
                sessions = dir(fullfile(output_path,mouse_files{m},Day,Cells(cc).name, stims(iii).name));
                sessions=sessions(~contains({sessions.name}, '.'));
                for s = 1:length(sessions)

                    if ~exist(fullfile(sessions(s).folder, sessions(s).name, 'mean_amp_ori_analysis.mat'))
                    if exist(fullfile(sessions(s).folder, sessions(s).name, 'normalized_data_by_stim.mat'))
                    
             
                        load(fullfile(sessions(s).folder, sessions(s).name, 'normalized_data_by_stim.mat'),'median_amplitude', 'mean_amplitude', 'std_amplitude','unique_orientations', 't_test')
                        if(size(mean_amplitude,3) > 1)
                            ind = find_preferred_SF(mean_amplitude, t_test);
                            mean_amplitude = mean_amplitude(:,:,ind);
                            std_amplitude = std_amplitude(:,:,ind);
                            median_amplitude = median_amplitude(:,:,ind);
                           
                        end
                         orientation_tuning_vector_longstim('mean_amp_ori_analysis',mean_amplitude, std_amplitude, unique_orientations, fullfile(sessions(s).folder,sessions(s).name),Day);
                         orientation_tuning_vector_longstim('median_amp_ori_analysis',median_amplitude, std_amplitude, unique_orientations, fullfile(sessions(s).folder,sessions(s).name),Day);
                         orientation_tuning_von_mises('mean_amp',mean_amplitude, unique_orientations, fullfile(sessions(s).folder,sessions(s).name));
                         orientation_tuning_von_mises('median_amp',mean_amplitude, unique_orientations, fullfile(sessions(s).folder,sessions(s).name));

     
                      end
                end
                end
            end
        end
   end
end
