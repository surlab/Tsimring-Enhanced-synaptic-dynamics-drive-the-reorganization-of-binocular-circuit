function analysis_spine_2p_data_function(input_path, mouse_files)
% Both Kyle and Katya worked on this script 

% program is designed to import 2P data. Select data, and then select
% timestamps. If anaylzing multiple recordings in a
% cocatenated stack, enter number of frames for each recording in the order
% in which they are in the stack as a vector. 

%scans: vector of frames in each scan in the order they are cocatenated



%% analyze processed data 
disp('Extract z scored data');
graph = 0;
for m = 1:length(mouse_files)
    disp(['Running mouse: ', mouse_files{m}])
   folders = dir(fullfile(input_path,mouse_files{m}));
   Days = {folders.name};
   Days=Days(~contains({folders.name}, '.'));
   for i = 1:length(Days)
        Day = Days{i};
        disp(['Running day: ',Day])
        Cells = dir(fullfile(input_path,mouse_files{m},Day));
        Cells=Cells(~contains({Cells.name}, '.'));
        for cc = 1:length(Cells)
            disp(['Running cell: ',Cells(cc).name])
            stims = dir(fullfile(input_path,mouse_files{m},Day,Cells(cc).name));
            stims=stims(~contains({stims.name}, '.'));
            for iii = 1:length(stims)
                disp(['Running segment: ',stims(iii).name])
                sessions = dir(fullfile(input_path,mouse_files{m},Day,Cells(cc).name, stims(iii).name));
                sessions=sessions(~contains({sessions.name}, '.'));
                for s = 1:length(sessions)
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

%%
%Analyze OSI and DSI for each FOV stim 
disp('Analyze vector based tuning properties');
for m = 1:length(mouse_files)
    disp(['Running mouse: ', mouse_files{m}])
   folders = dir(fullfile(input_path,mouse_files{m}));
   Days = {folders.name};
   Days=Days(~contains({folders.name}, '.'));
   for i = 1:length(Days)
       Day = Days{i}; 
       disp(['Running day: ',Day])
        
        Cells = dir(fullfile(input_path,mouse_files{m},Day));
        Cells=Cells(~contains({Cells.name}, '.'));
        for cc = 1:length(Cells)
            disp(['Running cell: ',Cells(cc).name])
            stims = dir(fullfile(input_path,mouse_files{m},Day,Cells(cc).name));
            stims=stims(~contains({stims.name}, '.'));
            for iii = 1:length(stims)
                 disp(['Running segment: ',stims(iii).name])
                
                sessions = dir(fullfile(input_path,mouse_files{m},Day,Cells(cc).name, stims(iii).name));
                sessions=sessions(~contains({sessions.name}, '.'));
                for s = 1:length(sessions)

                      if exist(fullfile(sessions(s).folder, sessions(s).name, 'normalized_data_by_stim.mat'))
                    
             
                        load(fullfile(sessions(s).folder, sessions(s).name, 'normalized_data_by_stim.mat'),'median_amplitude', 'mean_amplitude', 'std_amplitude','unique_orientations', 't_test')
                        if(size(mean_amplitude,3) > 1)
                            ind = find_preferred_SF(mean_amplitude, t_test);
                            mean_amplitude = mean_amplitude(:,:,ind);
                            std_amplitude = std_amplitude(:,:,ind);
                           
                        end
                         orientation_tuning_vector_longstim('mean_amp_ori_analysis',mean_amplitude, std_amplitude, unique_orientations, fullfile(sessions(s).folder,sessions(s).name));
     
                      end
                end
            end
        end
   end
end
