%% Run FOV alignment between each dendritic segment 
close all 
clear all
thresh_um = 1;
run_gui = 1;


mouse_path = '/Users/ktsimring/Dropbox (MIT)/Sur Lab/Projects/Development project/Binocular Matching/Spine imaging/Chronic Imaging/FOV_alignment';
mouse_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/FOV_alignment/to_do/';
mouse_path = 'G:/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/FOV_alignment/to_do/'
mice = dir(mouse_path);
mice = mice(~ismember({mice.name},{'.','..', '.DS_Store', 'ReadMe.txt'}));

for m = 1:length(mice)
    is_blinded = 0;
    cells = dir(fullfile(mice(m).folder, mice(m).name));
    cells = cells(~ismember({cells.name},{'.','..', '.DS_Store', 'ReadMe.txt'}));
    T = table;
    blind_folder = fullfile(mice(m).folder, mice(m).name, 'key_to_day_folders.csv');
    if exist(blind_folder,"file")
        is_blinded = 1;
        tbl = readtable(blind_folder);
       
    end

    for c = 3:length(cells)
        if contains(cells(c).name, 'Cell')
        disp(cells(c).name)
        days = dir(fullfile(cells(c).folder, cells(c).name));
        days = days(~ismember({days.name},{'.','..', '.DS_Store'}));
        days = days(~contains({days.name},{'Compare'}));
        x_t = [];
        cell_segments = [];
        if is_blinded
             tbl_filt = tbl(ismember(tbl.Var1_3, cells(c).name),{'Var1_5';'Var1_4'});
             real_date = tbl_filt.Var1_4;
             blind_date = tbl_filt.Var1_5;
             [~,inds_real] = sort(real_date);
             blind_date_sort = blind_date(inds_real);
             inds = arrayfun( @(x)( find(ismember( {days.name} ,x)) ), blind_date_sort, 'UniformOutput',false);
        else
            [~,inds] = sort({days.name}); 
        end
        if iscell(inds)
            inds = [inds{:}];
        end

         days = days(inds);
        for d = 1:length(days)-1
            day_comp = ['Compare  ', days(d).name, '-', days(d+1).name];
            segments1 = dir(fullfile(days(d).folder, days(d).name));
            segments1 = segments1(~ismember({segments1.name},{'.','..', '.DS_Store'}));
    
            segments2 = dir(fullfile(days(d+1).folder, days(d+1).name));
            segments2 = segments2(~ismember({segments2.name},{'.','..', '.DS_Store'}));

            unique_FOV = intersect({segments1.name}, {segments2.name});
            for s = 1:length(unique_FOV)
               if isfolder(fullfile(days(d).folder, days(d).name,...
                   unique_FOV{s}, 'dendritic_distance'))&...
                   isfolder(fullfile(days(d+1).folder, days(d+1).name,...
                   unique_FOV{s}, 'dendritic_distance'))
                   fid_file1 = readtable(fullfile(days(d).folder, days(d).name,...
                       unique_FOV{s}, 'dendritic_distance', 'fiducial_distances.csv'));
                   fid_file2 = readtable(fullfile(days(d+1).folder, days(d+1).name,...
                       unique_FOV{s}, 'dendritic_distance', 'fiducial_distances.csv'));
                   stats1 = readtable(fullfile(days(d).folder, days(d).name,...
                       unique_FOV{s}, 'dendritic_distance', 'spine_stats.csv'));
                   stats2 = readtable(fullfile(days(d+1).folder, days(d+1).name,...
                       unique_FOV{s}, 'dendritic_distance', 'spine_stats.csv'));
                   
                   stats1_seg = stats1.source_file;
                   stats2_seg = stats2.source_file;
                   unique_seg = intersect(stats1_seg, stats2_seg);
                   for ss = 1:length(unique_seg)
                       in = input(['Do you want to align ' day_comp ' ' unique_FOV{s}...
                           unique_seg{ss} ' Y/N [Y]: '],'s')
                       if in =='N'
                            continue;
                       end
                       seg_name = erase(unique_seg{ss}, 'RoiSet');
                       seg_name = erase(seg_name, '.zip');
                       seg_name = str2num(seg_name)-1;
                       sub_seg1_inds = find(ismember(stats1_seg, unique_seg{ss}));
                       sub_seg2_inds = find(ismember(stats2_seg, unique_seg{ss}));
                       fid1 = fid_file1(sub_seg1_inds,:);
                       fid2 = fid_file2(sub_seg2_inds,:);
                       fid1 = table2array(fid1);
                       fid2 = table2array(fid2);
%                        if size(fid1,2)>1
%                            fid1 = nanmean(fid1,2);
%                            fid2 = nanmean(fid2,2);
%                        end
                       nec_ang1 = table2array(stats1(sub_seg1_inds, 'relativeNeckAngle'));
                       nec_ang2 = table2array(stats2(sub_seg2_inds, 'relativeNeckAngle'));
                       [align_mat, squeeze_mat] = FOV_alignment(fid1, fid2, nec_ang1, nec_ang2, thresh_um);
                      

                       savepath = fullfile(cells(c).folder,cells(c).name, day_comp,[unique_FOV{s},'_', num2str(seg_name)]);
                        if ~isdir(savepath)
                            mkdir(savepath);
                        end
                        if run_gui
                           fov1_path = fullfile(days(d).folder, days(d).name,...
                       unique_FOV{s});
                           fov2_path = fullfile(days(d+1).folder, days(d+1).name,...
                       unique_FOV{s});
                           roi_file = unique_seg{ss};
                           ref_file = [unique_FOV{s}, '.tif'];
                           modify_mat = runGUI(savepath, fov1_path,fov2_path, roi_file, ref_file, fid1, fid2, squeeze_mat);
                       else
                           modify_mat = [];
                       end
                        save(fullfile(savepath, ['alignment_based_on_fiducial.mat']), 'squeeze_mat', 'modify_mat', 'fid_file1', 'fid_file2', 'stats1_seg', 'stats2_seg');
                   end
               end
               end
            end
        end
        
    end
end

%% Compare distances b/w 2 timepoints
function [alignment_mat, squeeze_mat] = FOV_alignment(fid1, fid2, nec_ang1, nec_ang2, thresh_um)
% alignment matrix (n by m matrix)
% n = spines in timepoint1
% m = spines in timepoint2
% if spine n = spine m then alignment matrix(n,m) = 1
% 
% squeeze matrix (2 column matrix), 
% 1st column is spine indices in timepoint1
% 2nd column is spine indices in timpoint2 
% if squeeze matrix(i, 1) = 0, spine not present in timepoint1
% if squeeze matrix(i, 2) = 0, spine not present in timepoint 2
% if squeeze matrix(i, 1) = -1, spine not present in timepoint1, but its
% imaged outside of timepoint2's FOV 

num_spines1 = height(fid1);
num_spines2 = height(fid2);
alignment_mat = zeros(num_spines1,num_spines2);
squeeze_mat = zeros(num_spines1+num_spines2,2);
squeeze_mat(1:num_spines1,1) = 1:num_spines1;
squeeze_mat(num_spines1+1:num_spines2+num_spines1,2) = 1:num_spines2;
fid1 =  fid1(:,any(~isnan(fid1))); 
fid2 =  fid2(:,any(~isnan(fid2))); 
for f1 = 1:height(fid1)
    dist = nanmean(abs(fid2 - fid1(f1,:)),2);
    ind = find(dist <= thresh_um);
    f2 = find_spine(nec_ang2, nec_ang1, dist(ind),f1, ind);
    if f2>0
        alignment_mat(f1,f2) = 1;
        squeeze_mat(f1,2) = f2;
        squeeze_mat(num_spines1+f2,2) = 0;
    end
end

for i = 1:size(squeeze_mat,1)
    ind1 = squeeze_mat(i,1);
    ind2 = squeeze_mat(i,2);
    if xor(ind1>0,ind2>0)
        if ind1 == 0

            if max(fid2(i-num_spines1,:))>max(fid1, [],"all")
                squeeze_mat(i,1) = -1;
            end
        else
            if max(fid1(i,:))>max(fid2, [],"all")
                squeeze_mat(i,2) = -1;
            end
        end
    end

end
end


%% Recursive function to find spine in timepoint2 that is within distance and has same neck angle
% Finds minimum distance within threshold --> if neck angle matches then
% spits out spine2 index, else goes to the next minimum distance
function f2 = find_spine(nec_ang2, nec_ang1, dist,f1, ind)
    if ~isempty(dist)
        [min_dist, min_ind] = min(dist);
        f2 = ind(min_ind);
        if sign(nec_ang2(f2)) == sign(nec_ang1(f1))
            return;
        else
            temp_dist = dist(dist~=min_dist);
            temp_ind = ind(dist~=min_dist);
            f2 = find_spine(nec_ang2, nec_ang1,temp_dist,f1,temp_ind);
        end
    else
        f2 = 0;
        return;
    end
end

%%
function modified_mat = runGUI(savepath, fov_1_path,fov_2_path, rois_file, ref_file, fid1, fid2, squeeze_matrix)
 close all
 rois1_seg = ReadImageJROI(fullfile(fov_1_path,rois_file));
        rois1_seg = rois1_seg(1:3:end);
        rois2_seg = ReadImageJROI(fullfile(fov_2_path, rois_file));
        rois2_seg = rois2_seg(1:3:end);
    
        
        fig1 = figure("Name",[ref_file], "Position", [37,293, 500, 1000]); 
        
        subplot(2,1,1)
        
        [f1,map] = imread(fullfile(fov_1_path, ref_file));
        %[f1, map] = imresize(f1,2);
        imshow(f1,map); hold on;
        title('T1');
        imcontrast
        axis('image', 'on'); 
        
        subplot(2,1,2);
        [f2,map] =  imread(fullfile(fov_2_path, ref_file));
        %[f2, map] = imresize(f2,2);
        imshow(f2,map); hold on;
        title('T2');
        imcontrast
        axis('image', 'on'); 
        fig2 = uifigure;
        modified_mat = uitable(fig2,  'Data' , squeeze_matrix, 'ColumnEditable' , [true, true],...
            'ColumnName' , {"T1"  ; "T2"  }, 'CellEditCallback', @(src,event) get(src, 'Data'));
        fig1;
        fid1 =  fid1(:,any(~isnan(fid1))); 
        fid2 =  fid2(:,any(~isnan(fid2))); 
        for rr = 1:size(squeeze_matrix,1)
            disp(['Number:',num2str(rr), 'out of ', num2str(length(squeeze_matrix))]);
            ind1 = squeeze_matrix(rr,1);
            ind2 = squeeze_matrix(rr,2);
            if ind1 > 0 && ind2 > 0
               coords = rois1_seg{ind1}.vfEllipsePoints;
               aspect_ratio = rois1_seg{ind1}.fAspectRatio;
               xc = (coords(1)+coords(3))/2;
               yc = (coords(2)+coords(4))/2;
               subplot(2,1,1);
               txt = ['\leftarrow ', num2str(ind1), ', dist: ', num2str(round(fid1(ind1),1))];
               text(xc,yc,txt, 'Color','c');

               coords = rois2_seg{ind2}.vfEllipsePoints;
               aspect_ratio = rois2_seg{ind2}.fAspectRatio;
               xc = (coords(1)+coords(3))/2;
               yc = (coords(2)+coords(4))/2;
               subplot(2,1,2);
               txt = ['\leftarrow ', num2str(ind2),', dist: ', num2str(round(fid2(ind2),1))];
               text(xc,yc,txt, 'Color','c');
            elseif ind2 > 0 
               coords = rois2_seg{ind2}.vfEllipsePoints;
               aspect_ratio = rois2_seg{ind2}.fAspectRatio;
               xc = (coords(1)+coords(3))/2;
               yc = (coords(2)+coords(4))/2;
               subplot(2,1,2);
               txt = ['\leftarrow ', num2str(ind2), ', dist: ', num2str(round(fid2(ind2),1))];
               text(xc,yc,txt, 'Color','g');
            elseif ind1 > 0 
               coords = rois1_seg{ind1}.vfEllipsePoints;
               aspect_ratio = rois1_seg{ind1}.fAspectRatio;
               xc = (coords(1)+coords(3))/2;
               yc = (coords(2)+coords(4))/2;
               subplot(2,1,1);
               txt = ['\leftarrow ', num2str(ind1), ', dist: ', num2str(round(fid1(ind1),1))];
               text(xc,yc,txt, 'Color','r');
            end
             pause;
        
        
        
        end
        
        modified_mat = modified_mat.Data;
        %uiwait(fig2)
        savefig(fullfile(savepath, 'Comparison_fig.fig'));
         disp('Program execution has resumed');


end 

             
