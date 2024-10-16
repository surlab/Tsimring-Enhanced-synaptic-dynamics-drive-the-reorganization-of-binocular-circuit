%% Compare ROIs between two time points
clear all
close all
path = 'C:/Users/Katya-PC/Dropbox (MIT)/Sur Lab/Projects/Development project/Binocular Matching/Spine imaging/';
savepath = fullfile(path, 'Chronic Imaging', 'FOV_alignment');

files = dir(fullfile(savepath));
files = files(~contains({files.name}, '.'));

for m = 1:length(files)
    cells = dir(fullfile(files(m).folder, files(m).name));
    cells = cells(~contains({cells.name},'.'));
    for c = 1:length(cells)
        days = dir(fullfile(cells(c).folder, cells(c).name));
        days = days(~contains({days.name},'.'));
        inds = randperm(length(days));
        
        for d = 1:length(inds)
            for dd = d+1:length(inds)
                day_comp = [ 'Compare ', days(inds(d)).name, '_', days(inds(dd)).name];
                fovs_1 = dir(fullfile(days(inds(d)).folder, days(inds(d)).name));
                fovs_1 = fovs_1(~contains({fovs_1.name},'.'));
          
                fovs_2 = dir(fullfile(days(inds(dd)).folder, days(inds(dd)).name));
                fovs_2 = fovs_2(~contains({fovs_2.name},'.'));

                save_fovs = intersect({fovs_1.name}, {fovs_2.name});
                for f = 1:length(save_fovs)

                    fov_1_path = fullfile(days(inds(d)).folder, days(inds(d)).name, save_fovs(f));
                    fov_2_path = fullfile(days(inds(dd)).folder, days(inds(dd)).name, save_fovs(f));

                    
                    rois1 = dir(fullfile(fov_1_path{1}, '*.zip'));
                    rois1 = rois1(~contains({rois1.name}, 'dend'));

                    rois2 = dir(fullfile(fov_2_path{1}, '*.zip'));
                    rois2 = rois2(~contains({rois2.name}, 'dend'));
                    
                    comp_rois = intersect({rois1.name}, {rois2.name});
                    savefile = fullfile(cells(c).folder, cells(c).name, day_comp);
                    if ~isfolder(savefile)
                        mkdir(savefile);
                    end
                    for r = 1:length(comp_rois)
                        
                        rois1_seg = ReadImageJROI(fullfile(fov_1_path{1}, comp_rois{r}));
                        rois1_seg = rois1_seg(1:3:end);
                        rois2_seg = ReadImageJROI(fullfile(fov_2_path{1}, comp_rois{r}));
                        rois2_seg = rois2_seg(1:3:end);
                    
                        [f1,map] = imread(fullfile(fov_1_path{1}, [save_fovs{f}, '.tif']));
                        figure("Name",[days(inds(d)).name, ': ', save_fovs{f}]) 
                        fig1 = imshow(f1,map); hold on;
                        imcontrast
                        axis('image', 'on'); 
                        for rr = 1:length(rois1_seg)
                            coords = rois1_seg{rr}.vfEllipsePoints;
                            aspect_ratio = rois1_seg{rr}.fAspectRatio;
                             xc = (coords(1)+coords(3))/2;
                             yc = (coords(2)+coords(4))/2;
                            draw_ellipse(xc,yc,aspect_ratio);
                            
                            txt = ['\leftarrow ' num2str(rr)];
                            text(xc,yc,txt, 'Color','y');
                        end
                        [f2,map] =  imread(fullfile(fov_2_path{1}, [save_fovs{f}, '.tif']));
                        figure("Name",[days(inds(dd)).name, ': ', save_fovs{f}])
                        fig2 = imshow(f2,map); hold on;
                        imcontrast;
                        axis('image', 'on'); 
                         for rr = 1:length(rois2_seg)
                            coords = rois2_seg{rr}.vfEllipsePoints;
                            aspect_ratio = rois2_seg{rr}.fAspectRatio;
                             xc = (coords(1)+coords(3))/2;
                             yc = (coords(2)+coords(4))/2;
                            draw_ellipse(xc,yc,aspect_ratio);
                            
                            txt = ['\leftarrow ' num2str(rr)];
                            text(xc,yc,txt, 'Color','y');
                         end
                         max_length = length(rois1_seg)+ length(rois2_seg);
                        data = zeros(max_length, 2);
                        data( [1:length(rois1_seg)],1) = [1:length(rois1_seg)];
                        data([1:length(rois2_seg)],2) = [1:length(rois2_seg)];
                    
                        fig3 = figure;
                        uit = uitable(fig3,  'Data' , data, 'ColumnEditable' , [true, true],...
                            'ColumnName' , {days(inds(d)).name  ; days(inds(dd)).name  });
                        set(uit, 'CellEditCallback', 'data_roi_table = get(uit,''Data'')')
                        uiwait(fig3)

                        disp('Program execution has resumed');

                        
                        
                        data_roi=table(data_roi_table);
                        if ~exist(fullfile(savefile, [save_fovs{f}, '_',comp_rois{r}, '.csv']), 'file')

                            writetable(data_roi, fullfile(savefile, [save_fovs{f}, '_',comp_rois{r}, '.csv']));
                        end
                        close all
                        clear data data_roi
                        
                    end
                end

            end
        end

    end
end
%%
function draw_ellipse(xc,yc,aspect_ratio)
    % compute points corresponding to axis-oriented ellipse
   
    r1 = aspect_ratio;
    r2 = 1/r1;
    t = linspace(0, 2*pi, 200);
    xt = r1 * cos(t) + xc;
    yt = r2 * sin(t) + yc;
    % aply rotation by angle theta
%     cot = cos(theta); sit = sin(theta);
%     x = xt * cot - yt * sit;
%     y = xt * sit - yt * cot;
    % draw the curbe
    plot(xt, yt, '-');
end
%%
function waitForReturnKey
    currkey=0;
    % do not move on until enter key is pressed
    while currkey~=1
        pause; % wait for a keypress
        currkey=get(gcf,'CurrentKey'); 
        if strcmp(currkey,'return');
            currkey==1;
        else
            currkey==0;
        end
    end
end
            