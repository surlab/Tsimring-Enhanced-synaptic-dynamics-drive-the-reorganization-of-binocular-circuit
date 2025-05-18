%% Run FOV alignment between each dendritic segment 
%close all 
clear all
thresh_um = 1;
run_gui = 1;

mouse_path = '/Volumes/GoogleDrive-108846495442099470486/My Drive/Sur Lab/Development project/Binocular_Matching/Spine_imaging/Chronic Imaging/FOV_alignment/';

mice = dir(mouse_path);
mice = mice(~ismember({mice.name},{'.','..', '.DS_Store', 'ReadMe.txt'}));
count = 0;
correct_label = [];
fraction_retained = [];
for m = 1:length(mice)
    is_blinded = 0;
    cells = dir(fullfile(mice(m).folder, mice(m).name));
    cells = cells(~ismember({cells.name},{'.','..', '.DS_Store', 'ReadMe.txt'}));
    T = table;
    T = table;
    blind_folder = fullfile(mice(m).folder, mice(m).name, 'key_to_day_folders.csv');
    if exist(blind_folder,"file")
        is_blinded = 1;
        tbl = readtable(blind_folder);
       
    end

    for c = 1:length(cells)
        %if contains(cells(c).name, 'Cell')
        disp(cells(c).name)
        days = dir(fullfile(cells(c).folder, cells(c).name));
        if ~isempty(find(contains({days.name}, {'Compare  '})))
            days = days(~ismember({days.name},{'.','..', '.DS_Store'}));
            days = days(~contains({days.name},{'Compare'}));
            inds = [];
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
            
            x_t = [];
            cell_segments = [];
           
            day_comp_d1_d5 = ['Compare  ', days(1).name, '-', days(length(days)).name];
            savepath1 = fullfile(cells(c).folder,cells(c).name, day_comp_d1_d5);
            dends = dir(savepath1);
            dends = dends(contains({dends.name},{'Dend'}));
            
            for i = 1:length(dends)
                dends_name = dends(i).name;
                d1_d5 = load(fullfile(savepath1, dends_name, ['alignment_based_on_fiducial.mat']), 'squeeze_mat', 'modify_mat', 'fid_file1', 'fid_file2', 'stats1_seg', 'stats2_seg');
                temp = d1_d5.modify_mat;
                num_retained_d1_d5_5d_interval = sum(arrayfun(@(x,y) x>0&y>0, temp(:,1), temp(:,2)));
                num_lost_d1_d5_5d_interval = sum(arrayfun(@(x,y) x>0&y==0, temp(:,1), temp(:,2)));
                num_added_d1_d5_5d_interval = sum(arrayfun(@(x,y) x==0&y>0, temp(:,1), temp(:,2)));
        
                d1_d5_5d_interval = temp(arrayfun(@(x,y) x>0&y>0, temp(:,1), temp(:,2)),:);
           
                %perform daily interval alignment
                day_comp = ['Compare: ', days(1).name, '-', days(2).name];
                num_lost = 0;
                num_added = 0;
                savepath2 = fullfile(cells(c).folder,cells(c).name, day_comp, dends_name);
                if isdir(savepath2)
                    d1_d2 = load(fullfile(savepath2, ['alignment_based_on_fiducial.mat']), 'squeeze_mat', 'modify_mat', 'fid_file1', 'fid_file2', 'stats1_seg', 'stats2_seg');
                    temp = d1_d2.modify_mat;
                    num_lost = sum(arrayfun(@(x,y) x>0&y==0, temp(:,1), temp(:,2)));
                    vars = ["D1", "D2"];
                    retained_spine_alignment = array2table(temp(arrayfun(@(x,y) x~=-1&y>0, temp(:,1), temp(:,2)),:), 'VariableNames',vars);
                    for d = 2:length(days)-1
                        day_comp = ['Compare: ', days(d).name, '-', days(d+1).name];
                        savepath3 = fullfile(cells(c).folder,cells(c).name, day_comp, dends_name);
                        if isdir(savepath3)
                            d3_d4 = load(fullfile(savepath3, ['alignment_based_on_fiducial.mat']), 'squeeze_mat', 'modify_mat', 'fid_file1', 'fid_file2', 'stats1_seg', 'stats2_seg');             
                            temp = d3_d4.modify_mat;
                            vars = [strcat("D",num2str(d)), strcat("D",num2str(d+1))];
                            feature = strcat("D",num2str(d));
                            retained_added_d3_d4 = array2table(temp(arrayfun(@(x,y) x~=-1&y>0, temp(:,1), temp(:,2)),:),'VariableNames',vars);
                            retained_spine_alignment = outerjoin(retained_spine_alignment, retained_added_d3_d4, "Keys",feature, "MergeKeys",true);
    
                             d1_d5_daily_interval = table2array(retained_spine_alignment(:,["D1", strcat("D",num2str(length(days)))]));
                            num_added_d1_d5_daily_interval =  sum(arrayfun(@(x,y) (x==0|isnan(x))&y>0, d1_d5_daily_interval(:,1), d1_d5_daily_interval(:,2)));
                            num_lost_d1_d5_daily_interval =  sum(arrayfun(@(x,y) x>0&(y==0|isnan(y)), d1_d5_daily_interval(:,1), d1_d5_daily_interval(:,2)));
                    
                            d1_d5_daily_interval= d1_d5_daily_interval(arrayfun(@(x,y) x>0&y>0,  d1_d5_daily_interval(:,1),    d1_d5_daily_interval(:,2)),:);
                            num_retained_d1_d5_daily_interval = size(d1_d5_daily_interval,1);
                            d1_d5_daily_interval = arrayfun(@(x,y) {strcat(num2str(x),'_',num2str(y))},  d1_d5_daily_interval(:,1),d1_d5_daily_interval(:,2));
                            d1_d5_5d_interval = arrayfun(@(x,y) {strcat(num2str(x),'_',num2str(y))},  d1_d5_5d_interval(:,1), d1_d5_5d_interval(:,2));
                            if xor(num_retained_d1_d5_daily_interval == 0, num_retained_d1_d5_5d_interval == 0) 
                                same = 0;
                                total = length(d1_d5_daily_interval) + length(d1_d5_5d_interval);
                            else
                           
                                same = length(intersect(d1_d5_daily_interval,d1_d5_5d_interval));
                                total = length(union(d1_d5_daily_interval,d1_d5_5d_interval));
                            end
                                                   %total = length(d1_d5_daily_interval);
                                                  
                            correct_label = [correct_label, same/total];
                            %figure(fig1) 
                            total_5d_interval = num_retained_d1_d5_5d_interval+num_lost_d1_d5_5d_interval;
                            total_daily_interval = num_retained_d1_d5_daily_interval+num_lost_d1_d5_daily_interval+num_lost;
                            fraction_retained = [fraction_retained; num_retained_d1_d5_daily_interval/total_daily_interval, num_retained_d1_d5_5d_interval/total_5d_interval ];
                            count = count + 1;
                                
                        end
                    end
                end
            end
        end
     end
end
%%
close all
figure
subplot(1,2,1)
bar(nanmean(correct_label)); hold on;
errorbar(nanmean(correct_label), nanstd(correct_label)/sqrt(length(correct_label)), 'k')
scatter(ones(size(correct_label)), correct_label, 10,'filled', 'k', 'jitter', 'on', 'jitterAmount', 0.05)
ylabel('Fraction')
title('Fraction of Correctly Labeled Spines (per dendritic segment)')

subplot(1,2,2)
 bar(nanmean(fraction_retained)); hold on;
 arrayfun(@(x,y) plot([1,2], [x,y],'-k', 'Marker', '.'), fraction_retained(:,1), fraction_retained(:,2))
xlabel('Dendrites')
 ylabel('Fraction of retained spines (D1 to D10)')
 xticklabels({'5-day', '10-day'})
 title('Fraction of Retained Spines (per dendritic segment)')