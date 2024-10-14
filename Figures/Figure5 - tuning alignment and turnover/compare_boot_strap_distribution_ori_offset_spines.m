function [tbl] = compare_boot_strap_distribution_ori_offset_spines(g,h,a,b,save_temp, days, numsubplots, n, session)
% perform boot strap of ori and osi offset by:
% 1) concatenating mean amp (8 by 10 trials) for d1 and d5
% 2) select 10 random trials for d1 and d5 per direction
% 3) find OSI and Ori pref for d1 and d5 and take offset
% 4) repeat N times per spine

unique_orientations = [0:45:315]; 
xtick_labels =[];
all_ori_difference = [];
all_dir_difference = [];
all_osi_difference = [];
all_dsi_difference = [];

all_days = [];
all_sessions = [];
count = 1;

for i = 1:length(save_temp)

        d1 = days{i};
        d5 = days{i+1};
        temp = save_temp{i};
        temp = splitvars(temp);
        
        
        d1_resp_inds = temp.([d1, '_resp'])==1;
        d5_resp_inds = temp.([d5, '_resp'])==1;
    
           
        % get visually responsive neurons on d1 and d5
        boot_temp = temp(d1_resp_inds&d5_resp_inds&strcmp(temp.session,session),:);
             
        subplot(2,numsubplots,n)
               
        all_ori_d1_boot =  boot_temp.( [d1, '_all_mean_amp_boot']);
        all_ori_d1_boot = cellfun(@(x) changed_to_0(x), all_ori_d1_boot, 'UniformOutput',false);
       
        all_ori_d5_boot = boot_temp.( [d5, '_all_mean_amp_boot']);
        all_ori_d5_boot = cellfun(@(x) changed_to_0(x), all_ori_d5_boot, 'UniformOutput',false);
        % permutation test to assess difference in orientation preference
        sig_inds_ori = [];
        ori_difference = [];

        sig_inds_tuning= [];
        tuning_corr = [];

        sig_inds_dir = [];
        dir_difference = [];

        high_sig_inds_osi = [];
        low_sig_inds_osi = [];
        osi_difference = [];

        high_sig_inds_dsi = [];
        low_sig_inds_dsi = [];
        dsi_difference = [];
        for ii = 1:length(all_ori_d1_boot)
            ori_d1_boot = all_ori_d1_boot{ii};
            ori_d5_boot = all_ori_d5_boot{ii};
            
            [Ori_pref_boot1,OSI_boot1,Dir_pref_boot1, DSI_boot1] = orientation_tuning_vector(ori_d1_boot', unique_orientations);
            [Ori_pref_boot5,OSI_boot5,Dir_pref_boot5, DSI_boot5] = orientation_tuning_vector(ori_d5_boot', unique_orientations);
%             best_shift_distribution = get_shifts_corr(ori_d1_boot,ori_d5_boot);
%             p = sum((best_shift_distribution-1)==0)/length(best_shift_distribution);
%             tuning_corr = [tuning_corr; corr(mean(ori_d1_boot,2), mean(ori_d5_boot,2))];
%             sig_inds_tuning = [sig_inds_tuning; p<0.05];
% 
            [true_diff,sig_inds]= permutation_test_ori(Ori_pref_boot1,Ori_pref_boot5);
            ori_difference = [ori_difference;true_diff];
            sig_inds_ori = [sig_inds_ori;sig_inds];

            [true_diff,sig_inds]= permutation_test_dir(Dir_pref_boot1,Dir_pref_boot5);
            dir_difference = [dir_difference;true_diff];
            sig_inds_dir = [sig_inds_dir;sig_inds];

            [true_diff,high_sig_inds,low_sig_inds]= permutation_test_selectivity(OSI_boot1,OSI_boot5);
            osi_difference = [osi_difference;true_diff];
            high_sig_inds_osi = [high_sig_inds_osi;high_sig_inds];
            low_sig_inds_osi = [low_sig_inds_osi;low_sig_inds];

            [true_diff,high_sig_inds,low_sig_inds]= permutation_test_selectivity(DSI_boot1,DSI_boot5);
            dsi_difference = [dsi_difference;true_diff];
            high_sig_inds_dsi = [high_sig_inds_dsi;high_sig_inds];
            low_sig_inds_dsi = [low_sig_inds_dsi;low_sig_inds];
            
        end
             
        figure(g)
        plot_angle_difference(numsubplots,n,ori_difference, sig_inds_ori,i, count,[0,90]), hold on;
        title([d1, ' to ', d5])
       
        figure(h)
        plot_angle_difference(numsubplots,n,dir_difference, sig_inds_dir,i, count,[0,180]), hold on;
        title([d1, ' to ', d5])
        
        

%         plot_angle_difference(numsubplots,n,dir_difference, sig_inds_dir,i,count,[0,180]),hold on;
%         title([d1, ' to ', d5])
              
        
        figure(a)
        plot_selectivity_difference(numsubplots,n,osi_difference, high_sig_inds_osi,low_sig_inds_osi,i, count), hold on;
        title([d1, ' to ', d5])

        figure(b)
        plot_selectivity_difference(numsubplots,n,dsi_difference, high_sig_inds_dsi,low_sig_inds_dsi,i, count), hold on;
        title([d1, ' to ', d5])

     
        all_ori_difference = [all_ori_difference;ori_difference];
        all_osi_difference = [all_osi_difference; osi_difference];
        all_dir_difference = [all_dir_difference; dir_difference];
        all_dsi_difference = [all_dsi_difference; dsi_difference];
        all_sessions = [all_sessions; repmat({session},size(ori_difference))];
        all_days = [all_days; repmat({[d1,d5]},size(ori_difference))];

        count = count + 3;
    end
subplot(2,numsubplots,n)
title(session)
xticks([1.5,4.5])
xticklabels(xtick_labels)
xlim([0,6])
%ylim([0,90])

if n == 1
    ylabel({'Tuning corr'})
end
tbl = table(all_dsi_difference, all_osi_difference, all_ori_difference,all_dir_difference,all_sessions,all_days);
end
%%
function x = changed_to_0(x)
  x(x<0) = 0;
end
%%
function mean_orientation_deg = calculate_mean_ori(Ori_pref_boot)
    orientations_rad = deg2rad(Ori_pref_boot);
    bin_edges = linspace(0,pi,10000);
    ori_hist = histcounts(orientations_rad,bin_edges);
    
    % Calculate the circular mean
    R_ori = sum(ori_hist.*exp(bin_edges(2:end)*2*1i));
    mean_orientation_deg = mod(atan2(imag(R_ori),real(R_ori)),2*pi)*(90/pi);
end
%%
function mean_orientation_deg = calculate_mean_dir(Dir_pref_boot)
    orientations_rad = deg2rad(Dir_pref_boot);
    bin_edges = linspace(0,2*pi,10000);
    ori_hist = histcounts(orientations_rad,bin_edges);
    
    % Calculate the circular mean
    R_ori = sum(ori_hist.*exp(bin_edges(2:end)*1i));
    mean_orientation_deg = mod(atan2(imag(R_ori),real(R_ori)),2*pi)*(180/pi);
end
%%
function offset = find_direction_offset(data1,data2)
    offset = data1-data2;
    if abs(offset)>180
        if offset<0
            offset = 360+offset;
        else
            offset = 360-offset;
        end
    end
end
%%
function offset = find_orientation_offset(data1,data2)
    offset = data1-data2;
    if abs(offset)>90
        if offset<0
            offset = 180+offset;
        else
            offset = 180-offset;
        end
    end
end

%%
function [true_diff,sig_inds]= permutation_test_ori(Ori_pref_boot1,Ori_pref_boot5)
        mean_deg1 = calculate_mean_ori(Ori_pref_boot1);
        mean_deg5 = calculate_mean_ori(Ori_pref_boot5);
        true_diff = abs(find_orientation_offset(mean_deg5,mean_deg1));
        diff1 = arrayfun(@(x) find_orientation_offset(mean_deg1,x), Ori_pref_boot5);
        diff5 = arrayfun(@(x) find_orientation_offset(mean_deg5,x), Ori_pref_boot1);

        p1 = min(sum(diff1>0),sum(diff1<0))/length(diff1);
        p2 = min(sum(diff5>0),sum(diff5<0))/length(diff5);

        if p1<0.05||p2<0.05
            sig_inds = 1;
        else
            sig_inds = 0;
        end
end

%%
function [true_diff,sig_inds]= permutation_test_dir(Dir_pref_boot1,Dir_pref_boot5)
        mean_deg1 = calculate_mean_dir(Dir_pref_boot1);
        mean_deg5 = calculate_mean_dir(Dir_pref_boot5);
        true_diff = abs(find_direction_offset(mean_deg5,mean_deg1));
        diff1 = arrayfun(@(x) find_direction_offset(mean_deg1,x), Dir_pref_boot5);
        diff5 = arrayfun(@(x) find_direction_offset(mean_deg5,x), Dir_pref_boot1);

        p1 = min(sum(diff1>0),sum(diff1<0))/length(diff1);
        p2 = min(sum(diff5>0),sum(diff5<0))/length(diff5);

        if p1<0.05||p2<0.05
            sig_inds = 1;
        else
            sig_inds = 0;
        end
end

%%
function [true_diff,high_sig_inds, low_sig_inds]= permutation_test_selectivity(boot1,boot5)
       
        high_sig_inds = 0;
        low_sig_inds = 0;
        mean_sel1 = nanmean(boot1);
        mean_sel5 = nanmean(boot5);
        true_diff = mean_sel5-mean_sel1;
        diff1 = mean_sel5-boot1;
        p_low = nansum(diff1>0)/sum(~isnan(diff1));
        p_high = nansum(diff1<0)/sum(~isnan(diff1));
        if p_low<0.05
            low_sig_inds = 1;
        elseif p_high<0.05
            high_sig_inds = 1;
        end
       
end

%% 
function plot_angle_difference(numsubplots,n,difference, sig_inds,i, count,ylims)
    colors = NaN(size(sig_inds));
    colors(sig_inds==1) = 1;
    colors(sig_inds==0) = 2;

    subplot(2,numsubplots,n)
    swarmchart(count*ones(size(difference)),  difference, [], colors,'filled','MarkerFaceAlpha',0.5,'SizeData',10); hold on;
    
    bar(count+1, mean(difference));
    errorbar(count+1, mean(difference), std(difference)/sqrt(length(difference)), '.k');
    ylim(ylims)
    subplot(2,2*numsubplots,(2*n-1)+(i-1)+ 2*numsubplots)
    pie([sum(sig_inds==1),sum(sig_inds==0)])

end
      

%%
function plot_selectivity_difference(numsubplots,n,difference, high_sig_inds,low_sig_inds,i, count)
    colors = NaN(size(high_sig_inds));    
    insig_inds = (high_sig_inds|low_sig_inds)==0;
    colors(insig_inds) = 1;
    colors(high_sig_inds==1) = 2;
    colors(low_sig_inds==1) = 3;
    subplot(2,numsubplots,n)
    yyaxis left
    swarmchart(count*ones(size(difference)), difference, [], colors,'filled','MarkerFaceAlpha',0.5,'SizeData',10); hold on;
    ylim([-0.6,0.6])

    yyaxis right
    bar(count+1, mean(difference));
    [p,a] = signrank(difference)
    plot_sig(count+1,0.2,p);
    errorbar(count+1, mean(difference), std(difference)/sqrt(length(difference)), '.k');
    ylim([-0.2,0.2])
    subplot(2,2*numsubplots,(2*n-1)+(i-1)+ 2*numsubplots)
    pie([sum(high_sig_inds),sum(low_sig_inds),sum(insig_inds)]);

end



