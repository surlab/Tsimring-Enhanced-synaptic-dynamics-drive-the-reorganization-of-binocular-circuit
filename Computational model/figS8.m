%%%%%%%
%
% Code for the paper 
% "Enhanced synaptic dynamics drive the reorganization 
% of binocular circuits in mouse visual cortex"
% for Fig S8
% - need to load results from main.m (adapted with different values)
%%%%%%%
clear all
%% Parameters
nrows = 3;
ncols = 3;
ntrials = 20;
N = 64; % Number of spines on each branch
n_ori = 8; % Number of preferred orientations
lin_angles = 0:360/n_ori:315; % linear angles in degrees
lin_angles_half = lin_angles(1:4);
rad_angles = deg2rad(lin_angles); % convert angles to radians
wmin = 0.05; wthresh = 1.5*wmin;
checkpoints = [1 1000]; % time checkpoints
timepoints = length(checkpoints);

%% Initialization
deltas = zeros(ntrials, 3, 2);  % Store deltas before and after plasticity (2nd dimension: 1 = ipsi vs bino view, 2 = contra vs bino, 3 = contra vs ipsi, 3rd dimension: 1 = before, 2 = after)
dists = cell(ntrials, 2);       % Store dists for before and after
corrs = cell(ntrials, 2);       % Store corrs for before and after
d_min = zeros(ntrials, 2);    % Store min distance for before and after
d_max = zeros(ntrials, 2);    % Store max distance for before and after
%% Analysis for mismatch and correlation
for row = 1:nrows
    for col = 1:ncols
        for trial= 1:ntrials
            % load example data here (for various sigma_het, sigma_heb)
            load(['data_' num2str(row) num2str(col) '_' num2str(trial) '.mat'])    
            % Loop over the two timepoints (before and after plasticity)
            for indtime = 1:timepoints
                timepoint = checkpoints(indtime);
                % Extract trial data for the current timepoint
                test_angles = [saving_angles{1}(timepoint, :); saving_angles{2}(timepoint, :)]; 
                test_W = [saving_W{1}(timepoint, :); saving_W{2}(timepoint, :)];
                test_eye = [saving_eye{1}(timepoint, :); saving_eye{2}(timepoint, :)];
                test_ampVM = [saving_ampVM{1}(timepoint, :); saving_ampVM{2}(timepoint, :)];
                test_zVM = [saving_zVM{1}(timepoint, :); saving_zVM{2}(timepoint, :)];
        
                 % test viewing
                [ipsi_left,ipsi_right,soma_ipsi] = test(1,test_angles,test_W,test_eye,test_ampVM,test_zVM,dMat,dSoma,1,randVM); % Ipsi viewing
                [contra_left,contra_right,soma_contra] = test(2,test_angles,test_W,test_eye,test_ampVM,test_zVM,dMat,dSoma,1,randVM); % Contra viewing
                [bino_left,bino_right,soma_bino] = test(3,test_angles,test_W,test_eye,test_ampVM,test_zVM,dMat,dSoma,1,randVM); % Binocular viewing
        
                % Normalize the firing rates
                fr_abs = [soma_ipsi; soma_contra; soma_bino];
                fr = fr_abs ./ max(fr_abs, [], 2); 
        
                % Combine results
                xtest_completeI = [ipsi_left; ipsi_right];
                xtest_completeC = [contra_left; contra_right];
                xtest_completeall = [bino_left, bino_right];
        
                % Calculate orientation preferences and deltas
                oripref = preferred_orientation(lin_angles, fr);
                delta = abs([oripref(1)-oripref(3), oripref(2)-oripref(3), oripref(1)-oripref(2)]);
                delta(delta > 90) = 180 - delta(delta > 90);
                deltas(trial, :, indtime) = delta / 90;
        
                cc_L = corr(bino_left); % Left branch correlations
                cc_R = corr(bino_right); % Right branch correlations
        
                dd_L = zeros(N);
                for i=1:N
                    [dd_L(i,:),sortIdx] = sort(abs(pos{1}(i)-pos{1}),'ascend');
                    cc_L(i,:) = cc_L(i,sortIdx);
                end
        
                dd_R = zeros(N);
                pair_dists = zeros(1,N);
                for i=1:N
                    [dd_R(i,:),sortIdx] = sort(abs(pos{2}(i)-pos{2}),'ascend');
                    cc_R(i,:) = cc_R(i,sortIdx);
                end
        
                dist_L = mean(dd_L,1);         corr_L = mean(cc_L,1);   % left branch
                dist_R = mean(dd_R,1);         corr_R = mean(cc_R,1);   % right branch
        
                [dd_tot_sort,sort_tot] = unique([dist_L, dist_R]);
                cc_tot = [corr_L corr_R];
                cc_tot_sort = cc_tot(sort_tot);
        
                dists{trial, indtime} = dd_tot_sort;
                corrs{trial, indtime} = cc_tot_sort;
                d_min(trial, indtime) = dd_tot_sort(1);
                d_max(trial, indtime) = dd_tot_sort(end);        
            end
        end
        %% Plotting
        % Mismatch 
        deltas_bef = squeeze(deltas(:, :, 1));  % Data before plasticity
        deltas_aft = squeeze(deltas(:, :, 2));  % Data after plasticity
        fig1 = figure;
        hold on
        errorbar(0,[mean(deltas_bef(:,1))],[],[std(deltas_bef(:,1))],'k','LineWidth',.75,'CapSize',8)
        errorbar(1.5,[mean(deltas_bef(:,2))],[], [std(deltas_bef(:,2))],'k','LineWidth',.75,'CapSize',8)
        errorbar(3,[mean(deltas_bef(:,3))],[],[std(deltas_bef(:,3))],'k','LineWidth',.75,'CapSize',8)
        errorbar(.6,[mean(deltas_aft(:,1))],[],[std(deltas_aft(:,1))],'k','LineWidth',.75,'CapSize',8)
        errorbar(2.1,[mean(deltas_aft(:,2))],[],[std(deltas_aft(:,2))],'k','LineWidth',.75,'CapSize',8)
        errorbar(3.6,[mean(deltas_aft(:,3))],[],[std(deltas_aft(:,3))],'k','LineWidth',.75,'CapSize',8)
        bar(0,mean(deltas_bef(:,1)),.5,'FaceColor',[51	112	181]/256);
        bar(1.5,mean(deltas_bef(:,2)),.5,'FaceColor', [199	92	44]/256);
        bar(3,mean(deltas_bef(:,3)),.5,'FaceColor',	[147, 54, 144]/256);
        b1 = bar(.6,mean(deltas_aft(:,1)),.5,'FaceColor',[15,114,186]/256); b1.FaceAlpha = 0.3;
        b2 = bar(2.1,mean(deltas_aft(:,2)),.5,'FaceColor',[215 85 39]/256); b2.FaceAlpha = 0.3;
        b3 = bar(3.6,mean(deltas_aft(:,3)),.5,'FaceColor',[147, 54, 144]/256); b3.FaceAlpha = 0.3;
        set(gca,'FontSize',22)
        set(gca,'color','none')
        axis square
        box off
        xticks([])
        % stat test
        [P1,H1,STATS1] = ranksum(deltas_bef(:,1),deltas_aft(:,1));
        [P2,H2,STATS2] = ranksum(deltas_bef(:,2),deltas_aft(:,2));
        [P3,H3,STATS3] = ranksum(deltas_bef(:,3),deltas_aft(:,3));
        % savefig(fig1,['mismatch' num2str(row) num2str(col) '.fig'])
        
        % Correlations vs spine distances
        d_min_avg = min(d_min(:, 1)); % if spine location changes through plasticity, adapt here
        d_max_avg = max(d_max(:, 1)); % if spine location changes through plasticity, adapt here
        dists_avg = d_min_avg:d_max_avg; 
        corrs_avg = zeros(ntrials,length(dists_avg),timepoints);
        for indtime = 1:timepoints
            for trial = 1:ntrials
                corrs_avg(trial,:,indtime) = interp1(dists{trial,indtime},corrs{trial,indtime},dists_avg);
                nan_vals = isnan(corrs_avg(trial,:,indtime));
                corrs_avg(trial,nan_vals,indtime) = mean(corrs_avg(trial, ~nan_vals, indtime));
            end
        end
        
        corr_mean_before = mean(corrs_avg(:,:,1),1);
        corr_std_before = std(corrs_avg(:,:,1), 0, 1);
        corr_mean_after = mean(corrs_avg(:,:,2),1);
        corr_std_after = std(corrs_avg(:,:,2), 0, 1);
        
        fig2 = figure; 
        hold on
        fill([dists_avg fliplr(dists_avg)], [corr_mean_before + corr_std_before fliplr(corr_mean_before - corr_std_before)], [0.5 0.8 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        fill([dists_avg fliplr(dists_avg)], [corr_mean_after + corr_std_after fliplr(corr_mean_after - corr_std_after)], [0. 0.3 0.], 'FaceAlpha', 0.2, 'EdgeColor', 'none')
        plot(dists_avg, corr_mean_before, '-o', 'color',[0.5 0.8 0.5], 'LineWidth', 2)
        plot(dists_avg, corr_mean_after, '-o', 'color', [0. 0.3 .0], 'LineWidth', 2)
        xlim([0 30])
        set(gca, 'FontSize', 20)
        set(gca,'color','none')
        axis square
        box off
        legend('Before plasticity', 'After plasticity')
        % savefig(fig2,['corrs_vs_dists' num2str(row) num2str(col) '.fig'])
    end
end
