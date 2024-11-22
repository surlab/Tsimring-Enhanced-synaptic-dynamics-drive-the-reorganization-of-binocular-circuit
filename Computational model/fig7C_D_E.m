clear all
ntrials = 20;
n_branch = 2; 
N = 64;
all_survtime = []; 
all_activity = [];
activity_lost =[]; 
activity_retained =[]; 
spinesomacorr = [];
lost_pertrial = cell(1,ntrials);
retained_pertrial = cell(1,ntrials);

for trial = 1:ntrials
    load(['main' num2str(trial) '.mat'])    
    lost_pertrial{trial} = [];
    retained_pertrial{trial} = [];
    for branch = 1:n_branch
        for spine_idx = 1:N
            % Get survival times and times when the spines that were present at t=0 were lost
            lost_time = lost_times{branch}{spine_idx}(1); % index 1 corresponds to spines present at t=0
            surv_time = surv_times{branch}{spine_idx}(1);
            activity_pre = stored_PRE{branch}{spine_idx}(1);
            all_survtime = [all_survtime surv_time];
            all_activity = [all_activity activity_pre];
            
            % Classify as lost or retained
            if lost_time(1)<T
                lost_pertrial{trial} = [lost_pertrial{trial} activity_pre];
            else
                retained_pertrial{trial} = [retained_pertrial{trial} activity_pre];
            end
            
            % Calculate spine-soma correlation
            % get normalized activity for the time spine was present
            surv_activity = saving_PRE{branch}(1:round(lost_time(1)/(T/1000)),spine_idx);
            norm_activity = surv_activity./max(surv_activity); 
            % normalize activity for soma
            surv_soma = saving_Ath(1:round(lost_time(1)/(T/1000)));
            norm_soma = surv_soma./max(surv_soma);
            spinesomacorr = [spinesomacorr corr(norm_activity,norm_soma)];
       end
    end
    % normalize and store activities
    all_spines = [lost_pertrial{trial} retained_pertrial{trial}];
    lost_pertrial{trial} = lost_pertrial{trial}./max(max(all_spines));
    retained_pertrial{trial} = retained_pertrial{trial}./max(max(all_spines));
    activity_lost = [activity_lost lost_pertrial{trial}];
    activity_retained = [activity_retained retained_pertrial{trial}];
end
% Calculate fraction of survival time
fraction_survtime = all_survtime./T;
survtime_lost = fraction_survtime(all_survtime<T-1);
survtime_retained = fraction_survtime(all_survtime>=T-1);

corr_lost = spinesomacorr(all_survtime<T-1);
corr_retained = spinesomacorr(all_survtime>=T-1);

figure()
scatter(activity_retained,survtime_retained,100,[129 129 129]./255,'o')
hold on
scatter(activity_lost,survtime_lost,100,[245 126 32]./255,'o')
set(gca,'color','none')
box off
axis square
activity_all = [activity_lost activity_retained];
survtime_all = [survtime_lost survtime_retained];
fitlm(activity_all,survtime_all)

figure
hold on
errorbar(0,[mean(corr_lost)],[],[std(corr_lost)],'k','LineWidth',.75,'CapSize',8)
errorbar(1,[mean(corr_retained)],[], [std(corr_retained)],'k','LineWidth',.75,'CapSize',8)
b1 = bar([0 1],[mean(corr_lost) mean(corr_retained)],.5);
[p_corr,h_corr,stats_corr] = ranksum(corr_lost,corr_retained);

figure
hold on
errorbar(0,[mean(activity_lost)],[],[std(activity_lost)],'k','LineWidth',.75,'CapSize',8)
errorbar(1,[mean(activity_retained)],[], [std(activity_retained)],'k','LineWidth',.75,'CapSize',8)
b1 = bar([0 1],[mean(activity_lost) mean(activity_retained)],.5);
[p_activity,h_activity,stats_activity] = ranksum(activity_lost,activity_retained);
