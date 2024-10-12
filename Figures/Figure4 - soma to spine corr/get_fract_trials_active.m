close all
unresp_cib = {'000'}; % unresponsive to all sessions
feature = 'all_active_spine_trials_smooth';
temp = all_stim_table(strcmp([all_stim_table.session{:}], 'binoc'),:);
temp.mouse_cell_day_session = strcat(temp.mouse_cell,temp.days, [temp.session{:}]');

unresp_soma = temp(temp.soma_resp==0,:);
resp_soma = temp(temp.soma_resp==1,:);
somas = {resp_soma};
figure("Position",[582,297,206,253])
all_soma_cas_active_trial = [];
all_soma_cas_sum_trial = [];

for s = 1:length(somas)
    soma = somas{s};
    soma_mean_resp = soma.soma_mean_resp;
    spine_trials = soma.(feature);
    inc_trials = soma.all_included_trial;
    peak_dir = cellfun(@(x) find(x==max(x)), soma_mean_resp, 'UniformOutput',false);

    sum_inc_trials = cellfun(@(x) sum(x,4), inc_trials, 'UniformOutput',false);
    spine_trials = cellfun(@(x,y) x&y, spine_trials,inc_trials, 'UniformOutput',false); %make sure active trial is an included one
    
    fract_active_trials_pref_dir = cellfun(@(x,y,z) sum(squeeze(x(1,y,:,:)))/z(y), spine_trials,peak_dir, sum_inc_trials);
    fract_active_trials = cellfun(@(x,y,z) sum(squeeze(x),2)./z', spine_trials,peak_dir, sum_inc_trials,'UniformOutput',false);

    active_spine_trials_pref_dir = cellfun(@(x,y) max(squeeze(x(1,y,:,:))), spine_trials,peak_dir);
    active_spine_trials = cellfun(@(x) max(squeeze(x),[],2), spine_trials, 'UniformOutput',false);

    sum_spine_trials_pref_dir = cellfun(@(x,y) sum(squeeze(x(1,y,:,:))), spine_trials,peak_dir);
    sum_spine_trials = cellfun(@(x) sum(squeeze(x),2), spine_trials, 'UniformOutput',false);

    % fraction of spines that were active for at least 1 trial 
    subplot(1,2,1)
    title('Fraction of spines were active (1 trial) to somas pref dir')
    G = findgroups(soma.mouse_cell_day_session);
    mean_ca_rate = splitapply(@mean, active_spine_trials_pref_dir,G);
    all_soma_cas_active_trial = [all_soma_cas_active_trial, {mean_ca_rate}];
    bar(s, mean(mean_ca_rate)); hold on;
    scatter(s*ones(size(mean_ca_rate)), mean_ca_rate,'k','jitter','on', 'jitterAmount', 0.05); 
    if s == length(somas)
       xticks([1:2]);
       xticklabels({'unresp soma', 'resp soma'})
       ylabel('Fraction')
      % ylim([0,1])
    end

    % fraction of trials that spines were active for      
    subplot(1,2,2)
    title('# of trials that spine was active to somas pref dir')
    [G,id] = findgroups(soma.mouse_cell_day_session);
    mean_ca_rate = splitapply(@mean, sum_spine_trials_pref_dir,G);
    all_soma_cas_sum_trial = [all_soma_cas_sum_trial, {mean_ca_rate}];
    bar(s, mean(mean_ca_rate)); hold on;
    scatter(s*ones(size(mean_ca_rate)), mean_ca_rate,'k','jitter','on', 'jitterAmount', 0.05); 
    if s == length(somas)
       xticks([1:2]);
       xticklabels({'unresp soma', 'resp soma'})
       ylabel('Count')
    end
    %ylim([0,5]);
    
    % plot soma's mean response against spines -> maybe try to make an
    sum_spine_trials_by_soma = splitapply(@mean,  cell2mat(sum_spine_trials')', G);
    soma_mean_resp_by_soma = splitapply(@mean,  cell2mat(soma_mean_resp), G);
    active_spine_trials_by_soma = splitapply(@mean, cell2mat(active_spine_trials')',G);
    fract_active_spine_trials_by_soma = splitapply(@mean, cell2mat(fract_active_trials')',G);

    %averaged tuning curve
    for i = 1:size(soma_mean_resp_by_soma,1)
        figure 
        
        subplot(2,1,1)
        yyaxis left
        plot(soma_mean_resp_by_soma(i,:)), hold on;
        ylabel('Soma mean amplitude')
        yyaxis right
        plot(fract_active_spine_trials_by_soma(i,:))
        ylabel('# of active trials/spine')
        title(id(i))

        
%         subplot(2,1,2)
%         yyaxis left
%         plot(soma_mean_resp_by_soma(i,:)), hold on;
%         ylabel('Soma mean amplitude')
%         yyaxis right
%         plot(active_spine_trials_by_soma(i,:))
%         ylabel('fraction of active spines')
    end
end
[p,a] = ranksum(all_soma_cas_sum_trial{1}, all_soma_cas_sum_trial{2}) %yes more spines active to soma's pref direction on responsive somas than unresponsive somas
[p,a] = ranksum(all_soma_cas_active_trial{1}, all_soma_cas_active_trial{2}) %yes spines active on more trials to soma's pref direction on responsive somas than unresponsive somas
