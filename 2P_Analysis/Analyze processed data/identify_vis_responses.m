%% identify visual responsiveness of spines
function identify_vis_responses(path,analyzefile, savefile, mean_amp_thresh, p_val_thresh)
savepath = fullfile(path, 'Analyzed Data');
load(fullfile(savepath,analyzefile));
savefile = fullfile(savepath,savefile);
all_vis_stims = [];
for i = 1:length(all_stims)
    t_test = all_stims(i).all_ttest;
    mean_amp = all_stims(i).all_mean_amp;
    large_pre_mean = all_stims(i).all_large_pre_mean;
    artifact_trial = all_stims(i).all_artifact_trial;
    gof = all_stims(i).all_GOF;
    ind = cellfun(@(x) find(x == max(x),1),mean_amp, 'UniformOutput',false);
    ind2 = cellfun(@(x) isempty(x),ind);
%     criteria = zeros(size(ind));
%     criteria(ind2==0) = cellfun(@(a,b,c,d,e)int8(a(e)>mean_amp_thresh) + int8(b(e)<p_val_thresh) +int8((c(e)+d(e)<=5)),...
%         mean_amp(ind2==0),t_test(ind2==0),large_pre_mean(ind2==0), artifact_trial(ind2==0),ind(ind2==0));
    
  
    criteria=cellfun(@(a,b,c,d) max(int8(a>mean_amp_thresh) + int8(b<p_val_thresh) +int8((c+d)<=5)), mean_amp,t_test,large_pre_mean, artifact_trial);
    shaft_corr = all_stims(i).all_shaft_corr;

    %stims.resp = criteria >=3 & shaft_corr<=0.5;
    stims.resp = criteria >=3;   
    criteria=cellfun(@(c,d) max(int8((c+d)<=5)),large_pre_mean, artifact_trial);
   
    stims.nonresp = criteria>=1;    

   non_resp = sum(stims.nonresp)
   resp_spines = sum(stims.resp)
   stims.session = all_stims(i).name;
   all_vis_stims = [all_vis_stims,stims];
end
save(savefile, 'all_vis_stims');
end
