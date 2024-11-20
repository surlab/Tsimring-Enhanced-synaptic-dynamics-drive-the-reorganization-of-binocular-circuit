%% identify visual responsiveness of spines
function identify_vis_responses(savepath,analyzefile, savefile, mean_amp_thresh, p_val_thresh)
load(fullfile(savepath,analyzefile));
savefile = fullfile(savepath,savefile);
all_vis_stims = [];
for i = 1:length(all_stims)
    t_test = all_stims(i).all_ttest;
    mean_amp = all_stims(i).all_mean_amp;
    large_pre_mean = all_stims(i).all_large_pre_mean;
    artifact_trial = all_stims(i).all_artifact_trial;
    ind = cellfun(@(x) find(x == max(x),1),mean_amp, 'UniformOutput',false);
     
    criteria=cellfun(@(a,b,c,d) max(int8(a>mean_amp_thresh) + int8(b<p_val_thresh) +int8((c+d)<=5)), mean_amp,t_test,large_pre_mean, artifact_trial);
    stims.resp = criteria >=3;   
    criteria=cellfun(@(c,d) max(int8((c+d)<=5)),large_pre_mean, artifact_trial);
   
    stims.nonresp = criteria>=1;    
   stims.session = all_stims(i).name;
   all_vis_stims = [all_vis_stims,stims];
end
save(savefile, 'all_vis_stims');
end
