%% calcium rate
function [ca_rate, smooth_data, pks, lks] = get_ca_rate(data, num_time,alpha, std_thresh, conseq_frames)
    smooth_data = exp_weight_filter(data, alpha);
    [pks,lks] = findpeaks(smooth_data, 'MinPeakHeight',std_thresh, 'MinPeakDistance',conseq_frames);  
    ca_rate = length(pks)/num_time;
end
