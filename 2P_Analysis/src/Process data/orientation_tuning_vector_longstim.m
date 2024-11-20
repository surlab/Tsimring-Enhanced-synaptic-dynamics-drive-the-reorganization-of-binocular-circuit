function orientation_tuning_vector_longstim(type,mean_data, std_data,unique_orientations, save_path)
% Orientation_analysis 

oris_deg_psycopy = (180-unique_orientations); %% based on psycho conversion
oris_rad = (oris_deg_psycopy.*pi)./180;
mean_data_responsive = mean_data;
mean_data_responsive(mean_data_responsive < 0) = 0;
num_oris = length(oris_rad)/2
mean_data_half = [mean_data_responsive(:,1:num_oris) + mean_data_responsive(:,num_oris+1:end)]/2;

oris_rad_half = oris_rad(num_oris+1:end);
R_ori = mean_data_half*exp(oris_rad_half*2*1i);
Ori_pref = mod(atan2(imag(R_ori),real(R_ori)),2*pi)*(90/pi);
OSI = abs(R_ori)./sum(mean_data_half,2);

oris_temp = oris_deg_psycopy;
oris_temp(oris_deg_psycopy<0) = oris_deg_psycopy(oris_deg_psycopy<0)+360;
[oris_deg_sort,inds] = sort(oris_temp);

    
save(fullfile(save_path, [type,'.mat']), 'OSI', 'DSI', 'Ori_pref', 'Dir_pref', 'mean_data_responsive', 'std_data' )


end
