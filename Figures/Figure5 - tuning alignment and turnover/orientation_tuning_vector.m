function [Ori_pref,OSI, Dir_pref, DSI] = orientation_tuning_vector(all_mean_data, unique_orientations)
    num_dir = length(unique_orientations);
    num_oris = num_dir/2;
    oris_deg_psycopy = (180-unique_orientations); %% based on psycho conversion
    oris_rad = (oris_deg_psycopy.*pi)./180;
    oris_rad_half = oris_rad(num_oris+1:end);
    all_mean_data_half = [all_mean_data(:,1:num_oris) + all_mean_data(:,num_oris+1:end)]/2;
       
    % find bootstrapped OSI and Ori pref
    R_ori = all_mean_data_half*exp(oris_rad_half'*2*1i);
    Ori_pref = mod(atan2(imag(R_ori),real(R_ori)),2*pi)*(90/pi);
    OSI = abs(R_ori)./sum(all_mean_data_half,2);
   
    % find bootstrapped DSI and direction preference
    R_dir = all_mean_data*exp(oris_rad'*1i);
    Dir_pref = mod(atan2(imag(R_dir),real(R_dir)),2*pi)*(180/pi);
    DSI = abs(R_dir)./sum(all_mean_data,2);

end