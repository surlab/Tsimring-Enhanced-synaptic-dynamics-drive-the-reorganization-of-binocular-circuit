function [OSI, Ori_pref, R_ori, mean_data_half, oris_rad_half] = get_orientation_tuning_vector(oris_rad, mean_amplitude)
    num_oris = length(oris_rad)/2
    mean_data_half = [mean_amplitude(:,1:num_oris) + mean_amplitude(:,num_oris+1:end)]/2;
    oris_rad_half = oris_rad(num_oris+1:end);
    R_ori = mean_data_half*exp(oris_rad_half'*2*1i);
    Ori_pref = mod(atan2(imag(R_ori),real(R_ori)),2*pi)*(90/pi);
    OSI = abs(R_ori)/sum(mean_data_half);
end