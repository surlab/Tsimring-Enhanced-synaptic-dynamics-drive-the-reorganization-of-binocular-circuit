function ori_pref = preferred_orientation(lin_angles, firing_rates)
    oris_rad = (lin_angles'.*pi)./180; % column vector
    num_oris = length(oris_rad)/2;
    fr_half = (firing_rates(:,1:num_oris) + firing_rates(:,num_oris+1:end))/2;
    oris_rad_half = oris_rad(1:num_oris); % oris_rad_half = oris_rad(num_oris+1:end);
    R_ori = fr_half*exp(oris_rad_half*2*1i);
    ori_pref = mod(atan2(imag(R_ori),real(R_ori)),2*pi)*(90/pi);
end