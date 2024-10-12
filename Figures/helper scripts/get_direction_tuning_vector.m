function [DSI, Dir_pref, R_dir] = get_direction_tuning_vector(oris_rad, mean_amplitude)
    R_dir = mean_amplitude*exp(oris_rad'*1i);
    Dir_pref = mod(atan2(imag(R_dir),real(R_dir)),2*pi)*(180/pi);
    DSI = abs(R_dir)./sum(mean_amplitude,2);
end