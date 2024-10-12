function smooth_data = exp_weight_filter(data, alpha)
    smooth_data = zeros(size(data));
    smooth_data(1) = data(1);
    for t = 2:length(data)
        smooth_data(t) = alpha*data(t) + (1-alpha)*smooth_data(t-1);
    end
end