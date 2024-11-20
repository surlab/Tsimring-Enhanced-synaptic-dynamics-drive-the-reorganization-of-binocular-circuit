
% find pair-wise offset between orientation preferences
% for each spine in fov
% input = array of orientation preferences 
% output = matrix of offsets between each pair of spines
function offset = get_orientation_offset(oris)

offset = zeros(length(oris));
for i = 1:length(oris)
    for ii = 1:length(oris)
        diff1 = abs(oris(i)- oris(ii));
        if diff1 > 90
            diff2 = oris(i) + oris(ii);
            diff1 = abs(180 - diff2);
        end
        
        offset(i,ii) = diff1;
    end
end