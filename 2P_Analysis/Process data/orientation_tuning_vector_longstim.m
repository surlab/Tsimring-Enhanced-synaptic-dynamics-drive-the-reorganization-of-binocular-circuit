function orientation_tuning_vector_longstim(type,mean_data, std_data,unique_orientations, save_path,Day)
% Orientation_analysis 
if contains(save_path, 'contra')
    c = 'r';
elseif contains(save_path, 'ipsi')
    c = 'b';
else
    c = 'g';
end

%only look at 8 directions
if str2double(Day(2:end)) > 50
    index = 1:2:16;
else
    index = 1:8;
end

unique_orientations = unique_orientations(index);


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

% find DSI and direction preference
R_dir = mean_data_responsive*exp(oris_rad*1i);
Dir_pref = mod(atan2(imag(R_dir),real(R_dir)),2*pi)*(180/pi);
DSI = abs(R_dir)./sum(mean_data_responsive,2);
    

% for i = 1:size(mean_data_responsive,1)
%     %i = responsive_cell_index(k);
%     figure;
%     subplot(1,2,1)
%     polarplot([oris_rad_half*2;oris_rad_half(1)*2],[mean_data_half(i,:), mean_data_half(i,1)],c,'LineWidth',2);
%     hold on,polarplot([0,mod(atan2(imag(R_ori(i)),real(R_ori(i))),2*pi)], [0,OSI(i)], 'k','LineWidth',1);
%     thetaticks([0:90:360])
%     thetaticklabels([0:45:180])
%     title(['OTI: ', num2str(OSI(i)), ' Orientation pref: ', num2str(Ori_pref(i))])
% 
%     subplot(1,2,2)
%     polarplot([oris_rad;oris_rad(1)],[mean_data_responsive(i,:), mean_data_responsive(i,1)],c,'LineWidth',2);
%     hold on,polarplot([0,mod(atan2(imag(R_dir(i)),real(R_dir(i))),2*pi)], [0,DSI(i)], 'k','LineWidth',1);
%     title(['DSI: ', num2str(DSI(i)), ' Direction pref: ', num2str(Dir_pref(i))])
%     thetaticks([0:45:360])
%     thetaticklabels([0:45:360])
%     
%     if ~isfolder(fullfile(save_path, 'polarplots'))
%         mkdir(fullfile(save_path, 'polarplots'))
%     end
% 
%     
% 
%     savefig(fullfile(save_path, 'polarplots',['ROI',num2str(i), '.fig']))
%     close all
% end
oris_temp = oris_deg_psycopy;
oris_temp(oris_deg_psycopy<0) = oris_deg_psycopy(oris_deg_psycopy<0)+360;
[oris_deg_sort,inds] = sort(oris_temp);
% for i = 1:size(mean_data_responsive,1)
%     %i = responsive_cell_index(k);
%     figure;
%     
%     errorbar(oris_deg_sort,mean_data_responsive(i,inds), std_data(i,inds)/sqrt(size(std_data,2)),c,'LineWidth',2);
%     title(['ROI: ', num2str(i)]);
%     ylabel('Mean trial amplitude (Z)')
%     xlabel('Degrees');
%      if ~isfolder(fullfile(save_path, 'ori_figs'))
%         mkdir(fullfile(save_path, 'ori_figs'))
%     end
% 
%     savefig(fullfile(save_path, 'ori_figs',['ROI',num2str(i),'.fig']))
%     close all
% 
% end    
%    
    
save(fullfile(save_path, [type,'.mat']), 'OSI', 'DSI', 'Ori_pref', 'Dir_pref', 'mean_data_responsive', 'std_data' )


end
