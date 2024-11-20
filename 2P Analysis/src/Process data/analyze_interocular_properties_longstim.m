function analyze_interocular_properties_longstim(inputpath,FOV_name,savename, savepath)


%% load vis stim files 
contra_file = fullfile(inputpath, [FOV_name, '_contra'], 'ori_analysis.mat');
ipsi_file = fullfile(inputpath, [FOV_name, '_ipsi'], 'ori_analysis.mat');
binoc_file = fullfile(inputpath, [FOV_name, '_binoc'], 'ori_analysis.mat');
if exist(contra_file,"file") & exist(ipsi_file, "file") & exist(binoc_file, "file")
    contra =  load(contra_file);
    ipsi = load(ipsi_file);
    binoc = load(binoc_file);


%% preallocate matrices
mean_data_resp = NaN([size(contra.mean_data_responsive),3]);
mean_data_resp(:,:,1) = contra.mean_data_responsive;
mean_data_resp(:,:,2) = ipsi.mean_data_responsive;
mean_data_resp(:,:,3) = binoc.mean_data_responsive;

Ori_resp = NaN([size(contra.Ori_pref),3]);
Ori_resp(:,:,1) = contra.Ori_pref;
Ori_resp(:,:,2) = ipsi.Ori_pref;
Ori_resp(:,:,3) = binoc.Ori_pref;

%% orientation difference 
diff = (abs(Ori_resp(:,:,1) - Ori_resp(:,:,3)));
diff2 = (abs(Ori_resp(:,:,1) + Ori_resp(:,:,3)));
C_B_ori_diff = diff;
if ~isempty(find(diff  > 90))
    C_B_ori_diff(diff > 90) =  abs(180 - diff2(diff > 90));
end

diff = (abs(Ori_resp(:,:,2) - Ori_resp(:,:,3)));
diff2 = (abs(Ori_resp(:,:,2) + Ori_resp(:,:,3)));
B_I_ori_diff = diff;
if ~isempty(find(diff> 90))
    B_I_ori_diff(diff > 90) = abs(180 - diff2(diff > 90));
end

diff = (abs(Ori_resp(:,:,2) - Ori_resp(:,:,1)));
diff2 = (abs(Ori_resp(:,:,2) + Ori_resp(:,:,1)));
C_I_ori_diff = diff;
if ~isempty(find(diff>90))
    C_I_ori_diff(diff > 90) = abs(180 - diff2(diff > 90));
end
%% Direction congruency
[R,pval] = corr(mean_data_resp(:,:,1)',mean_data_resp(:,:,3)');
C_B_dir_cong = diag(R);

[R,pval] = corr(mean_data_resp(:,:,1)',mean_data_resp(:,:,2)');
C_I_dir_cong = diag(R);

[R,pval] = corr(mean_data_resp(:,:,3)',mean_data_resp(:,:,2)');
B_I_dir_cong = diag(R);

%% ODI and BRM calc

ipsi_resp = max(mean_data_resp(:,:,2),[],2);
contra_resp = max(mean_data_resp(:,:,1),[],2);
binoc_resp = max(mean_data_resp(:,:,3),[],2);
ODI = (contra_resp - ipsi_resp)./(contra_resp + ipsi_resp);

BRM = (binoc_resp - max(ipsi_resp, contra_resp))./(binoc_resp + max(ipsi_resp, contra_resp))


%% save ODI
if ~exist(savepath, 'dir')
   mkdir(savepath)
end
save(fullfile(savepath, [savename, '.mat']), 'ODI','BRM','C_B_dir_cong','C_I_dir_cong','B_I_dir_cong','C_B_ori_diff','C_I_ori_diff','B_I_ori_diff', 'mean_data_resp')
end