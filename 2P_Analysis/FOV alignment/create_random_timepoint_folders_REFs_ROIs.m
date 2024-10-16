%% Renaming Chronic imaging FOV alignment

clear all

%path = '/Users/ktsimring/Dropbox (MIT)/Sur Lab/Projects/Development project/Binocular Matching/Spine imaging/';
path = 'C:/Users/Katya-PC/Dropbox (MIT)/Sur Lab/Projects/Development project/Binocular Matching/Spine imaging/';
input_path = 'D:BinocularMatching\Spines\2P_imaging';
savepath = fullfile(path, 'Analyzed Data');
load(fullfile(savepath,'all_analyzed_spines_BM014-BM021.mat'));
load(fullfile(savepath,'visually_responsive_spines_BM014-BM021.mat'));
savefile = fullfile(savepath, 'soma_spine_all_alignment_all_resp_BM014-BM021.mat')


%%
stim_binoc = all_stims(strcmp({all_stims.name},'binoc')==1);
stim_contra = all_stims(strcmp({all_stims.name}, 'contra')==1);
stim_ipsi = all_stims(strcmp({all_stims.name}, 'ipsi')==1);
stim_binoc = rmfield(stim_binoc,'name');
stim_ipsi = rmfield(stim_ipsi,'name');
stim_contra = rmfield(stim_contra, 'name');

vis_stim_binoc = all_vis_stims(strcmp({all_vis_stims.session},'binoc')==1);
vis_stim_contra = all_vis_stims(strcmp({all_vis_stims.session}, 'contra')==1);
vis_stim_ipsi = all_vis_stims(strcmp({all_vis_stims.session}, 'ipsi')==1);
vis_stim_binoc = rmfield(vis_stim_binoc,'session');
vis_stim_ipsi = rmfield(vis_stim_ipsi,'session');
vis_stim_contra = rmfield(vis_stim_contra, 'session');

stim_table_binoc = struct2table(stim_binoc);
stim_table_contra = struct2table(stim_contra);
stim_table_ipsi = struct2table(stim_ipsi);

vis_stim_table_binoc = struct2table(vis_stim_binoc);
vis_stim_table_contra = struct2table(vis_stim_contra);
vis_stim_table_ipsi = struct2table(vis_stim_ipsi);
%%
combine_stim_binoc = [vis_stim_table_binoc,stim_table_binoc];
combine_stim_contra = [vis_stim_table_contra, stim_table_contra];
combine_stim_ipsi = [vis_stim_table_ipsi, stim_table_ipsi];