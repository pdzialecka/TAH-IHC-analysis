%% IHC analysis pipeline for TI-AD-Hipp project
% @author: pdzialecka

%% Deconvolve all files inside a folder
main_folder = 'C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC';
folder = fullfile(main_folder,'images');

files = dir(fullfile(folder,'*.svs'));

deconvolve_full(files);

%%
