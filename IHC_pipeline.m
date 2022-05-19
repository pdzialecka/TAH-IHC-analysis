%% IHC analysis pipeline for TI-AD-Hipp project
% @author: pdzialecka

%% Deconvolve all files inside a folder
cohort_folder = 'Cohort_2_AD-Hipp28_29_32_35';

main_folder = 'C:\Users\Pat\Desktop';
folder = fullfile(main_folder,cohort_folder,'IHC\Images');
files = dir(fullfile(folder,'*.svs'));

deconvolve_full(files);

%%
