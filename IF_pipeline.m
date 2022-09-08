%% IF analysis pipeline for TI-AD-Hipp project
% @author: pdzialecka

%%
clc;
close all;
clear all;

%% Analysis folder
analysis_folder = 'C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC';
addpath(genpath(analysis_folder));

%% Root data folder
base_folder = 'C:\Users\Pat\Desktop\TAH';
% base_folder = 'K:\TAH';

%% Analysis steps
analyse = 1;
summarise = 1;

%% Cohort case
cohort_case = 2; % 1 = 13mo (cohort 1), 2 = 6mo (cohorts 2-5)

if cohort_case == 1
    cohort_folders = {'Cohort_1'};
    
elseif cohort_case == 2
    cohort_folders = {'Cohort_2','Cohort_3','Cohort_4','Cohort_5'};
end


%% Within animal analysis inside cohort folders
for cohort_idx = 1:length(cohort_folders)
    
    %% File directories
    cohort_folder = cohort_folders{cohort_idx};
    data_folder = fullfile(base_folder,'Data',cohort_folder,'IF\Images');
    roi_images_folder = fullfile(base_folder,'Data',cohort_folder,'IF\ROI_images');
    
    %%
    if exist(roi_images_folder,'dir')
        %%
        if analyse
            files = dir(fullfile(roi_images_folder,'*','*','*.tif'));
            files(contains({files.name}','Merged')) = [];
            analyse_data_IF(files);
        end
    end

end

%% Summarise results between animals
if summarise
    results_folder = fullfile(base_folder,'IF_results');
    save_cohort_info(results_folder);

    summarise_results_IF(base_folder,cohort_case);
end
