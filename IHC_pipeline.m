%% IHC analysis pipeline for TI-AD-Hipp project
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

%% Analysis steps
deconvolve = 1;
select_rois = 1;
analyse = 1;
analyse_all = 1;

all_image_types = {'moc23','cfos','GFAP','Iba1','12f4','ki67','ct695'};
image_type = 'moc23'; % specific analysis
summarise = 1;

%% Cohort case
cohort_case = 1; % 1 = 13mo (cohort 1), 2 = 6mo (cohorts 2-5)

if cohort_case == 1
    cohort_folders = {'Cohort_1'};
elseif cohort_case == 2
    cohort_folders = {'Cohort_2','Cohort_3','Cohort_4','Cohort_5'};
end


%% Within animal analysis inside cohort folders
for cohort_idx = 1:length(cohort_folders)
    %% File directories
    cohort_folder = cohort_folders{cohort_idx};
    data_folder = fullfile(base_folder,'Data',cohort_folder,'IHC\Images');
    processed_folder = fullfile(base_folder,'Data',cohort_folder,'IHC\Images_processed');
    
    %%
    if exist(data_folder)
        %% Deconvolve all files inside a folder
        if deconvolve
            files = dir(fullfile(data_folder,'*','*.svs'));
            deconvolve_full(files);

            cd(analysis_folder);
        end

        %% Pre-select all ROIs
        if select_rois
            magnification = 20;

            if analyse_all
                files = dir(fullfile(processed_folder,'**',strcat('*deconv.tif')));
            else
                files = dir(fullfile(processed_folder,'**',strcat('*',image_type,'*deconv.tif')));
            end

%             files = files(1); % test file

            for idx = 1:length(files)
                file_ = files(idx);
                select_roi(file_,magnification);
            end
        end

        %% Analyse all data or that from one antibody only
        if analyse
            load_rois = 1;

            if analyse_all
                files = dir(fullfile(processed_folder,'**',strcat('*deconv.tif')));
            else
                files = dir(fullfile(processed_folder,'**',strcat('*',image_type,'*deconv.tif')));
            end

%             files = files(1); % test file
            analyse_data(files,load_rois);

        end
    end

end


%% Summarise results between animals within age groups
% TODO: test this works for cohort_case == 2

if summarise
    results_folder = fullfile(base_folder,'IHC_results');
    save_cohort_info(results_folder);

    if analyse_all
        for i = 1:length(all_image_types)
            image_type_i = all_image_types{i};
            summarise_results(base_folder,cohort_case,image_type_i);
        end

    else
        summarise_results(base_folder,cohort_case,image_type);
    end
end
