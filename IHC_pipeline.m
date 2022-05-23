%% IHC analysis pipeline for TI-AD-Hipp project
% @author: pdzialecka

%%
clc;
close all;
clear all;

%% Analysis folder
analysis_folder = 'C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC';
addpath(genpath(analysis_folder));

%% File directories
cohort_folder = 'Cohort_1_AD_Hipp18-23';

main_folder = 'C:\Users\Pat\Desktop';
data_folder = fullfile(main_folder,cohort_folder,'IHC\Images');
% processed_folder = fullfile(data_folder,'Processed');

%% Analysis steps
deconvolve = 1;
select_rois = 1;
analyse = 1;
analyse_all = 1;

%% Deconvolve all files inside a folder
if deconvolve
    files = dir(fullfile(data_folder,'**','*.svs'));
    deconvolve_full(files);
    
    cd(analysis_folder);
end

%% Pre-select all ROIs
if select_rois
    magnification = 20;
    files = dir(fullfile(data_folder,'**',strcat('*deconv.tif')));
    
    files = files(3); % test file

    for idx = 1:length(files)
        file_ = files(idx);
        select_roi(file_,magnification);
    end
end

%% Analyse all data or that from one antibody only
if analyse
    load_rois = 1;
    image_type = 'moc23'; % specific analysis

    if analyse_all
        files = dir(fullfile(data_folder,'**',strcat('*deconv.tif')));

    else
        files = dir(fullfile(data_folder,'**',strcat('*',image_type,'*deconv.tif')));
    end
    
    files = files(3); % test file
    analyse_data(files,load_rois);
    
end
