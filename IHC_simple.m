%% Simple IHC analysis for TI-AD-Hipp project
% @author: pdzialecka

%%
clc;
close all;
clear all;

%% Analysis folder
analysis_folder = 'C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC'; % ADJUST
addpath(genpath(analysis_folder));

%% Find data folder within cohort folder
% here, use example image inside Data subfolder of analysis folder
base_folder = fullfile(analysis_folder);
cohort_folder = 'Cohort_2';
data_folder = fullfile(base_folder,'Data',cohort_folder,'Images');

%% Deconvolve all files inside the data folder
files = dir(fullfile(data_folder,'**','*.svs'));
deconvolve_full(files);

cd(analysis_folder);

%% Pre-select all ROIs
% select magnification of choice for ROI images, e.g. 10 or 20x
magnification = 20;

% find all deconvolved files
files = dir(fullfile(data_folder,'**',strcat('*deconv.tif')));

% select an example file to test - comme
files = files(1);

for idx = 1:length(files)
    file_ = files(idx);
    select_roi(file_,magnification);
end

%% Analyse data to get % area covered by antibody
load_rois = 1; % load previously selected ROIs

% find all files
files = dir(fullfile(data_folder,'**',strcat('*deconv.tif')));

% select an example file to test
files = files(1);

analyse_data(files,load_rois);

