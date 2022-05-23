%% IHC analysis pipeline for TI-AD-Hipp project
% @author: pdzialecka

%% File directories
cohort_folder = 'Cohort_1_AD_Hipp18-23';

main_folder = 'C:\Users\Pat\Desktop';
data_folder = fullfile(main_folder,cohort_folder,'IHC\Images');
processed_folder = fullfile(data_folder,'Processed');

%% Analysis steps
deconvolve = 1;
select_rois = 1;
analyse_moc = 1;

%% Deconvolve all files inside a folder
if deconvolve
    files = dir(fullfile(data_folder,'*.svs'));
    deconvolve_full(files);
end

%% Pre-select all ROIs
if select_rois
    magnification = 10;
    files = dir(fullfile(processed_folder,strcat('*deconv.tif')));

    for idx = 1:length(files)
        file_ = files(idx);
        select_roi(file_,magnification)
    end
end

%% Analyse moc23 data
load_rois = 1;

if analyse_moc
    image_type = 'moc23';
    files = dir(fullfile(processed_folder,strcat('*',image_type,'*deconv.tif')));
end

