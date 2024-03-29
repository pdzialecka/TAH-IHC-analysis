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
base_folder = 'D:\TAH';

%% Analysis steps
deconvolve = 0;
select_rois = 0;
create_mask = 0;
analyse = 0;
analyse_all = 0;
summarise = 1;
summarise_all = 0;
summarise_tables = 1;

%% Analysis settings
all_img_types = {'moc23','12f4','ct695','iba1','gfap','cfos','ki67','dcx','sox2'};
img_type = 'iba1'; % specific analysis

%% Cohort case
cohort_case = 2; % 1 = 13mo (cohort 1), 2 = 6mo (cohorts 2-5)

if cohort_case == 1
    cohort_folders = {'Cohort_1'};
    
elseif cohort_case == 2
    cohort_folders = {'Cohort_2','Cohort_3','Cohort_4','Cohort_5','Cohort_6'};
end

%% Within animal analysis inside cohort folders
for cohort_idx = 1:length(cohort_folders)
    
    %% File directories
    cohort_folder = cohort_folders{cohort_idx};
    data_folder = fullfile(base_folder,'Data',cohort_folder,'IHC\Images');
    processed_folder = fullfile(base_folder,'Data',cohort_folder,'IHC\Images_processed');
%     roi_images_folder = fullfile(base_folder,'Data',cohort_folder,'IHC\ROI_images');
    
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
            
            % or select all ROIs per mouse
            mouse_pfolders = dir(processed_folder);
            mouse_pfolders(ismember({mouse_pfolders.name}, {'.', '..'})) = [];
            
            for m_idx = 1:length(mouse_pfolders)
                mouse_pfolder = fullfile(mouse_pfolders(m_idx).folder,mouse_pfolders(m_idx).name);
                files = dir(fullfile(mouse_pfolder,strcat('*deconv.tif')));
                
                % rearrange files so moc23 is first
                files = rearrange_files(files,all_img_types);
                
                for idx = 1:length(files)
                    file_ = files(idx);
                    select_roi_semi(file_);
%                     select_roi_semi(file_,0); % to enforce manual selection
                end
            end
            
        end
        
        %% Create slice mask (run before + after artefact removal)
        if create_mask
            
            % mask made of composite image
            all_files = dir(fullfile(processed_folder,'*',strcat('*.tif')));
            keep_idxs = ~contains({all_files.name}','deconv');
            files = all_files(keep_idxs);
            
            create_slice_masks(files);
        end

        %% Analyse all data or that from one antibody only
        if analyse
            load_rois = 1;

            if analyse_all
                files = dir(fullfile(processed_folder,'*',strcat('*deconv.tif')));
            else
                files = dir(fullfile(processed_folder,'*',strcat('*',img_type,'*deconv.tif')));
            end

%             files = files(1); % test file
            analyse_data(files,load_rois);

        end
    end

end


%% Summarise results between animals within age groups
mice_to_exclude = {'AD-Hipp41','AD-Hipp44','AD-Hipp45','AD-Hipp53'};

if summarise
    results_folder = fullfile(base_folder,'IHC_results');
    save_cohort_info(results_folder);

    if summarise_all
        for i = 1:length(all_img_types)
            img_type_i = all_img_types{i};
            summarise_results(base_folder,cohort_case,img_type_i,mice_to_exclude);
        end

    else
        summarise_results(base_folder,cohort_case,img_type,mice_to_exclude);
    end
end

%% Create one summary table
if summarise_tables
    summary_folder = fullfile(base_folder,'IHC_results');
    make_one_summary_file(summary_folder);
end
