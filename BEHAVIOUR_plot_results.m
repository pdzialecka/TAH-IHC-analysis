%% Behavioural results plotting for TI-AD-Hipp project
% @author: pdzialecka

% This function only extracts presaved results, plots them & computes stats
% consistently with IHC plots
% For full behavioural analysis, check Behaviour analysis functions

%%
clc;
close all;
clear all;

%% Analysis folder
analysis_folder = 'C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC';
addpath(genpath(analysis_folder));

%% Define root folders
root_data_folder = 'D:\TAH\Behaviour_results\from_Maria\Results\';
root_results_folder = 'D:\TAH\Behaviour_results\Cohorts_2-5_6mo';

%% Settings
% cond_names = {'Sham','Delta','Theta','Gamma'};

save_results = 1;
roi_idxs = 1;
img_type = 'Behaviour';

%% Load cohort info
% cohort_idxs = 2:5;
% cohort_info_folder = fullfile('D:\TAH\IHC_results','Cohort_info');
% cohort_infos = {};
% 
% mouse_ids = [];
% mouse_cond_idxs = [];
% 
% for i = 1:length(cohort_idxs)
%     cohort_infos{i} = load(fullfile(cohort_info_folder,...
%         sprintf('Cohort_%d_info.mat',cohort_idxs(i)))).cohort;
%     mouse_ids = [mouse_ids cohort_infos{i}.mouse_ids];
%     mouse_cond_idxs = [mouse_cond_idxs cohort_infos{i}.mouse_cond_idxs];
% end
% 
% mouse_names = mouse_ids_to_names(mouse_ids);
% [~,m_idxs] = sort(mouse_names);
% mouse_names = mouse_names(m_idxs);
% mouse_cond_idxs = mouse_cond_idxs(m_idxs);
% mouse_no = length(mouse_ids);

%% %%%%%%%%%%%%%%%%%%%%% OLM %%%%%%%%%%%%%%%%%%%%%
%% Load data
load(fullfile(root_data_folder,'OLM\OLM_data.mat'));

%% Extract treatment groups
% same as above but safer to keep this
treatment_groups = dataset_OLM_data.Treatment;
mouse_cond_idxs = nan(size(treatment_groups))';

mouse_cond_idxs(treatment_groups == 'Sham') = 1;
mouse_cond_idxs(treatment_groups == 'LTD') = 4;
mouse_cond_idxs(treatment_groups == '8 Hz') = 3;
mouse_cond_idxs(treatment_groups == '40 Hz') = 2;

%% Specify results folders
results_folder = fullfile(root_results_folder,'OLM');
if ~exist(results_folder)
    mkdir(results_folder)
end  

stats_folder = fullfile(results_folder,'Stats');

%% Extract data
Training_expTime = dataset_OLM_data.Total_expTime_training';
Training_DI = dataset_OLM_data.Training_DI';
Test_DI = dataset_OLM_data.Test_DI';

%% Exploration time
[results_expTime] = plot_results(Training_expTime,'expTime',mouse_cond_idxs,...
                        img_type,results_folder,roi_idxs,save_results);

[p_expTime,h_expTime] = compute_stats_behaviour(results_expTime,'expTime',img_type,...
                               stats_folder,save_results);

%% Training DI
[results_training_DI] = plot_results(Training_DI,'Training_DI',mouse_cond_idxs,...
                        img_type,results_folder,roi_idxs,save_results);

[p_training_DI,h_training_DI] = compute_stats_behaviour(results_training_DI,'Training_DI',img_type,...
                               stats_folder,save_results);
                           
%% Test DI
[results_test_DI] = plot_results(Test_DI,'Test_DI',mouse_cond_idxs,...
                        img_type,results_folder,roi_idxs,save_results);

[p_test_DI,h_test_DI] = compute_stats_behaviour(results_test_DI,'Test_DI',img_type,...
                               stats_folder,save_results);
                       

%% %%%%%%%%%%%%%%%%%%%%% Y-maze %%%%%%%%%%%%%%%%%%%%%
%% Load data
load(fullfile(root_data_folder,'Y_maze\YMaze_data.mat'));

%% Specify results folders
results_folder = fullfile(root_results_folder,'Y_Maze');
if ~exist(results_folder)
    mkdir(results_folder)
end

stats_folder = fullfile(results_folder,'Stats');

%% Extract data
SAP_index = [YMaze_data.SAP_index];
AAR_index = [YMaze_data.AAR_index];
SAR_index = [YMaze_data.SAR_index];

%% SAP index
[results_SAP_idx] = plot_results(SAP_index,'SAP_index',mouse_cond_idxs,...
                        img_type,results_folder,roi_idxs,save_results);

[p_SAP_idx,h_SAP_idx] = compute_stats_behaviour(results_SAP_idx,'SAP_index',img_type,...
                               stats_folder,save_results);

%% AAR index
[results_AAR_idx] = plot_results(AAR_index,'AAR_index',mouse_cond_idxs,...
                        img_type,results_folder,roi_idxs,save_results);

[p_AAR_idx,h_AAR_idx] = compute_stats_behaviour(results_AAR_idx,'AAR_index',img_type,...
                               stats_folder,save_results);
                           
%% SAR index
[results_SAR_idx] = plot_results(SAR_index,'SAR_index',mouse_cond_idxs,...
                        img_type,results_folder,roi_idxs,save_results);

[p_SAR_idx,h_SAR_idx] = compute_stats_behaviour(results_SAR_idx,'SAR_index',img_type,...
                               stats_folder,save_results);
