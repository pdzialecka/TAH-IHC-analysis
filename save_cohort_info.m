function [] = save_cohort_info(results_folder)
    %% Save cohort info for easy retrieval later
    % @author: pdzialecka
    
    %% Cohort info subfolder
    cohort_info_folder = fullfile(results_folder,'Cohort_info');
    if ~exist(cohort_info_folder)
        mkdir(cohort_info_folder);
    end
    
    %% Save info
    % cohort 1
    cohort = [];
    cohort.cohort_id = 1;
    cohort.cohort_age = 13; % months
    cohort.conds = {'Sham','40 Hz', '8 Hz', 'LTD'};
    cohort.mouse_ids = [18,19,20,21,22,23];
    cohort.mouse_names = mouse_ids_to_names(cohort.mouse_ids);
    cohort.mouse_sex = {'F','F','F','F','F','F'};
    cohort.mouse_cond_idxs = [2,1,1,2,1,2];
    cohort.mouse_conds = mouse_ids_to_conds(cohort.mouse_cond_idxs);
    save(fullfile(cohort_info_folder,'cohort_1_info.mat'),'cohort');

    % cohorts 2
    cohort = [];    
    cohort.cohort_id = 2;
    cohort.cohort_age = 6;
    cohort.conds = {'Sham','40 Hz','8 Hz','LTD'};
    cohort.mouse_ids = [29,32,28,35];
    cohort.mouse_names = mouse_ids_to_names(cohort.mouse_ids);
    cohort.mouse_sex = {'F','M','F','M'};
    cohort.mouse_cond_idxs = [1,2,2,3];
    cohort.mouse_conds = mouse_ids_to_conds(cohort.mouse_cond_idxs);
    save(fullfile(cohort_info_folder,'cohort_2_info.mat'),'cohort');


    % cohort 3
    cohort = [];
    cohort.cohort_id = 3;
    cohort.cohort_age = 6;
    cohort.conds = {'Sham','40 Hz','8 Hz','LTD'};
    cohort.mouse_ids = [38,33,37,34,39,36,40];
    cohort.mouse_names = mouse_ids_to_names(cohort.mouse_ids);
    cohort.mouse_sex = {'F','M','F','M','F','M','M'};
    cohort.mouse_cond_idxs = [1,3,3,1,2,1,2];
    cohort.mouse_conds = mouse_ids_to_conds(cohort.mouse_cond_idxs);
    save(fullfile(cohort_info_folder,'cohort_3_info.mat'),'cohort');
    
    % cohort 4
    cohort = [];
    cohort.cohort_id = 4;
    cohort.cohort_age = 6;
    cohort.conds = {'Sham','40 Hz','8 Hz','LTD'};
    cohort.mouse_ids = [];
    cohort.mouse_names = mouse_ids_to_names(cohort.mouse_ids);
    cohort.mouse_sex = {};
    cohort.mouse_cond_idxs = [];
    cohort.mouse_conds = mouse_ids_to_conds(cohort.mouse_cond_idxs);
    save(fullfile(cohort_info_folder,'cohort_4_info.mat'),'cohort');
    
    % cohort 5
    cohort = [];
    cohort.cohort_id = 5;
    cohort.cohort_age = 6;
    cohort.conds = {'Sham','40 Hz','8 Hz','LTD'};
    cohort.mouse_ids = [];
    cohort.mouse_names = mouse_ids_to_names(cohort.mouse_ids);
    cohort.mouse_sex = {};
    cohort.mouse_cond_idxs = [];
    cohort.mouse_conds = mouse_ids_to_conds(cohort.mouse_cond_idxs);
    save(fullfile(cohort_info_folder,'cohort_5_info.mat'),'cohort');

end
