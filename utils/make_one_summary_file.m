function [] = make_one_summary_file(summary_folder)

    %%
    all_img_types = {'moc23','12f4','ct695','iba1','gfap','cfos','ki67','dcx','sox2'};
    
    cohort_subfolder = 'Cohorts_2-5_6mo';
    fname_1 = fullfile(summary_folder,cohort_subfolder,'IHC_density_results.xlsx');
    fname_2 = fullfile(summary_folder,cohort_subfolder,'IHC_count_results.xlsx');
    fname_3 = fullfile(summary_folder,cohort_subfolder,'IHC_ratio_results.xlsx');

%     cond_names = {'Sham','40 Hz','8 Hz','LTD'};
%         var_names = {'Sham vs Gamma','Sham vs Theta','Sham vs LTD'};
        
%         empty_table = cell2table({'' '' '' '';'' '' '' ''; '' '' '' ''},'VariableNames',cond_names);
%         gap_row_table = array2table(nan(3,4),'VariableNames',cond_names);

    
    for img_type_idx = 1:length(all_img_types)
        img_type = all_img_types{img_type_idx};
        
        result_files = dir(fullfile(summary_folder,cohort_subfolder,img_type,'Stats','*results.xlsx'));
        stats_files = dir(fullfile(summary_folder,cohort_subfolder,img_type,'Stats','*stats.xlsx'));
        
        
        % DENSITY RESULTS
        [t_results_1,t_stats_1] = extract_results(result_files,stats_files,'density');
        if ~isempty(t_results_1)
            writetable(t_results_1,fname_1,'Sheet',strcat(img_type,'_results'));
            writetable(t_stats_1,fname_1,'Sheet',strcat(img_type,'_stats'));
        end
        
        % CELL COUNT RESULTS
        [t_results_2,t_stats_2] = extract_results(result_files,stats_files,'count');
        if ~isempty(t_results_2)
            writetable(t_results_2,fname_2,'Sheet',strcat(img_type,'_results'));
            writetable(t_stats_2,fname_2,'Sheet',strcat(img_type,'_stats'));
        end
        
        % PERCENTAGE POSITIVE (CFOS)
        [t_results_3,t_stats_3] = extract_results(result_files,stats_files,'ratio');
        if ~isempty(t_results_3)
            writetable(t_results_3,fname_3,'Sheet',strcat(img_type,'_results'));
            writetable(t_stats_3,fname_3,'Sheet',strcat(img_type,'_stats'));
        end
        
    end
    
end
