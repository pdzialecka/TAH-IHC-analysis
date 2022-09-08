function [t_results,t_stats] = extract_results(result_files,stats_files,result_type)
    %%
    
    %%
    result_files_ = result_files(contains({result_files.name}',result_type));
    stats_file_ = stats_files(contains({stats_files.name}',result_type));

    t = {};
    t_results = [];
    t_stats = [];
    var_names_h = {'Gamma vs Sham h','Theta vs Sham h','LTD vs Sham h'};
%     gap_row_table = array2table(nan(3,4),'VariableNames',cond_names);

    if ~isempty(result_files_)
        for i = 1:length(result_files_)
            t{i} = readtable(fullfile(result_files_(i).folder,result_files_(i).name),...
                'VariableNamingRule','preserve');

            t_results = [t_results; t{i}]; % gap_row_table];
        end


        t_stats = readtable(fullfile(stats_file_(1).folder,stats_file_(1).name),...
                'VariableNamingRule','preserve');
            
        % add a matrix with significance test
        h_matrix = t_stats{1:end,2:end}<0.05;
        h_table = array2table(uint16(h_matrix),'VariableNames',var_names_h);
        
        t_stats = [t_stats, h_table];
        
    end
    
end
