function [p,h] = compute_stats_behaviour(quantity_to_plot_all,quantity_name,img_type,...
                               stats_folder,save_results)

    %% Compute and save stats for behavioural results
    % @author: pdzialecka

    % not all data normally distributed; use wilcoxon test instead of ttest

    %%
    var_names = {'Sham','Delta','Theta','Gamma'};
    row_names = {'Cond vs sham','Cond vs chance'};

    %%    
    if ~exist('save_results','var')
        save_results = 1;
    end

    %%
    p = nan(2,4);
    h = nan(2,4);
    
    %% Calculate stats: compare each to sham
    roi_idx = 1;
    roi_results = quantity_to_plot_all{roi_idx};
%     [h,p] = kstest(results_density); % h = 0 if normally distributed

    % compare all stims to sham
    [p(roi_idx,1),h(roi_idx,1)] = ranksum(roi_results(:,1),roi_results(:,1));
    [p(roi_idx,2),h(roi_idx,2)] = ranksum(roi_results(:,1),roi_results(:,2));
    [p(roi_idx,3),h(roi_idx,3)] = ranksum(roi_results(:,1),roi_results(:,3));
    [p(roi_idx,4),h(roi_idx,4)] = ranksum(roi_results(:,1),roi_results(:,4));
    
    %% Calculate stats: compare each to chance
    chance_level = nan;
    
    if contains(quantity_name,'DI')
        chance_level = 0.5;
    elseif contains(quantity_name,'SAP')
        chance_level = 6/27*100;
    end
    
    if ~isnan(chance_level)
        [p(2,:),h(2,:)] = ranksum_matrix(roi_results,chance_level);
    end
    
    %% Save results
    if save_results
        file_name = sprintf('%s_%s_stats',img_type,quantity_name);
        save(fullfile(stats_folder,strcat(file_name,'.mat')),'p','h');

        stats_T = array2table(p,'VariableNames',var_names,'RowNames',row_names);
        table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
        writetable(stats_T,table_name,'WriteRowNames',true);
    end
    
end
