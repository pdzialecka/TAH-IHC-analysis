function [p,h] = compute_stats(quantity_to_plot_all,quantity_name,img_type,...
                               stats_folder,roi_idxs,save_results)

    %% Compute and save stats
    % @author: pdzialecka

    % not all data normally distributed; use wilcoxon test instead of ttest

    %%
    if length(quantity_to_plot_all) == 6
        [roi_names,~,roi_no] = get_roi_list_IF();
    else
        [roi_names,~,roi_no] = get_roi_list();
    end
    var_names = {'Gamma vs Sham','Theta vs Sham','LTD vs Sham'};

    %%
    if ~exist('roi_idxs','var')
        roi_idxs = 1:roi_no;
    end
    
    if ~exist('save_results','var')
        save_results = 1;
    end

    %%
    p = nan(length(roi_idxs),3);
    h = nan(length(roi_idxs),3);
    
    %% Calculate stats
    for roi_idx = roi_idxs
        roi_results = quantity_to_plot_all{roi_idx};
    %     [h,p] = kstest(results_density); % h = 0 if normally distributed

        % compare all stims to sham
        [p(roi_idx,1),h(roi_idx,1)] = ranksum(roi_results(:,1),roi_results(:,2));
        [p(roi_idx,2),h(roi_idx,2)] = ranksum(roi_results(:,1),roi_results(:,3));
        [p(roi_idx,3),h(roi_idx,3)] = ranksum(roi_results(:,1),roi_results(:,4));
    end
    
    %% Save results
    if save_results
        file_name = sprintf('%s_%s_stats',img_type,quantity_name);
        save(fullfile(stats_folder,strcat(file_name,'.mat')),'p','h');

        stats_T = array2table(p,'VariableNames',var_names,'RowNames',roi_names(roi_idxs));
        table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
        writetable(stats_T,table_name,'WriteRowNames',true);
    end
    
end
