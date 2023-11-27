function [p,h,stats_T] = compute_stats(quantity_to_plot_all,quantity_name,img_type,...
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
    
%     var_names = {'Delta vs Sham','Theta vs Sham','Gamma vs Sham'};
    var_names = {'Test',...
                 'Delta vs Sham p','Theta vs Sham p','Gamma vs Sham p',...
                 'Delta vs Sham p_c','Theta vs Sham p_c','Gamma vs Sham p_c',...
                 'Delta vs Sham h','Theta vs Sham h','Gamma vs Sham h'};
             
    mc_method = 'bc-h';
    
    %%
    if ~exist('roi_idxs','var')
        roi_idxs = 1:roi_no;
    end
    
    if ~exist('save_results','var')
        save_results = 1;
    end

    %%
%     p = nan(length(roi_idxs),3);
%     h = nan(length(roi_idxs),3);
    p = nan(2,3);
    p_c = nan(2,3);
    hh = nan(2,3);
    h = cell(2,3);
    h(:) = {NaN};
    tests = {};
    
    %% Calculate stats
    for roi_idx = roi_idxs
        roi_results = quantity_to_plot_all{roi_idx};
        [hn,pn] = test_normality(roi_results);
        
        % compare all stims to sham
        if all(~hn)
            [hh(roi_idx,1),p(roi_idx,1)] = ttest2(roi_results(:,1),roi_results(:,2));
            [hh(roi_idx,2),p(roi_idx,2)] = ttest2(roi_results(:,1),roi_results(:,3));
            [hh(roi_idx,3),p(roi_idx,3)] = ttest2(roi_results(:,1),roi_results(:,4));

            tests{roi_idx} = 'ttest2';

        else
            [p(roi_idx,1),hh(roi_idx,1)] = ranksum(roi_results(:,1),roi_results(:,2));
            [p(roi_idx,2),hh(roi_idx,2)] = ranksum(roi_results(:,1),roi_results(:,3));
            [p(roi_idx,3),hh(roi_idx,3)] = ranksum(roi_results(:,1),roi_results(:,4));

            tests{roi_idx} = 'ranksum';
        end
        
        [h(roi_idx,:),p_c(roi_idx,:)] = correct_significance(p(roi_idx,:),mc_method);
    end
    
        
    %% Save results
    tests = tests';
    
    if save_results
        file_name = sprintf('%s_%s_stats',img_type,quantity_name);
        save(fullfile(stats_folder,strcat(file_name,'.mat')),'p','p_c','h','tests');

%         stats_T = array2table(p,'VariableNames',var_names,'RowNames',roi_names(roi_idxs));
        stats_T = array2table([tests,num2cell(p),num2cell(p_c),h],'VariableNames',var_names,'RowNames',roi_names(roi_idxs));
        table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
        writetable(stats_T,table_name,'WriteRowNames',true);
    end
    
end
