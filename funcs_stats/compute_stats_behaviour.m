function [p,hh] = compute_stats_behaviour(quantity_to_plot_all,quantity_name,img_type,...
                               stats_folder,save_results)

    %% Compute and save stats for behavioural results
    % @author: pdzialecka

    % not all data normally distributed; use wilcoxon test instead of ttest

    %%    
    if ~exist('save_results','var')
        save_results = 1;
    end
    
    %%
%     var_names = {'Sham','Delta','Theta','Gamma'};
    var_names = {'Sham p','Delta p','Theta p','Gamma p',...
                 'Sham p_c','Delta p_c','Theta p_c','Gamma p_c',...
                 'Sham h','Delta h','Theta h','Gamma h',...
                 'Test'};
             
    row_names = {'Cond vs sham','Cond vs chance'};
    
    mc_method = 'bc-h';

    %%
    p = nan(2,4);
    p_c = nan(2,4);
    hh = nan(2,4);
    h = cell(2,4);
    h(:) = {NaN};

    %% Calculate stats: compare each to sham
    roi_idx = 1;
    roi_results = quantity_to_plot_all{roi_idx};
    [hn,pn] = test_normality(roi_results);

    % compare all stims to sham
    if all(~hn)
%         [hh(roi_idx,1),p(roi_idx,1)] = ttest2(roi_results(:,1),roi_results(:,1));
        [hh(roi_idx,2),p(roi_idx,2)] = ttest2(roi_results(:,1),roi_results(:,2));
        [hh(roi_idx,3),p(roi_idx,3)] = ttest2(roi_results(:,1),roi_results(:,3));
        [hh(roi_idx,4),p(roi_idx,4)] = ttest2(roi_results(:,1),roi_results(:,4));

        test1 = 'ttest2';

    else
%         [p(roi_idx,1),hh(roi_idx,1)] = ranksum(roi_results(:,1),roi_results(:,1));
        [p(roi_idx,2),hh(roi_idx,2)] = ranksum(roi_results(:,1),roi_results(:,2));
        [p(roi_idx,3),hh(roi_idx,3)] = ranksum(roi_results(:,1),roi_results(:,3));
        [p(roi_idx,4),hh(roi_idx,4)] = ranksum(roi_results(:,1),roi_results(:,4));

        test1 = 'ranksum';
    end
    
    [h(roi_idx,2:end),p_c(roi_idx,2:end)] = correct_significance(p(roi_idx,2:end),mc_method);
    
    %% Calculate stats: compare each to chance
    roi_idx = 2;
    chance_level = nan;
    test2 = nan;
    
    if contains(quantity_name,'DI')
        chance_level = 0.5;
    elseif contains(quantity_name,'SAP')
        chance_level = 6/27*100;
    end
    
    if ~isnan(chance_level)
        if all(~hn)
            [hh(roi_idx,:),p(roi_idx,:)] = ttest(roi_results,chance_level);
            test2 = 'ttest';
        else
            [p(roi_idx,:),hh(roi_idx,:)] = signrank_matrix(roi_results,chance_level);
            test2 = 'signrank';
        end
        
        [h(roi_idx,:),p_c(roi_idx,:)] = correct_significance(p(roi_idx,:),mc_method);
    end
    
    %% Save results
    tests = {test1;test2};
    
    if save_results
        file_name = sprintf('%s_%s_stats',img_type,quantity_name);
        save(fullfile(stats_folder,strcat(file_name,'.mat')),'p','p_c','h','tests');

%         stats_T = array2table(p,'VariableNames',var_names,'RowNames',row_names);
        stats_T = array2table([num2cell(p),num2cell(p_c),h,tests],'VariableNames',var_names,'RowNames',row_names);
        table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
        writetable(stats_T,table_name,'WriteRowNames',true);
    end
    
end
