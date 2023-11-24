function [p,p_c,h_c,test] = test_significance_2s(data_1,data_2,tail,mc)
    % @author: pdzialecka
    
    %%
    if ~exist('tail','var')
        tail = 'both';
    end
    
    if ~exist('mc','var')
        mc = 1;
    end

    %% Test data normality
    [hn_1,pn_1] = test_normality(data_1);
    [hn_2,pn_2] = test_normality(data_2);

    hn = [hn_1, hn_2];
    
    %%
    if mc
        mc_method = 'bc-h';
    else
        mc_method = '';
    end
    
    %%
%     if isempty(idxs)
%         data_1 = data;
%         data_2 = 0;
%     else
%         data_1 = data(idxs{1},:); % data(:,idxs{1});
%         data_2 = data(idxs{2},:); % data(:,idxs{2});
%     end

    %% Calculate p value
    if all(~hn)
        [~,p] = ttest2(data_1,data_2,'Tail',tail);
        test = 'ttest2';
    else
        [p,~] = ranksum_matrix(data_1,data_2,tail); % data signrank_matrix(data_1,data_2,tail);
        test = 'ranksum';
    end
    
    %% MC correction
    [h_c,p_c] = correct_significance(p,mc_method);

end
