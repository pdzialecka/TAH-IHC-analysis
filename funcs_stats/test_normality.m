function [h,p] = test_normality(data)
    %%
    % @author: pdzialecka

    % Normality test done per column
    
    % kstest: needs mean = 0, std = 1 -> (to_plot-nanmean(to_plot)./nanstd(to_plot)
    % jbtest: ok
    % lillietest: extension of kstest, no need to normalise data
    
    % in all tests:
    % h = 1 -> hypothesis than normal dist is rejected
    % h = 0 -> hypothesis cannot be rejected; approx normal
    
    %%
    cond_no = size(data,2);
    
    h = zeros(1,cond_no);
    p = nan(1,cond_no);
    
    for i = 1:cond_no
        if sum(~isnan(data(:,i))) >= 4 % at least 4 points needed for lillietest
            if ~(all(data(:,i)==0)) % disregard normalisation day 1 case
                [h(i),p(i)] = lillietest(data(:,i));
            end
%             if all(data(:,i)==0) % normalisation case; disregard
%                 h(i) = 0;
%                 p(i) = nan;
%             end
        end
    end
    
end
