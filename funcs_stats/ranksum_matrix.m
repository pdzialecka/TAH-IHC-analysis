function [p,h] = ranksum_matrix(X,Y,tail)
    % @author: pdzialecka
    %%
    if ~exist('tail','var')
        tail = 'both';
    end
    
    %%
    p = nan;
    h = nan;
    
    if size(Y,1) == 1 && size(Y,2) == 1
        Y_val = Y;
        Y = ones(size(X))*Y_val;
    end
    
    for i = 1:size(X,2)
        [p(i),h(i)] = ranksum(X(:,i),Y(:,i),'Tail',tail);
    end
    
end