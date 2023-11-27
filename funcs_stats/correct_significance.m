function [h_cor,p_cor] = correct_significance(p,method,idxs) % half_n)
    % @author: pdzialecka

    if ~exist('method','var')
        method = 'bc-h';
    end
    
    if ~exist('idxs','var')
        idxs = [];
    end

    %%
    if ~isempty(idxs)
        sine_idxs = idxs{1};
        square_idxs = idxs{2};
        half_n = 1;
    else
        sine_idxs = 1:length(p);
        square_idxs = [];
        half_n = 0;
    end


    n = length(p(~isnan(p))); % length(p);
    threshs = [0.05,0.01,0.001];
    mc_thresh = threshs(1);
    p_cor = nan(n,1);
    
    %%
    % bonferroni correct p value
    if strcmp(method,'bc') % bc
        if half_n
%             threshs = threshs./(n/2);
            p_cor = p.*(n/2);
        else
%             threshs = threshs./n;
            p_cor = p.*n;
        end
        
    elseif strcmp(method,'bc-h')
        p_cor_sine = bonf_holm(p(sine_idxs));
        p_cor(sine_idxs) = p_cor_sine;
        
%         if ~isempty(idxs)
        p_cor_square = bonf_holm(p(square_idxs));
        p_cor(square_idxs) = p_cor_square;
%         end
    end
    
    p_cor = p_cor';
    
    if all(isnan(p_cor)) % no correction case (simply find significance)
        p_cor = p;
    end
    
    %%
    for i = 1:length(p_cor)
%         p_i = p_cor(i); % simple: use corrected p values
        p_i = p(i); % use original p values + check corrected below mc_thresh
        p_c_i = p_cor(i);
        stars = [];
        
        if p_i < threshs(3) && p_c_i < mc_thresh % added && p_c_i < mc_thresh
            stars = '***';
        elseif p_i < threshs(2) && p_c_i < mc_thresh
            stars = '**';
        elseif p_i < threshs(1) && p_c_i < mc_thresh
            stars = '*';
        end
        
        h_cor{i} = stars;
    end
    
    if any(isnan(p))
        h_cor{isnan(p)} = nan;
    end
    
end
