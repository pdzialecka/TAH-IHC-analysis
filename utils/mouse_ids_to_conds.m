function [cond_names] = mouse_ids_to_conds(cond_idxs)
    %% Convert condition idxs (numbers) to condition names
    % @author: pdzialecka
    
    % e.g. [1,3] -> {'Sham','8 Hz'}
    
    %%
    conds = {'Sham','40 Hz','8 Hz','LTD'};
    cond_names = {};
    
    for idx = 1:length(cond_idxs)
        cond_names{idx} = conds{cond_idxs(idx)};
    end
    
end
