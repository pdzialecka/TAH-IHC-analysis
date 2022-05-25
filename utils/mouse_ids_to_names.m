function [mouse_names] = mouse_ids_to_names(mouse_ids)
    %% Convert animal IDs (numbers) to full names
    % @author: pdzialecka
    
    % e.g. [18,31] -> {'AD-Hipp18','AD-Hipp31'}
    
    %%
    mouse_names = {};
    
    for idx = 1:length(mouse_ids)
        mouse_names{idx} = sprintf('AD-Hipp%d',mouse_ids(idx));
    end
    
end
