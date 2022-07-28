function [img_type] = find_img_type(fname)
    %% Find image type based on fname
    % @author: pdzialecka

    %%
    idxs_ = strfind(fname,'_');
    img_type = fname(idxs_(1)+1:idxs_(2)-1);
    
end
