function [img_type] = find_img_type(fname)
    %% Find image type based on fname
    % @author: pdzialecka

    %%
    idxs_ = strfind(fname,'_');
    
    if length(idxs_)>1
        img_type = fname(idxs_(1)+1:idxs_(2)-1);
    elseif length(idxs_) == 1
        idxs_2 = strfind(fname,'.');
        img_type = fname(idxs_(1)+1:idxs_2(1)-1);
    end
end
