function [files_] = rearrange_files(files,all_image_types)
    %% Rearrange files within a mouse folder
    % @author: pdzialecka
    
    % Rearranged to match the order from all_image_types
    
    %%
    files_ = files;
    img_types = {};
    
    for file_idx = 1:length(files)
        fname = files(file_idx).name;
        img_type = find_img_type(fname);
        img_idx = find(contains(all_image_types,img_type));
        
        files_(img_idx) = files(file_idx);
        img_types{img_idx} = img_type;
    end
    
    % ensure no repeats present
    if any(cellfun(@isempty,img_types))
        missing_img_types = cellfun(@isempty,img_types);
        files_(missing_img_types) = [];
    end
end
