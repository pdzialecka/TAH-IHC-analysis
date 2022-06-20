function [roi_names,roi_fnames,roi_no] = get_roi_list()
    %% Get a list of ROIs
    % @author: pdzialecka

    %%
    roi_names = {'Right Hipp','Left Hipp','Right Cortex','Left Cortex','Right Auditory Cortex'};
    roi_no = length(roi_names);
    
    roi_fnames = {'hipp_R','hipp_L','cortex_R','cortex_L','cortex_AU'};
    
end
