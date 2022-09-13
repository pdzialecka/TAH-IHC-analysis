function [roi_names,roi_fnames,roi_no] = get_roi_list_IF()
    %% Get a list of ROIs (IF images)
    % @author: pdzialecka

    %%
    roi_names = {'Left DG','Right DG','Left CA1','Right CA1','Left Cortex',...
        'Right Cortex'};
    
    roi_no = length(roi_names);
    
    roi_fnames = {'DG_L','DG_R','CA1_L','CA1_R','cortex_L','cortex_R'};
                    
end
