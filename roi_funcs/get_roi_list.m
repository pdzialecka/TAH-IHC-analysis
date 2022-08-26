function [roi_names,roi_fnames,roi_no,roi_sizes_um] = get_roi_list()
    %% Get a list of ROIs
    % @author: pdzialecka

    %%
    roi_names = {'Left DG','Right DG','Left CA1','Right CA1',...
        'Left CA3','Right CA3','Left Cortex','Right Cortex'};
    
    roi_no = length(roi_names);
    
    roi_fnames = {'hipp_DG_L','hipp_DG_R','hipp_CA1_L','hipp_CA1_R',...
                  'hipp_CA3_L','hipp_CA3_R','cortex_L','cortex_R'};
            
    roi_sizes_um = {[1300,500],[1300,500],[1300,400],[1300,400],...
                    [800,600],[800,600],[1300,500],[1300,500]};
    
end
