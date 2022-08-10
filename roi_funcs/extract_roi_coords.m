function [coords,this_roi,coords_field] = extract_roi_coords(rois,roi_fname)
    %% 
    % @author: pdzialecka
    
    %% Extract ROI info
    dg_rois = rois.dg_rois;
    ca1_rois = rois.ca1_rois;
    ca3_rois = rois.ca3_rois;
    cortex_rois = rois.cortex_rois;

    %% Find the correct ROI + coords
    if contains(roi_fname,'R')
        coords_field = 'R_coords';
    elseif contains(roi_fname,'L')
        coords_field = 'L_coords';
    end

    if contains(roi_fname,'DG')
        this_roi = dg_rois;
    elseif contains(roi_fname,'CA1')
        this_roi = ca1_rois;
    elseif contains(roi_fname,'CA3')
        this_roi = ca3_rois;
    elseif contains(roi_fname,'cortex')
        this_roi = cortex_rois;
    end

    coords = getfield(this_roi,coords_field);
end
