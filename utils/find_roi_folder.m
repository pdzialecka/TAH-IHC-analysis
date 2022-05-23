function [roi_folder] = find_roi_folder(data_folder)
    %% Find ROI folder associated with a mouse data folder
    % @author: pdzialecka

    %%
    roi_folder = fullfile(fileparts(fileparts(data_folder)),'ROIs');

    % create subfolder with mouse name
%     mouse_name = file(1:9);
    [~,mouse_name] = fileparts(data_folder);
    roi_folder = fullfile(roi_folder,mouse_name);

    if ~exist(roi_folder)
        mkdir(roi_folder);
    end
    
end
