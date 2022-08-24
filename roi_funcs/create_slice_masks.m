function [] = create_slice_masks(files)
    %% Create a slice mask for given files
    % @author: pdzialecka
    
    %%
    for file_idx = 1:length(files)
        
        try
            %% Directory info
            file_ = files(file_idx);
            file = file_.name;
            folder = file_.folder;

            %% Load image
            file_path = fullfile(folder,file);
    %         [h_image,dab_image,res_image] = load_deconvolved_images(file_path,1);

            image = imread(file_path);
            image = imrotate(image,-90);

            %% Convert to grayscale
            image = rgb2gray(image);

            %% Load theta
            [roi_folder,~] = find_roi_folder(folder);
            dg_fname = strcat(file(1:end-4),'_dg_points.mat');
            dg_fname_full = fullfile(roi_folder,dg_fname);

            dg_points = load(dg_fname_full);
            theta = dg_points.theta;

            %% Calculate transform matrix
            tform = create_tform(theta);

            %% Rotate image
            image = imwarp(image,tform,'interp','cubic','FillValues',255);

            %% Create slice mask
            [slice_mask,slice_mask_filled] = create_slice_mask(image,file_);
            
        catch
            fprintf('Could not create a slice mask for %s file',file)
        end
    end
end

