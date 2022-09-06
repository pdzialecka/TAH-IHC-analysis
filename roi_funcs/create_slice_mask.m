function [slice_mask,slice_mask_filled,slice_region,mask_thresh] = create_slice_mask(image,file_,mask_thresh)
    %% Create slice mask based on image
    % @author: pdzialecka
    
    % The function should be used on rotated h_image (if auto rois selected)
    
    %%
    if ~exist('mask_thresh','var')
        mask_thresh = 240;
    end
    
    %%
    show_figs = 0;
    
    %%
    file = file_.name;
    [roi_folder,~] = find_roi_folder(file_.folder);
    
    %% Two exceptions - manual thresholds
    if strcmp(file,'AD-Hipp48_iba1.tif') || strcmp(file,'AD-Hipp51_12f4.tif')
        % 244 (mask_thresh+10) a bit better out of slice for 51 but taking
        % too much out of dg
        mask_thresh = 245;
    end

    %% Find the slice mask
    slice_mask = image < mask_thresh;
    k3 = 20;
    kernel3 = 1/(k3*k3)*ones([k3,k3]);
    slice_mask = imfilter(slice_mask,kernel3,'replicate');
    
    min_con_pixels = 10e4;
    slice_mask = bwareaopen(slice_mask,min_con_pixels,8);

    %% Find boundaries of the slice mask
%     h_image_2 = h_image_;
%     h_image_2(~slice_mask) = nan;
%     
%     boundaries = bwboundaries(h_image_2);
%     
%     if show_figs
%         figure,imshow(h_image_2),hold on
%         for k = 1:length(boundaries)
%            b = boundaries{k};
%            plot(b(:,2),b(:,1),'g','LineWidth',3);
%         end
%     end

    %% Visualise inverse of the slice mask found
    h_image_slice_mask = labeloverlay(image,~slice_mask,...
        'Colormap',[0,0,1],'Transparency',0.7);
    fig = figure;
    imshow(h_image_slice_mask)
    
%     fname = strcat(file(1:end-11),'_slice_mask.tif'); % for h_image
    fname = strcat(file(1:end-4),'_slice_mask.tif');
    saveas(fig,fullfile(roi_folder,fname));
    close(fig);
    
    %% Find the filled slice mask
    slice_mask_filled = imfill(slice_mask,'holes');
     
    if show_figs
        figure,imshow(slice_mask_filled)
    end

    %% Find the slice region
%     slice_mask = bwpropfilt(slice_mask,'Area',[1e7,inf]);
    slice_region = regionprops(slice_mask_filled);
    slice_region = slice_region([slice_region.Area] == max([slice_region.Area]));
    
    if show_figs
        figure,imshow(slice_mask),hold on
        h = rectangle('Position',slice_region.BoundingBox);
        set(h,'EdgeColor','r','LineWidth',1.5);
    end
    
    %% Remove manually selected parts of the mask
    img_type = find_img_type(file);
    remove_files = dir(fullfile(roi_folder,strcat('*',img_type,'*remove_mask.mat')));
    
    if ~isempty(remove_files)
        
        new_slice_mask = slice_mask;
        total_remove_mask = false(size(slice_mask));
        for i = 1:length(remove_files)
            remove_mask = load(fullfile(remove_files(i).folder,remove_files(i).name)).remove_mask;
            total_remove_mask(remove_mask) = 1;
        end
        new_slice_mask(total_remove_mask) = 0;
        
        % Replot the slice mask figure
        h_image_slice_mask = labeloverlay(image,~slice_mask,...
            'Colormap',[0,0,1],'Transparency',0.7);
        h_image_slice_mask_2 = labeloverlay(h_image_slice_mask,total_remove_mask,...
            'Colormap',[1,0,0],'Transparency',0.7);
        
        fig = figure;
        imshow(h_image_slice_mask_2)
%         fname = strcat(file(1:end-11),'_slice_mask.tif');
        fname = strcat(file(1:end-4),'_slice_mask.tif');
        saveas(fig,fullfile(roi_folder,fname));
        close(fig);
        
        
        slice_mask = new_slice_mask;
%         figure,imshow(slice_mask)
    end
    
        
    %% Save masks
%     mask_fname = strcat(file(1:end-11),'_slice_mask.mat');
    mask_fname = strcat(file(1:end-4),'_slice_mask.mat');
    save(fullfile(roi_folder,mask_fname),'slice_mask','slice_mask_filled',...
        'slice_region','mask_thresh');
    
end
