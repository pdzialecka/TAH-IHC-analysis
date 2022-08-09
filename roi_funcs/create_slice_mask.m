function [slice_mask,slice_mask_filled,slice_region] = create_slice_mask(h_image,file_)
    %%
    % @author: pdzialecka
    
    %%
    show_figs = 0;
    
    %%
    file = file_.name;
    [roi_folder,~] = find_roi_folder(file_.folder);

    %% Find the slice mask
    slice_mask = h_image < 240; % ones(size(I));
    k3 = 20;
    kernel3 = 1/(k3*k3)*ones([k3,k3]);
    slice_mask = imfilter(slice_mask,kernel3);
    
    min_con_pixels = 10e4;
    slice_mask = bwareaopen(slice_mask,min_con_pixels,8);
    
    % visualise inverse of the mask found
    h_image_slice_mask = labeloverlay(h_image,~slice_mask,...
        'Colormap',[0,0,1],'Transparency',0.7);
    fig = figure;
    imshow(h_image_slice_mask)
    
    fname = strcat(file(1:end-11),'_slice_mask.tif');
    saveas(fig,fullfile(roi_folder,fname));
    close(fig);

    %% Find the filled slice mask
    slice_mask_filled = imfill(slice_mask,'holes');
     
    if show_figs
        figure,imshow(slice_mask_filled)
    end
    
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
    
    %% Keep only the largest part of the mask
%     slice_mask = bwpropfilt(slice_mask,'Area',[1e7,inf]);
    slice_region = regionprops(slice_mask_filled);
    
    if show_figs
        figure,imshow(slice_mask),hold on
        h = rectangle('Position',slice_region.BoundingBox);
        set(h,'EdgeColor','r','LineWidth',1.5);
    end
    
    %% Save masks
    mask_fname = strcat(file(1:end-11),'_slice_mask.mat');
    save(fullfile(roi_folder,mask_fname),'slice_mask','slice_mask_filled',...
        'slice_region');
    
end
