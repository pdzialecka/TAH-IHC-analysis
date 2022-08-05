function [roi_mask,slice_mask] = create_roi_h_masks(h_image)
    %% Create morphological roi + slice masks
    % @author: pdzialecka
    
    %%
    show_figs = 0;
    
    %% ********* ROI MASK *********
    %% Enhance image contrast
    enhance_contrast = 1;

    if enhance_contrast
        k1 = 10;
        kernel1 = 1/(k1*k1)*ones([k1,k1]);
        h_image_ = imfilter(h_image,kernel1);

%                         h_imadjust = imadjust(h_image_);
%                         h_histeq = histeq(h_image_);
        h_adapthisteq = adapthisteq(h_image_);

        h_image_ = h_adapthisteq;
        if show_figs
            figure,imshow(h_image_)
        end
    end

    %% Create binary dab image
    if enhance_contrast
        thresh = (0.7*mean(mean(h_image_))); % 150;
    else
        thresh = 200; % (0.9*mean(mean(h_image_)));
    end

    I = h_image_<thresh;
    if show_figs
        figure,imshow(I)
    end

    k2 = 10;
    kernel2 = 1/(k2*k2)*ones([k2,k2]);
    I2 = imfilter(I,kernel2);

    if show_figs
        figure,imshow(I2)
    end
    
    %% Keep only large regions in the mask
    min_con_pixels = 10e3; % 10e3
    connectivity = 8; % default: 4
    I3 = bwareaopen(I2,min_con_pixels,connectivity);
    
    if show_figs
        figure,imshow(I3)
    end
    
    %% ROI mask
    roi_mask = I3;
    
    %% ********* SLICE MASK *********
    %% Find outline of the whole slice
    slice_mask = h_image < 240; % ones(size(I));
    k3 = 100;
    kernel3 = 1/(k3*k3)*ones([k3,k3]);
    slice_mask = imfilter(slice_mask,kernel3);
    slice_mask = imfill(slice_mask,'holes');
    
    if show_figs
        figure,imshow(slice_mask)
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
    slice_mask = bwpropfilt(slice_mask,'Area',[1e7,inf]);
    slice_region = regionprops(slice_mask);
    
    if show_figs
        figure,imshow(slice_mask),hold on
        h = rectangle('Position',slice_region.BoundingBox);
        set(h,'EdgeColor','r','LineWidth',1.5);
    end

end
