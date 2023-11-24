function [dg_regions,dg_centroids,regions] = find_mask_dg(roi_mask)
    %% Find dendate gyri within regions of a binary mask
    % @author: pdzialecka
    
    %%
    % consider only upper half of the slice + middle 3/5
%     height_size = round([1,size(roi_mask,1)/2]);
% %     width_size = round([size(roi_mask,2)*1/5,size(roi_mask,2)*4/5]);
% %     roi_mask_ = roi_mask(height_size(1):height_size(2),width_size(1):width_size(2));
%     roi_mask_ = roi_mask(height_size(1):height_size(2),1:end);

    % exclude too big or small regions ('filled' area)
    roi_mask_ = bwpropfilt(roi_mask,'ConvexArea',[0.5e6,4e6]);
    
    % exclude too long or short regions
    roi_mask_ = bwpropfilt(roi_mask_,'MajorAxisLength',[1e3,3e3]);
    
    regions = regionprops(roi_mask_,'Area','Centroid',...
        'BoundingBox','Circularity','Eccentricity',...
        'MajorAxisLength','MinorAxisLength','ConvexArea',...
        'Perimeter','EquivDiameter');

    % consider only upper half of the slice
    centroids = cat(1,regions.Centroid);
    y_lim =  size(roi_mask,1)/2;
    regions = regions(centroids(:,2)<y_lim);

    % middle region only
    x_lim = round([size(roi_mask,2)*1/5, size(roi_mask,2)*4/5]);
    centroids = cat(1,regions.Centroid);
    regions = regions(centroids(:,1)<x_lim(2) & centroids(:,1)>x_lim(1));

    % sort regions by area (total pixels)
    % DG = l
    [b,i] = sort([regions.Area],'descend');
    regions = regions(i);
    
    % find regions within the right height only
    h = [];
    for idx = 1:length(regions)
        h(idx) = regions(idx).BoundingBox(4);
    end
    h_exclude = h < 500 | h > 1100;
    regions(h_exclude) = [];
   
    
%     plot_regions(roi_mask,regions);

    
    % ensure only regions with major/minor axis ratio > 2
    bb_regions = cat(1,regions.BoundingBox);
    width_to_height_ratio = bb_regions(:,3)./bb_regions(:,4);
    ratio_exclude = width_to_height_ratio < 1 | width_to_height_ratio > 5;
%     regions_axis_ratio = [regions.MajorAxisLength]./[regions.MinorAxisLength];
%     ratio_exclude = regions_axis_ratio < 2 | regions_axis_ratio > 5;
    regions(ratio_exclude) = [];
    
    % accept only regions with large perimeter relative to major axis
    per_ratio_exclude = [regions.MajorAxisLength]./[regions.Perimeter] > 0.3;
    regions(per_ratio_exclude) = [];
    
    % correlation with template? complex
%     for idx = 1:length(regions)
%         region_img = regions(idx).Image;
%         cc = xcorr2(single(region_img),single(dg_template));
%         C = normxcorr2(dg_template,region_img);
%         figure,imagesc(cc)
%         figure,imagesc(C)
%     end

    
    % 2 biggest regions = dg
    dg_regions = regions(1:2);
    dg_centroids = cat(1,dg_regions.Centroid);
    
    % rearrange so 1 = left dg, 2 = right dg
    [~,i] = sort(dg_centroids(:,1));
    dg_regions = dg_regions(i);
    dg_centroids = round(cat(1,dg_regions.Centroid));
    
end
