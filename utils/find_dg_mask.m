function [dg_regions,dg_centroids,regions] = find_dg_mask(I_mask)
    %% Find dendate gyri within regions of a binary mask
    % @author: pdzialecka
    
    %%
    regions = regionprops(I_mask,'Area','Centroid',...
        'BoundingBox','Circularity','Eccentricity',...
        'MajorAxisLength','MinorAxisLength','ConvexArea');

    % consider only upper half of the slice
    centroids = cat(1,regions.Centroid);
    y_lim =  size(I_mask,1)/2;
    regions = regions(centroids(:,2)<y_lim);

    % middle region only
    x_lim = [size(I_mask,2)*1/3, size(I_mask,2)*2/3];
    centroids = cat(1,regions.Centroid);
    regions = regions(centroids(:,1)<x_lim(2) & centroids(:,1)>x_lim(1));

    % sort regions by area (total pixels)
    % DG = l
    [b,i] = sort([regions.Area],'descend');
    regions = regions(i);

    % exclude too big or small regions ('filled' area)
    size_exclude = [regions.ConvexArea] > 4e6 | [regions.ConvexArea] < 1e6;
    regions(size_exclude) = [];
    
    % ensure only regions with major/minor axis ratio > 2
    axis_ratio = 2;
    bb_regions = cat(1,regions.BoundingBox);
    width_to_height_ratio = bb_regions(:,3)./bb_regions(:,4);
    ratio_exclude = width_to_height_ratio < 2 | width_to_height_ratio > 5;
%     regions_axis_ratio = [regions.MajorAxisLength]./[regions.MinorAxisLength];
%     ratio_exclude = regions_axis_ratio < 2 | regions_axis_ratio > 5;
    regions(ratio_exclude) = [];
    
    % 2 biggest regions = dg
    dg_regions = regions(1:2);
    dg_centroids = cat(1,dg_regions.Centroid);
    
    % rearrange so 1 = left dg, 2 = right dg
    [~,i] = sort(dg_centroids(:,1));
    dg_regions = dg_regions(i);
    dg_centroids = round(cat(1,dg_regions.Centroid));
    
end
