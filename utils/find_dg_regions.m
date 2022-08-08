function [dg_regions,dg_centroids,offset] = find_dg_regions(h_image,file_)
    %% Find DG regions in each image
    % @author: pdzialecka
    
    %%
    file = file_.name;
    [roi_folder,~] = find_roi_folder(file_.folder);
    
    pixel_size = 0.504;

    %% Find existing dg_regions file (moc23)
    base_dg_file = dir(fullfile(roi_folder,'*moc23*dg_regions*'));
    
    if ~isempty(base_dg_file) && ~contains(file,'moc23')
        base_dg = load(fullfile(base_dg_file(1).folder,base_dg_file(1).name));
        base_dg_centroids = base_dg.dg_centroids;
        
        base_compare = 1;
    else
        base_compare = 0;
    end
    
    %%
%     show_figs = 0;
%     
%     %% Enhance image contrast
%     enhance_contrast = 1;
% 
%     if enhance_contrast
%         k1 = 10;
%         kernel1 = 1/(k1*k1)*ones([k1,k1]);
%         h_image_ = imfilter(h_image,kernel1);
% 
% %                         h_imadjust = imadjust(h_image_);
% %                         h_histeq = histeq(h_image_);
%         h_adapthisteq = adapthisteq(h_image_);
% 
%         h_image_ = h_adapthisteq;
%         if show_figs
%             figure,imshow(h_image_)
%         end
%     end
% 
%     %% Create binary dab image
%     if enhance_contrast
%         thresh = (0.7*mean(mean(h_image_))); % 150;
%     else
%         thresh = 200; % (0.9*mean(mean(h_image_)));
%     end
% 
%     I = h_image_<thresh;
%     if show_figs
%         figure,imshow(I)
%     end
% 
%     k2 = 10;
%     kernel2 = 1/(k2*k2)*ones([k2,k2]);
%     I2 = imfilter(I,kernel2);
% 
%     if show_figs
%         figure,imshow(I2)
%     end
% 
%     %% Keep only large regions in the mask
%     min_con_pixels = 10e3; % 10e3
%     connectivity = 8; % default: 4
%     I3 = bwareaopen(I2,min_con_pixels,connectivity);
%     
%     if show_figs
%         figure,imshow(I3)
%     end

    %%
%     boundaries = bwboundaries(I);
%     figure,imshow(I),hold on
%     for k=1:length(boundaries)
%        b = boundaries{k};
%        plot(b(:,2),b(:,1),'g','LineWidth',3);
%     end

    %% Create morphological binary mask
%     [I_mask] = create_roi_h_mask(h_image);
% 
%     %% Find outline of the whole slice
%     slice_mask = h_image < 240; % ones(size(I));
%     k3 = 100;
%     kernel3 = 1/(k3*k3)*ones([k3,k3]);
%     slice_mask = imfilter(slice_mask,kernel3);
%     slice_mask = imfill(slice_mask,'holes');
%     
%     figure,imshow(slice_mask)
%     
%     %%
%     h_image_2 = h_image_;
%     h_image_2(~slice_mask) = nan;
%     
%     boundaries = bwboundaries(h_image_2);
%     figure,imshow(h_image_2),hold on
%     for k = 1:length(boundaries)
%        b = boundaries{k};
%        plot(b(:,2),b(:,1),'g','LineWidth',3);
%     end
% 
%     %%
%     % keep only the largest part of the mask
%     slice_mask_2 = bwpropfilt(slice_mask,'Area',[1e7,inf]);
% 
%     slice_region = regionprops(slice_mask_2);
%     
%     %%
%     figure,imshow(slice_mask_2),hold on
%     h = rectangle('Position',slice_region.BoundingBox);
%     set(h,'EdgeColor','r','LineWidth',1.5);
    
    %% Find right and left hemispheres (not in use)
%     slice_xs = round([slice_region.BoundingBox(1),slice_region.BoundingBox(1)+slice_region.BoundingBox(3)]);
%     slice_ys = round([slice_region.BoundingBox(2),slice_region.BoundingBox(1)+slice_region.BoundingBox(4)]);
%     
%     hem_r_xs = round([slice_xs(1),slice_xs(2)/2]);
%     hem_l_xs = round([hem_r_xs(2)+1,slice_xs(2)]);
%     
%     hem_r_mask = slice_mask_2(slice_ys(1):slice_ys(2),hem_r_xs(1):hem_r_xs(2));
%     figure,imshow(hem_r_mask)
%     
%     hem_l_mask = slice_mask_2(slice_ys(1):slice_ys(2),hem_l_xs(1):hem_l_xs(2));
%     figure,imshow(hem_l_mask)

    %% Find dips to calculate rotation angle (not in use)
%     % consider middle fourth of slice only
%     x_thickness_half = round(((slice_xs(2)-slice_xs(1)+1)/4)/2);
%     x_middle = round((slice_xs(2)-slice_xs(1)+1)/2);
%     middle_xs = [x_middle-x_thickness_half,x_middle+x_thickness_half];
%     
%     middle_slice_mask = slice_mask_2(slice_ys(1):slice_ys(2),middle_xs(1):middle_xs(2));
%     figure,imshow(middle_slice_mask)
%     
%     %%
%     y_thickness = round(pixel_size*1000*2); % 5 mm
% %     slice_top = slice_mask_2(slice_ys(1):slice_ys(1)+y_thickness,slice_xs(1):slice_xs(2));
%     slice_top = middle_slice_mask(1:y_thickness,:);
%     figure,imshow(slice_top)
%     
%     top_boundary = bwboundaries(slice_top);
%     tb = top_boundary{1};
% 
%     figure,imshow(slice_top),hold on
%     plot(tb(:,2),tb(:,1),'g','LineWidth',3);
%     
%     %%
% %     figure,plot(tb(:,2),tb(:,1))
%     
%     [pks,locs,w,p] = findpeaks(tb(:,1),'MinPeakProminence',100,...
%         'SortStr','descend');
% %     [b,idxs] = sort(pks,'descend');
% 
% %     dip_idx = locs(idxs(1));
%     dip_idx = locs(1);
%     mid_p1 = tb(dip_idx,:);
%     plot(mid_p1(2),mid_p1(1),'r.','MarkerSize',20);
%     
%     %%
%     thicknes_b = round(pixel_size*1000*5); % 5 mm
%     slice_bottom = slice_mask_2(slice_ys(2)-y_thickness:slice_ys(2),slice_xs(1):slice_xs(2));
%     
%     bottom_boundary = bwboundaries(slice_bottom);
%     bb = bottom_boundary{1};
% 
%     figure,imshow(slice_bottom),hold on
%     plot(bb(:,2),bb(:,1),'g','LineWidth',3);
% 
%     %%
%     
%     boundaries = bwboundaries(slice_mask_2);
%     b = boundaries{1};
%     
%     figure,imshow(slice_mask_2),hold on
% %     for k = 1:length(boundaries)
% %        b = boundaries{k};
%    plot(b(:,2),b(:,1),'g','LineWidth',3);
%     end

    %% Find initial binary mask
    [roi_mask,slice_mask] = create_roi_h_masks(h_image);
    
    %% Find initial DG regions
    [dg_regions,dg_centroids] = find_dg_mask(roi_mask);
    fig = plot_regions(roi_mask,dg_regions);

    %% Find angle based on dg
    dg_dist = dg_centroids(2,:)-dg_centroids(1,:);
    mid_point = dg_centroids(1,:)+dg_dist/2;
    
    plot(mid_point(1),mid_point(2),'g.',...
        'MarkerSize',10)
    
    x1 = dg_centroids(1,1);
    y1 = dg_centroids(1,2);
    x2 = dg_centroids(2,1);
    y2 = dg_centroids(2,2);
    
    a = y2-y1;
    b = x2-x1;

    theta = -atand(a/b); %  - 90*sign(b)
    
    %%
    trans = [0,0];
    rot = [cosd(theta),sind(theta);...
        -sind(theta),cosd(theta);];
    tform = rigid2d(rot,trans);

    h_image_rot = imwarp(h_image,tform,'interp','cubic','FillValues',255);    
    figure,imshow(h_image_rot)
        
    roi_mask_rot = imwarp(roi_mask,tform,'interp','cubic','FillValues',0);
    figure,imshow(roi_mask_rot)
    
    slice_mask_rot = imwarp(slice_mask,tform,'interp','cubic','FillValues',0);

    [x,y] = transformPointsForward(tform,mid_point(1),mid_point(2));
    mid_point_rot = round([x,y]);
    mid_x = mid_point_rot(1);
    
    xline(mid_x,'g')

    %% Offset image so midline in the middle
    new_size = size(h_image_rot);
    offset_size = new_size(2)/2-mid_x;
    
    if offset_size > 0 % right side too large
        h_image_rot = h_image_rot(:,1:end-offset_size);
        roi_mask_rot = roi_mask_rot(:,1:end-offset_size);
        
    else % left side too large
        h_image_rot = h_image_rot(:,abs(offset_size)+1:end);
        roi_mask_rot = roi_mask_rot(:,abs(offset_size)+1:end);
        mid_x = mid_x-offset_size;
        mid_point_rot(1) = mid_x;
    end
    
    figure,imshow(roi_mask_rot)
    xline(mid_x,'g')
    
    %%
%                     BWoutline = bwperim(slice_mask);
%                     Segout = zeros(size(I)); 
%                     Segout(BWoutline) = 255;
%                     figure,imshow(Segout)

    %%
%                     bw = activecontour(I,mask,20);
%                     figure,imshow(bw)
    
    %% Find DG regions
    [dg_regions,dg_centroids] = find_dg_mask(roi_mask_rot);
    
%     regions = regionprops(I3,'Area','Centroid',...
%         'BoundingBox','Circularity','Eccentricity',...
%         'MajorAxisLength','MinorAxisLength','ConvexArea');
% 
%     % consider only upper half of the slice
%     centroids = cat(1,regions.Centroid);
%     y_lim =  size(I3,1)/2;
%     regions = regions(centroids(:,2)<y_lim);
% 
%     % middle region only
%     x_lim = [size(I3,2)*1/3, size(I3,2)*2/3];
%     centroids = cat(1,regions.Centroid);
%     regions = regions(centroids(:,1)<x_lim(2) & centroids(:,1)>x_lim(1));
% 
%     % sort regions by area (total pixels)
%     % DG = l
%     [b,i] = sort([regions.Area],'descend');
%     regions = regions(i);
% 
%     % exclude too big or small regions ('filled' area)
%     size_exclude = [regions.ConvexArea] > 4e6 | [regions.ConvexArea] < 1e6;
%     regions(size_exclude) = [];
%     
%     % ensure only regions with major/minor axis ratio > 2
%     axis_ratio = 2;
%     bb_regions = cat(1,regions.BoundingBox);
%     width_to_height_ratio = bb_regions(:,3)./bb_regions(:,4);
%     ratio_exclude = width_to_height_ratio < 2 | width_to_height_ratio > 5;
% %     regions_axis_ratio = [regions.MajorAxisLength]./[regions.MinorAxisLength];
% %     ratio_exclude = regions_axis_ratio < 2 | regions_axis_ratio > 5;
%     regions(ratio_exclude) = [];
%     
%     % 2 biggest regions = dg
%     dg_regions = regions(1:2);
%     dg_centroids = cat(1,dg_regions.Centroid);
%     
%     % rearrange so 1 = left dg, 2 = right dg
%     [~,i] = sort(dg_centroids(:,1));
%     dg_regions = dg_regions(i);
%     dg_centroids = round(cat(1,dg_regions.Centroid));

    
    
%     plot_regions = dg_regions;
%     plot_centroids = dg_centroids;
%     plot_idxs = 1:2;
    
%     plot_regions = regions;
%     plot_centroids = centroids;
%     plot_idxs = 1:length(regions);

    fig1 = plot_regions(roi_mask_rot,dg_regions);

%     fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
%     imshow(I3),hold on
% 
%     for idx = plot_idxs
%         h = rectangle('Position',plot_regions(idx).BoundingBox);
%         set(h,'EdgeColor','r','LineWidth',1.5);
%         plot(plot_centroids(idx,1),plot_centroids(idx,2),'r.',...
%             'MarkerSize',10)
%         1;
%     end
    
    if base_compare
        for idx = 1:2
            h = rectangle('Position',base_dg.dg_regions(idx).BoundingBox,...
                'LineStyle',':');
            set(h,'EdgeColor','g','LineWidth',1.25);
            plot(base_dg_centroids(idx,1),base_dg_centroids(idx,2),'g.',...
                'MarkerSize',10)
            1;
        end
    end
    
    fname = strcat(file(1:end-11),'_dg_regions.tif');
    saveas(fig1,fullfile(roi_folder,fname));
    close(fig1);
    
    %%
%     figure,imshow(slice_mask_rot)
%     slice_region = regionprops(slice_mask_rot,'Area','Centroid',...
%         'BoundingBox','Circularity','Eccentricity',...
%         'MajorAxisLength','MinorAxisLength','ConvexArea');
%     h = rectangle('Position',slice_region(1).BoundingBox,...
%                 'LineStyle',':');
%     set(h,'EdgeColor','g','LineWidth',1.25);
    
    %% Find all ROIs
%     dg_x_start = dg_regions(1).BoundingBox(1);
%     dg_y_start = dg_regions(1).BoundingBox(2);
%     dg_w = dg_regions(1).BoundingBox(3);
%     dg_h = dg_regions(1).BoundingBox(4);
    
    dg_dims = round([1300,600]/pixel_size); % in um
    dg_centroid_L = dg_centroids(1,:);
    dg_centroid_R = dg_centroids(2,:);

    % coordinates
    % start = upper left corner on the pic
    dg_start_L = [dg_centroid_L(1)-dg_dims(1)/2 dg_centroid_L(2)-dg_dims(2)/2];
    coords_dg_L = [dg_start_L(1) dg_start_L(2) dg_dims(1) dg_dims(2)];
    
    dg_start_R = [dg_centroid_R(1)-dg_dims(1)/2 dg_centroid_R(2)-dg_dims(2)/2];
    coords_dg_R = [dg_start_R(1) dg_start_R(2) dg_dims(1) dg_dims(2)];
    
    % ca1 x higher than dg
    ca1_dims = round([1300,500]/pixel_size);
    h_above_dg = round(100/pixel_size);
    coords_ca1_L = [dg_start_L(1) dg_start_L(2)-ca1_dims(2)-h_above_dg ca1_dims(1) ca1_dims(2)];
    coords_ca1_R = [dg_start_R(1) dg_start_R(2)-ca1_dims(2)-h_above_dg ca1_dims(1) ca1_dims(2)];
    
    % cortex x higher than ca1
    cortex_dims = round([1300,700]/pixel_size);
    h_above_ca1 = round(100/pixel_size);
    coords_cortex_L = [coords_ca1_L(1) coords_ca1_L(2)-cortex_dims(2)-h_above_ca1 cortex_dims(1) cortex_dims(2)];
    coords_cortex_R = [coords_ca1_R(1) coords_ca1_R(2)-cortex_dims(2)-h_above_ca1 cortex_dims(1) cortex_dims(2)];

    % ca3? tricky as varies between images
    ca3_dims = round([750,800]/pixel_size);
    w_from_dg = round(100/pixel_size);
    
    % OPTION 1: includes ca2
    coords_ca3_L = [dg_start_L(1)-ca3_dims(1)-w_from_dg dg_start_L(2) ca3_dims(1) ca3_dims(2)];
    coords_ca3_R = [dg_start_R(1)+dg_dims(1)+w_from_dg dg_start_R(2) ca3_dims(1) ca3_dims(2)];
    
    % OPTION 2: height only up to dg centroid
    ca3_dims = round([750,500]/pixel_size);
    coords_ca3_L = [dg_start_L(1)-ca3_dims(1)-w_from_dg dg_centroid_L(2) ca3_dims(1) ca3_dims(2)];
    coords_ca3_R = [dg_start_R(1)+dg_dims(1)+w_from_dg dg_centroid_R(2) ca3_dims(1) ca3_dims(2)];


    fig2 = plot_regions(h_image_rot,dg_regions);colormap(h_colormap)
    dg_roi_L = drawrectangle('Position',coords_dg_L,'Color','g');
    ca1_roi_L = drawrectangle('Position',coords_ca1_L,'Color','m');
    cortex_roi_L = drawrectangle('Position',coords_cortex_L,'Color','y');
    ca3_roi_L = drawrectangle('Position',coords_ca3_L,'Color','b');

    dg_roi_R = drawrectangle('Position',coords_dg_R,'Color','g');
    ca1_roi_R = drawrectangle('Position',coords_ca1_R,'Color','m');
    cortex_roi_R = drawrectangle('Position',coords_cortex_R,'Color','y');
    ca3_roi_R = drawrectangle('Position',coords_ca3_R,'Color','b');

    %% Calculate offset
    if base_compare
        offset = dg_centroids - base_dg_centroids;
    else
        offset = [0 0; 0 0];
    end
    
    %% Save DG regions
    dg_fname = strcat(file(1:end-11),'_dg_regions.mat');
    save(fullfile(roi_folder,dg_fname),'dg_regions','dg_centroids','offset');
    
end
    