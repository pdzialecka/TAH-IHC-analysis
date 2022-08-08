function [dg_regions,dg_centroids,offset] = find_regions(h_image,file_,...
    find_all_rois)

    %% Find regions in each image: DG and (optional) all other ROIs
    % @author: pdzialecka
    
    %%
%     if ~exist('rotate_slice','var')
%         rotate_slice = 1;
%     end
    
    if ~exist('find_all_rois','var')
        find_all_rois = 1;
    end
    
    %%
    pixel_size = 0.504;
    theta = 0;
    
    if find_all_rois
        rotate_slice = 1;
    else
        rotate_slice = 0;
    end

    %%
    file = file_.name;
    [roi_folder,~] = find_roi_folder(file_.folder);
    
    %% Find existing dg_regions file (moc23)
    base_dg_file = dir(fullfile(roi_folder,'*moc23*dg_regions*'));
    
    if ~isempty(base_dg_file) && ~contains(file,'moc23')
        base_dg = load(fullfile(base_dg_file(1).folder,base_dg_file(1).name));
        base_dg_centroids = base_dg.dg_centroids;
        
        base_compare = 1;
    else
        base_compare = 0;
    end

    
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

    %% Find roi binary mask
    [roi_mask,slice_mask] = create_roi_h_masks(h_image);
    
    %% Find initial DG regions    
    if rotate_slice
        [dg_regions,dg_centroids] = find_mask_dg(roi_mask);
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

        theta = -atand(a/b);

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
        roi_mask = roi_mask_rot;
    end
    
    %% Find DG regions
    [dg_regions,dg_centroids] = find_mask_dg(roi_mask);
    
    fig1 = plot_regions(roi_mask,dg_regions);

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
    
    %% Calculate offset
    if base_compare
        offset = dg_centroids - base_dg_centroids;
    else
        offset = [0 0; 0 0];
    end
    
    %% Save DG regions
    dg_fname = strcat(file(1:end-11),'_dg_regions.mat');
    save(fullfile(roi_folder,dg_fname),'dg_regions','dg_centroids','theta','offset');
    
    %% Find all ROIs
    if find_all_rois
        dg_dims = round([1300,600]/pixel_size); % in um
        dg_centroid_L = dg_centroids(1,:);
        dg_centroid_R = dg_centroids(2,:);

        % coordinates. start = upper left corner on the pic
        dg_start_L = [dg_centroid_L(1)-dg_dims(1)/2 dg_centroid_L(2)-dg_dims(2)/2];
        dg_L_coords = [dg_start_L(1) dg_start_L(2) dg_dims(1) dg_dims(2)];

        dg_start_R = [dg_centroid_R(1)-dg_dims(1)/2 dg_centroid_R(2)-dg_dims(2)/2];
        dg_R_coords = [dg_start_R(1) dg_start_R(2) dg_dims(1) dg_dims(2)];

        % ca1 100um higher than dg
        ca1_dims = round([1300,500]/pixel_size);
        ca1_dg_h = round(100/pixel_size);
        ca1_L_coords = [dg_start_L(1) dg_start_L(2)-ca1_dims(2)-ca1_dg_h ca1_dims(1) ca1_dims(2)];
        ca1_R_coords = [dg_start_R(1) dg_start_R(2)-ca1_dims(2)-ca1_dg_h ca1_dims(1) ca1_dims(2)];

        % cortex 100um higher than ca1
        cortex_dims = round([1300,700]/pixel_size);
        cortex_ca1_h = round(100/pixel_size);
        cortex_L_coords = [ca1_L_coords(1) ca1_L_coords(2)-cortex_dims(2)-cortex_ca1_h cortex_dims(1) cortex_dims(2)];
        cortex_R_coords = [ca1_R_coords(1) ca1_R_coords(2)-cortex_dims(2)-cortex_ca1_h cortex_dims(1) cortex_dims(2)];

        % ca3 is tricky, varies between images
        ca3_dims = round([750,800]/pixel_size);
        ca3_dg_w = round(100/pixel_size);

        % OPTION 1: includes ca2
%         coords_ca3_L = [dg_start_L(1)-ca3_dims(1)-w_from_dg dg_start_L(2) ca3_dims(1) ca3_dims(2)];
%         coords_ca3_R = [dg_start_R(1)+dg_dims(1)+w_from_dg dg_start_R(2) ca3_dims(1) ca3_dims(2)];

        % OPTION 2: height only up to dg centroid
        ca3_dims = round([750,500]/pixel_size);
        coords_ca3_L = [dg_start_L(1)-ca3_dims(1)-ca3_dg_w dg_centroid_L(2) ca3_dims(1) ca3_dims(2)];
        coords_ca3_R = [dg_start_R(1)+dg_dims(1)+ca3_dg_w dg_centroid_R(2) ca3_dims(1) ca3_dims(2)];


        fig2 = plot_regions(h_image_rot,dg_regions);colormap(h_colormap)
        dg_L_roi = drawrectangle('Position',dg_L_coords,'Color','g');
        ca1_L_roi = drawrectangle('Position',ca1_L_coords,'Color','m');
        cortex_L_roi = drawrectangle('Position',cortex_L_coords,'Color','y');
        ca3_L_roi = drawrectangle('Position',coords_ca3_L,'Color','b');

        dg_R_roi = drawrectangle('Position',dg_R_coords,'Color','g');
        ca1_R_roi = drawrectangle('Position',ca1_R_coords,'Color','m');
        cortex_R_roi = drawrectangle('Position',cortex_R_coords,'Color','y');
        ca3_R_roi = drawrectangle('Position',coords_ca3_R,'Color','b');
    end
    
end
    