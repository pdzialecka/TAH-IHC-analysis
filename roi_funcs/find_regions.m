function [dg_regions,dg_centroids,offset,theta,rois] = find_regions(h_image,file_,...
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
%     pixel_size = 0.504;
    [h_colormap,dab_colormap] = create_hdab_colormaps();
    theta = 0;
    rois = {};
    
    if find_all_rois
        rotate_slice = 1;
    else
        rotate_slice = 0;
    end
    
    show_figs = 0;

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

        %% Find angle based on dg
        mid_point = find_dg_mid_point(dg_centroids);
        
        if show_figs
            plot_regions(roi_mask,dg_regions);
            plot(mid_point(1),mid_point(2),'g.',...
                'MarkerSize',10)
        end

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

        %%
        h_image_rot = imwarp(h_image,tform,'interp','cubic','FillValues',255);    
        if show_figs
            figure,imshow(h_image_rot)
        end
        
        roi_mask_rot = imwarp(roi_mask,tform,'interp','cubic','FillValues',0);
        slice_mask_rot = imwarp(slice_mask,tform,'interp','cubic','FillValues',0);

        % didn't work in some cases - not sure why
%         [x,y] = transformPointsForward(tform,mid_point(1),mid_point(2));
%         mid_point_rot = round([x,y]);

        [dg_regions_rot,dg_centroids_rot] = find_mask_dg(roi_mask_rot);
        mid_point_rot = find_dg_mid_point(dg_centroids_rot);
        mid_x = mid_point_rot(1);

        if show_figs
            figure,imshow(roi_mask_rot),hold on
            plot(mid_point_rot(1),mid_point_rot(2),'g.',...
                'MarkerSize',10)
            xline(mid_x,'g')
        end

        %% Offset image so midline in the middle
        new_size = size(h_image_rot);
        mid_offset = new_size(2)/2-mid_x;

        if mid_offset > 0 % right side too large
            h_image_rot = h_image_rot(:,1:end-mid_offset);
            roi_mask_rot = roi_mask_rot(:,1:end-mid_offset);

        else % left side too large
            h_image_rot = h_image_rot(:,abs(mid_offset)+1:end);
            roi_mask_rot = roi_mask_rot(:,abs(mid_offset)+1:end);
            mid_x = mid_x-mid_offset;
            mid_point_rot(1) = mid_x;
        end

        %% Final version
        if show_figs
            figure,imshow(h_image_rot); hold on
            colormap(h_colormap)
            plot(mid_point_rot(1),mid_point_rot(2),'g.',...
                'MarkerSize',10)
            xline(mid_x,'g')
        end
        
        %%
        h_image = h_image_rot;
        roi_mask = roi_mask_rot;
    end
    
    %% Find DG regions
    [dg_regions,dg_centroids] = find_mask_dg(roi_mask);
    mid_point = find_dg_mid_point(dg_centroids);

    % overlay dg regions on mask
    fig1 = plot_regions(roi_mask,dg_regions);
    plot(mid_point_rot(1),mid_point_rot(2),'r.',...
            'MarkerSize',10)
    xline(mid_point(1),'r','LineWidth',4)

    
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
    
    fname = strcat(file(1:end-11),'_dg_regions_mask.tif');
    saveas(fig1,fullfile(roi_folder,fname));
    close(fig1);
    
    % overlay dg regions on h image
    fig2 = plot_regions(h_image,dg_regions);
    colormap(h_colormap)
    plot(mid_point_rot(1),mid_point_rot(2),'r.',...
            'MarkerSize',10)
    xline(mid_point(1),'r','LineWidth',3)
    
    fname = strcat(file(1:end-11),'_dg_regions.tif');
    saveas(fig2,fullfile(roi_folder,fname));
    close(fig2);
    
    %% Calculate offset
    if base_compare
        offset = dg_centroids - base_dg_centroids;
    else
        offset = [0 0; 0 0];
    end
    
    %% Save DG regions
    dg_fname = strcat(file(1:end-11),'_dg_regions.mat');
    save(fullfile(roi_folder,dg_fname),'dg_regions','dg_centroids','theta',...
        'offset','mid_offset');
    
    %% Find all ROIs
    if find_all_rois
        dg_dims_um = [1300,600];
        dg_dims = um_to_pixel(dg_dims_um);
        dg_centroid_L = dg_centroids(1,:);
        dg_centroid_R = dg_centroids(2,:);

        % coordinates. start = upper left corner on the pic
        dg_start_L = round([dg_centroid_L(1)-dg_dims(1)/2 dg_centroid_L(2)-dg_dims(2)/2]);
        dg_L_coords = [dg_start_L(1) dg_start_L(2) dg_dims(1) dg_dims(2)];

        dg_start_R = round([dg_centroid_R(1)-dg_dims(1)/2 dg_centroid_R(2)-dg_dims(2)/2]);
        dg_R_coords = [dg_start_R(1) dg_start_R(2) dg_dims(1) dg_dims(2)];

        % ca1 100um higher than dg
        ca1_dims_um = [1300,500];
        ca1_dims = um_to_pixel(ca1_dims_um);
        ca1_dg_h_um = 100;
        ca1_dg_h = um_to_pixel(ca1_dg_h_um);
        ca1_L_coords = [dg_start_L(1) dg_start_L(2)-ca1_dims(2)-ca1_dg_h ca1_dims(1) ca1_dims(2)];
        ca1_R_coords = [dg_start_R(1) dg_start_R(2)-ca1_dims(2)-ca1_dg_h ca1_dims(1) ca1_dims(2)];

        % cortex 100um higher than ca1
        cortex_dims_um = [1300,700];
        cortex_dims = um_to_pixel(cortex_dims_um);
        cortex_ca1_h_um = 100;
        cortex_ca1_h = um_to_pixel(cortex_ca1_h_um);
        cortex_L_coords = [ca1_L_coords(1) ca1_L_coords(2)-cortex_dims(2)-cortex_ca1_h cortex_dims(1) cortex_dims(2)];
        cortex_R_coords = [ca1_R_coords(1) ca1_R_coords(2)-cortex_dims(2)-cortex_ca1_h cortex_dims(1) cortex_dims(2)];

        % ca3 is tricky, varies between images
        ca3_dg_w_um = 100;
        ca3_dg_w = um_to_pixel(ca3_dg_w_um);

        % OPTION 1: includes ca2
        ca3_dims_um = [750,800];
        ca3_dims = um_to_pixel(ca3_dims_um);
        ca3_L_coords = [dg_start_L(1)-ca3_dims(1)-ca3_dg_w dg_start_L(2) ca3_dims(1) ca3_dims(2)];
        ca3_R_coords = [dg_start_R(1)+dg_dims(1)+ca3_dg_w dg_start_R(2) ca3_dims(1) ca3_dims(2)];

        % OPTION 2: height only up to dg centroid. may be not good enough
%         ca3_dims = round([750,500]/pixel_size);
%         ca3_L_coords = [dg_start_L(1)-ca3_dims(1)-ca3_dg_w dg_centroid_L(2) ca3_dims(1) ca3_dims(2)];
%         ca3_R_coords = [dg_start_R(1)+dg_dims(1)+ca3_dg_w dg_centroid_R(2) ca3_dims(1) ca3_dims(2)];


        %% Plot auto ROIs
        fig3 = plot_regions(h_image,dg_regions);colormap(h_colormap)
        
        lw = 2.5;
        col = 'k';
        ls = '--';
        
        dg_L_roi = rectangle('Position',dg_L_coords,'EdgeColor',col,'LineWidth',lw,'LineStyle',ls);
        ca1_L_roi = rectangle('Position',ca1_L_coords,'EdgeColor',col,'LineWidth',lw,'LineStyle',ls);
        cortex_L_roi = rectangle('Position',cortex_L_coords,'EdgeColor',col,'LineWidth',lw,'LineStyle',ls);
        ca3_L_roi = rectangle('Position',ca3_L_coords,'EdgeColor',col,'LineWidth',lw,'LineStyle',ls);

        dg_R_roi = rectangle('Position',dg_R_coords,'EdgeColor',col,'LineWidth',lw,'LineStyle',ls);
        ca1_R_roi = rectangle('Position',ca1_R_coords,'EdgeColor',col,'LineWidth',lw,'LineStyle',ls);
        cortex_R_roi = rectangle('Position',cortex_R_coords,'EdgeColor',col,'LineWidth',lw,'LineStyle',ls);
        ca3_R_roi = rectangle('Position',ca3_R_coords,'EdgeColor',col,'LineWidth',lw,'LineStyle',ls);
        
        fname = strcat(file(1:end-11),'_rois_auto.tif');
        saveas(fig3,fullfile(roi_folder,fname));
        close(fig3);
    
        %% Store ROIs
        dg_rois.name = 'DG';
        dg_rois.dims_um = dg_dims_um;
        dg_rois.dims = dg_dims;
        dg_rois.L_coords = dg_L_coords;
        dg_rois.R_coords = dg_R_coords;

        ca1_rois.name = 'CA1';
        ca1_rois.dims_um = ca1_dims_um;
        ca1_rois.dims = ca1_dims;
        ca1_rois.L_coords = ca1_L_coords;
        ca1_rois.R_coords = ca1_R_coords;
        ca1_rois.ca1_dg_h_um = ca1_dg_h_um;
        ca1_rois.ca1_dg_h = ca1_dg_h;

        ca3_rois.name = 'CA3';
        ca3_rois.dims_um = ca3_dims_um;
        ca3_rois.dims = ca3_dims;
        ca3_rois.L_coords = ca3_L_coords;
        ca3_rois.R_coords = ca3_R_coords;
        ca3_rois.ca3_dg_w_um = ca3_dg_w_um;
        ca3_rois.ca3_dg_w = ca3_dg_w;

        cortex_rois.name = 'Cortex';
        cortex_rois.dims_um = cortex_dims_um;
        cortex_rois.dims = cortex_dims;
        cortex_rois.L_coords = cortex_L_coords;
        cortex_rois.R_coords = cortex_R_coords;
        cortex_rois.cortex_ca1_h_um = cortex_ca1_h_um;
        cortex_rois.cortex_ca1_h = cortex_ca1_h;


        rois = {dg_rois,ca1_rois,ca3_rois,cortex_rois};
        
        %% Save ROIs found
        dg_roi_fname = strcat(file(1:end-11),'_rois_auto.mat');
        save(fullfile(roi_folder,dg_roi_fname),'dg_rois','ca1_rois',...
            'ca3_rois','cortex_rois');

    end
    
end
    