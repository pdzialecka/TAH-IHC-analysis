function [dg_regions,dg_centroids,offset] = find_dg_regions(h_image,file_)
    %% Find DG regions in each image
    % @author: pdzialecka
    
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
    
    %%
    show_figs = 0;
    
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
        thresh = 150;
    else
        thresh = 200; % (0.9*mean(mean(h_image_))); % 200
    end

    I = h_image_<thresh;

    k2 = 50;
    kernel2 = 1/(k2*k2)*ones([k2,k2]);
    I2 = imfilter(I,kernel2);

    if show_figs
        figure,imshow(I2)
    end

    %%
%                     boundaries = bwboundaries(I);
%                     figure,imshow(I),hold on
%                     for k=1:length(boundaries)
%                        b = boundaries{k};
%                        plot(b(:,2),b(:,1),'g','LineWidth',3);
%                     end

    %% Find outline of the whole slice
%                     slice_mask = h_image < 240; % ones(size(I));
%                     slice_mask = imfilter(mask,kernel2);
%                     slice_mask = imfill(slice_mask,'holes');

    %%
%                     BWoutline = bwperim(slice_mask);
%                     Segout = zeros(size(I)); 
%                     Segout(BWoutline) = 255;
%                     figure,imshow(Segout)

    %%
%                     bw = activecontour(I,mask,20);
%                     figure,imshow(bw)

    %% Keep only large regions in the mask
    min_con_pixels = 10e3; % max ~10 um long
    connectivity = 8; % default: 4
    I3 = bwareaopen(I2,min_con_pixels,connectivity);
    
    if show_figs
        figure,imshow(I3)
    end
    
    %% Find DG regions
    regions = regionprops(I3,'Area','Centroid',...
        'BoundingBox','Circularity','Eccentricity',...
        'MajorAxisLength','MinorAxisLength','ConvexArea');

    % consider only upper half of the slice
    centroids = cat(1,regions.Centroid);
    y_lim =  size(I3,1)/2;
    regions = regions(centroids(:,2)<y_lim);

    % middle region only
    x_lim = [size(I3,2)*1/3, size(I3,2)*2/3];
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

    
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(I3),hold on

    for idx = 1:2
        col = 'b';
        h = rectangle('Position',dg_regions(idx).BoundingBox);
        set(h,'EdgeColor',col,'LineWidth',1);
        plot(dg_centroids(idx,1),dg_centroids(idx,2),'r.')
        1;
    end
    
    fname = strcat(file(1:end-11),'_dg_regions.tif');
    saveas(fig1,fullfile(roi_folder,fname));
    close(fig1);
    
    %% Calculate offset
    if base_compare
        offset = base_dg_centroids - dg_centroids;
    else
        offset = [0 0; 0 0];
    end
    
    %% Save DG regions
    dg_fname = strcat(file(1:end-11),'_dg_regions.mat');
    save(fullfile(roi_folder,dg_fname),'dg_regions','dg_centroids','offset');
    
end
    