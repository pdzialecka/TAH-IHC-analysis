function [] = analyse_data(files,load_rois)
    %% Analyse IHC data across selected ROIs
    % @author: pdzialecka
    
    %% Default input options
    if ~exist('load_rois','var')
        load_rois = 0;
    end
    
%     if ~exist('magnification','var')
%         magnification = 10; % 10 or 20x
%     end
    
    %% H DAB colormaps
    % for easier visualisation
    dab_col = (1-[0.26814753,0.57031375,0.77642715]);
    h_col = (1-[0.6500286,0.704031,0.2860126]);

    dab_colormap = [(linspace(dab_col(1),1,256)'),(linspace(dab_col(2),1,256)'),(linspace(dab_col(3),1,256)')];
    h_colormap = [(linspace(h_col(1),1,256)'),(linspace(h_col(2),1,256)'),(linspace(h_col(3),1,256)')];
    
    %% ROIs
%     roi_names = {'Right Hipp','Left Hipp','Right Cortex','Left Cortex','Control Area'};
%     roi_no = length(roi_names);
%     roi_fnames = {'hipp_R','hipp_L','cortex_R','cortex_L','control'};
%     roi_order_no = 1:roi_no;
%     
%     rois_x = {};
%     rois_y = {};
    
    %%
    idx = 1;
    file = files(idx).name;
    folder = files(idx).folder;
        
    %% Image type
    if contains(file,'moc23')
        image_type = 'moc23';
        pixel_thresh = 160;
        
    elseif contains(file,'cfos')
        image_type = 'cfos';
        pixel_thresh = 150;
        
    elseif contains(file,'GFAP')
        image_type = 'GFAP';
        pixel_thresh = 150;
        
    elseif contains(file,'Iba1')
        image_type = 'Iba1';
        pixel_thresh = 150;
        
    end
    
    %% ROI folder
%     roi_folder = fullfile(fileparts(fileparts(folder)),'ROIs');
% 
%     if ~exist(roi_folder)
%         mkdir(roi_folder);
%     end
    
    %% Load deconvolved images
    image = read_file(fullfile(folder,file));
    h_image = image(:,:,1);
    dab_image = image(:,:,2);
%     image_res = image(:,:,3);

    %% Rotate images
    h_image = imrotate(h_image,-90);
    dab_image = imrotate(dab_image,-90);
    
    %% Find ROIs
    roi_not_found = 0;
%     roi_size = round(size(dab_image)/magnification);

    for roi_idx = 1:roi_no
% %         fname = strcat(file(1:end-11),'_roi_',roi_fnames{roi_idx},'.mat');
        fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'.mat');
%         
        %% ROI folder
        roi_folder = fullfile(fileparts(fileparts(folder)),'ROIs');

        % create subfolder with mouse name
        mouse_name = file(1:9);
        roi_folder = fullfile(roi_folder,mouse_name);

        if ~exist(roi_folder)
            mkdir(roi_folder);
        end

        %% Load preselected ROI
        if load_rois
%             fname = strcat(file(1:end-11),'_roi_',roi_fnames{roi_idx},'.mat');
            roi_fname = fullfile(roi_folder,fname);

            if exist(roi_fname)
                roi_file = load(roi_fname).roi;
                rois_x{roi_idx} = roi_file.x;
                rois_y{roi_idx} = roi_file.y;
                fprintf('Successfully loaded %s ROI\n',fname);

            else
                roi_not_found = 1; % flag for at least one missing ROI
                fprintf('ROI file %s not found; ROI selection required\n',fname)
            end
        end

        %% Select ROIs
%         if ~load_rois || roi_not_found
%             [rois_x,rois_y] = select_roi(files(idx));
%             
%             roi_accepted = 0;
% 
%             while ~roi_accepted
%                 %% Select ROI on H image
%                 fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
%                 imshow(h_image),colormap(h_colormap)
% 
%                 title(sprintf('Select %s ROI',roi_names{roi_idx})) 
% 
%                 % roi = drawrectangle();
%                 % 
%                 % roi.Position = round(roi.Position);
%                 % roi_x = roi.Position(2):roi.Position(2)+roi.Position(4);
%                 % roi_y = roi.Position(1):roi.Position(1)+roi.Position(3);
% 
%                 roi_point = drawpoint();
% 
%                 %% Diplay ROI on DAB image
%                 roi_point.Position = round(roi_point.Position);
%                 roi_x_1 = round(roi_point.Position(2)-roi_size(1)/2);
%                 roi_y_1 = round(roi_point.Position(1)-roi_size(2)/2);
% 
%                 roi_x = roi_x_1:roi_x_1+roi_size(1);
%                 roi_y = roi_y_1:roi_y_1+roi_size(2);
% 
%                 roi_rect = drawrectangle('Position',[roi_y_1,roi_x_1,roi_size(2),roi_size(1)]);
%                 title(sprintf('%s ROI',roi_names{roi_idx}));
% 
%                 fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
%                 imshow(dab_image),colormap(dab_colormap)
%                 roi_rect = drawrectangle('Position',[roi_y_1,roi_x_1,roi_size(2),roi_size(1)]);
%                 title(sprintf('%s ROI',roi_names{roi_idx}));
% 
% 
%                 %% Display ROI selected
%                 h_image_roi = h_image(roi_x,roi_y);
%                 dab_image_roi = dab_image(roi_x,roi_y);
% 
%                 fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
%                 ax(1) = subplot(121); imshow(h_image_roi)
%                 ax(2) = subplot(122); imshow(dab_image_roi)
%                 colormap(ax(1),h_colormap)
%                 colormap(ax(2),dab_colormap)
% 
%                 %% Accept or reject ROI
%                 answer = questdlg('Do you want to accept this ROI?',...
%                     'Accept ROI?','Yes','No','Yes');
% 
%                 % Handle response
%                 switch answer
%                     case 'Yes'
%                         roi_accepted = 1;
%                         rois_x{roi_idx} = roi_x;
%                         rois_y{roi_idx} = roi_y;
% 
%                     case 'No'
%                         roi_accepted = 0;
%                         close(fig1);
%                         close(fig2);
%                         close(fig3);
%                 end
% 
%                 %% Save accepted ROI coordinates
%                 if roi_accepted
% %                     fname = strcat(file(1:end-11),'_roi_',roi_fnames{roi_idx});
% 
%                     % figures
%                     fname1 = strcat(fname,'_1_H_full.tif');
%                     saveas(fig1,fullfile(roi_folder,fname1));
% 
%                     fname2 = strcat(fname,'_2_DAB_full.tif');
%                     saveas(fig2,fullfile(roi_folder,fname2));
% 
%                     fname3 = strcat(fname,'_3_H_DAB_',num2str(magnification),'x.tif');
%                     saveas(fig3,fullfile(roi_folder,fname3));
%                     
%                     close(fig1);
%                     close(fig2);
%                     close(fig3);
%                     
%                     
%                     % roi details
%                     roi.name = roi_names{roi_idx};
%                     roi.fname = roi_fnames{roi_idx};
%                     roi.x = roi_x;
%                     roi.y = roi_y;
% 
%                     fname = strcat(file(1:end-11),'_roi_',roi_fnames{roi_idx},'.mat');
%                     save(fullfile(roi_folder,fname),'roi');
%                 end
%             end


    end
    
    %% Else select ROIs
    if ~load_rois || roi_not_found
        [rois_x,rois_y] = select_roi(files(idx));
    end
    
    %% Analysis per ROI
    for roi_idx = 1:roi_no
        use_roi = 1;
        
        dab_image_roi = dab_image(rois_x{roi_idx},rois_y{roi_idx});
        h_image_roi = h_image(rois_x{roi_idx},rois_y{roi_idx});

        if use_roi
            dab_image_ = dab_image_roi;
            h_image_ = h_image_roi;
        else
            dab_image_ = dab_image;
            h_image_ = h_image;
        end


        dab_image_orig = dab_image_;
        h_image_orig = h_image_;


        % Spatial smoothing
        spatial_avg = 1;
        if spatial_avg
            k3 = 3;

            % both filters comparable
            % box blur applied many times = approximation of Gaussian blur

    %         box_kernel = ones(k3,k3)*1/(k3*k3);
    %         h_image_f = imfilter(h_image_,box_kernel,'replicate');
    %         dab_image_f = imfilter(dab_image_,box_kernel,'replicate');

            h_image_f = imgaussfilt(h_image_,'FilterSize',k3,'Padding','replicate');
            dab_image_f = imgaussfilt(dab_image_,'FilterSize',k3,'Padding','replicate');

        %     dab_image_f = imsharpen(image_dab_,'Radius',9,'Amount',1);

            figure('units','normalized','outerposition',[0 0 1 1])
            subplot(121),imshow(h_image_),colormap(h_colormap)
            subplot(122),imshow(h_image_f),colormap(h_colormap)

            figure('units','normalized','outerposition',[0 0 1 1])
            subplot(121),imshow(dab_image_),colormap(dab_colormap)
            subplot(122),imshow(dab_image_f),colormap(dab_colormap)

    %         dab_image_orig = dab_image_;
            dab_image_ = dab_image_f;

    %         h_image_orig = h_image_;
            h_image_ = h_image_f;
        end

        %% Threshold
        % white background: 255
        % areas of interest: lower pixel values

        dab_roi_mask = dab_image_<pixel_thresh;
        dab_image_mask = imoverlay(dab_image_orig,dab_roi_mask,'red');
        figure,imshow(dab_image_mask)

        %% Density (% area covered)
        positive_pixels = sum(dab_roi_mask,[1,2]);
        total_pixels = size(dab_roi_mask,1)*size(dab_roi_mask,2);

        density = positive_pixels/total_pixels*100

        %%
        disp('* BW conditioning');
        figure,imshow(dab_roi_mask)

        % Morphological opening using a disc kernel (smoothes contours, breaks
        % narrow bridges, eliminates small islands)
        strel_disk_size = 1;
        se = strel('disk',strel_disk_size);
        dab_roi_mask_2 = imopen(dab_roi_mask, se);
        disp([' Morphological opening. Disc size used = ', num2str(strel_disk_size)]);
        figure,imshow(dab_roi_mask_2)

        % Remove samll regions of connected pixels
        min_con_pixels = 30;
        connectivity = 4; % 4
        dab_roi_mask_3 = bwareaopen(dab_roi_mask_2, min_con_pixels, connectivity);
        disp([' Removed connected regions smaller than ', num2str(min_con_pixels),' with connectivity of ', num2str(connectivity), ' pixels']);

        figure,imshow(dab_roi_mask_3)
        dab_image_mask_2 = imoverlay(dab_image_orig,dab_roi_mask_3,'red');
        figure,imshow(dab_image_mask_2)

    end
    
end
