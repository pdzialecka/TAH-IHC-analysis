function [] = analyse_data(files,load_rois)
    %% Analyse IHC data across selected ROIs
    % @author: pdzialecka
    
    % Variables to explore:
    % * pixel_thresh_factor (default: 1)
    % * min_con_pixels (default: 25 at 20x, 100 at 10x)
    % * k3 (spatial filter variable) - currently no difference after filter?
    % * magnification (may affect other parameters)
    
    %% Default input options
    if ~exist('load_rois','var')
        load_rois = 1;
    end
        
    %% H DAB colormaps
    % for easier visualisation
    [h_colormap,dab_colormap] = create_hdab_colormaps();
    
    %% ROIs
    [~,roi_fnames,roi_no] = get_roi_list();
    roi_order_no = 1:roi_no;
    
    %%
    for idx = 1:length(files)
        
        file = files(idx).name;
        folder = files(idx).folder;
        
        %% ROI folder
        data_folder = (folder);
        [roi_folder,mouse_name] = find_roi_folder(data_folder);

        %% ROI images folder
        roi_img_folder = fullfile(fileparts(fileparts(folder)),'ROI_images',mouse_name);

        if ~exist(roi_img_folder)
            mkdir(roi_img_folder);
        end

        %% Image type & threshold
        if contains(file,'moc23')
            image_type = 'moc23';
%             pixel_thresh = 0.65; % 160;

        elseif contains(file,'cfos')
            image_type = 'cfos';
%             pixel_thresh = 0.2; % decrease to detect only v dark cells

        elseif contains(file,'GFAP')
            image_type = 'GFAP';
%             pixel_thresh = 0.5;

        elseif contains(file,'Iba1')
            image_type = 'Iba1';
%             pixel_thresh = 0.35;
        end
        
        pixel_thresh = get_antibody_threshold(image_type);
        
    
        %% Results folder
        results_folder = fullfile((fileparts(fileparts(folder))),'Results');

        % create subfolder with mouse name
        [~,mouse_name] = fileparts((folder));
        results_folder = fullfile(results_folder,mouse_name);

        if ~exist(results_folder)
            mkdir(results_folder);
        end

        %% Load deconvolved images
        file_path = fullfile(folder,file);
        [h_image,dab_image] = load_deconvolved_images(file_path);
        
        %% Find ROIs
        roi_not_found = 0;

        for roi_idx = 1:roi_no
    % %         fname = strcat(file(1:end-11),'_roi_',roi_fnames{roi_idx},'.mat');
            fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'.mat');

            %% Load preselected ROI
            if load_rois
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

        end

        %% Else select ROIs
        if ~load_rois || roi_not_found
            [rois_x,rois_y] = select_roi(files(idx));
        end

        %% Analysis per ROI
        for roi_idx = 1:roi_no
            
            %%
            fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'');

            %% Extract ROIs + denoise
            dab_image_roi = dab_image(rois_x{roi_idx},rois_y{roi_idx});
            h_image_roi = h_image(rois_x{roi_idx},rois_y{roi_idx});

            dab_image_ = dab_image_roi;
            h_image_ = h_image_roi;


            % Spatial smoothing
            spatial_avg = 1;
            if spatial_avg
                k3 = 3;

                % gaussian and box filters comparable
                % box blur applied many times = approximation of Gaussian blur
                % gaussian filter requires sigma estimation and potentially
                % larger filter size

                box_kernel = ones(k3,k3)*1/(k3*k3);
                h_image_f = imfilter(h_image_,box_kernel,'replicate');
                dab_image_f = imfilter(dab_image_,box_kernel,'replicate');

%                 sigma = 10;
%                 h_image_f = imgaussfilt(h_image_,sigma,'FilterSize',k3,'Padding','replicate');
%                 dab_image_f = imgaussfilt(dab_image_,sigma,'FilterSize',k3,'Padding','replicate');

            %     dab_image_f = imsharpen(image_dab_,'Radius',9,'Amount',1);
            
            
                % from Nir's scripts: contrast enhancing + adaptive filtering
                % doesn't improve moc23 images
%                 n_tiles = 8;
%                 dab_image_ = adapthisteq(dab_image_, 'NumTiles',[n_tiles n_tiles],'ClipLimit',0.001);   % Default: NumTiles=[8 8], ClipLimit=0.01
% 
%                 flt_loc_pixs_size = k3; % 5
%                 [dab_image_f, noise] = wiener2(dab_image_, [flt_loc_pixs_size flt_loc_pixs_size]);
%                 disp([' estimated noise level before filtering: ', num2str(noise)]);

               
                fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
                subplot(121),imshow(h_image_),colormap(h_colormap),title('Original')
                subplot(122),imshow(h_image_f),colormap(h_colormap),title('Filtered')

                fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
                subplot(121),imshow(dab_image_),colormap(dab_colormap),title('Original')
                subplot(122),imshow(dab_image_f),colormap(dab_colormap),title('Filtered')
                
                fname1 = strcat(fname,'_1_H_filtered.tif');
                saveas(fig1,fullfile(results_folder,fname1));
                close(fig1);

                fname2 = strcat(fname,'_2_DAB_filtered.tif');
                saveas(fig2,fullfile(results_folder,fname2));
                close(fig2);
                

                dab_image_ = dab_image_f;
                h_image_ = h_image_f;
            end

            %% Threshold
            % white background: 255
            % areas of interest: lower pixel values
            
            % constant threshold needed for all images with the same
            % antibody (or just one image?)

            % method 1: manually chosen threshold
%             dab_roi_mask = dab_image_ < pixel_thresh;
            dab_roi_mask = ~imbinarize(dab_image_,pixel_thresh);
            
            % method 2: automated Otsu's threshold, can use to find approx value
%             pixel_thresh = graythresh(dab_image_);
%             pixel_thresh_factor = 1;
%             pixel_thresh = pixel_thresh*pixel_thresh_factor
%             dab_roi_mask = ~imbinarize(dab_image_,pixel_thresh);
            

%             figure,imshow(dab_image_),colormap(dab_colormap)
            
            dab_image_mask = imoverlay(dab_image_,dab_roi_mask,[1,0,0]); % red
%             figure,imshow(dab_image_mask)

            %%
            % Mask filtering from Nir's scripts
% 
            % Morphological opening using a disc kernel (smoothes contours, breaks
            % narrow bridges, eliminates small islands)
%             strel_disk_size = 1;
%             se = strel('disk',strel_disk_size);
%             dab_roi_mask_2 = imopen(dab_roi_mask, se);
%             disp([' Morphological opening. Disc size used = ', num2str(strel_disk_size)]);
%             figure,imshow(dab_roi_mask_2)

            % Remove small regions of connected pixels
            min_con_pixels = 25;
            connectivity = 4; % default: 4
            dab_roi_mask_f = bwareaopen(dab_roi_mask,min_con_pixels,connectivity);
            dab_image_mask_f = imoverlay(dab_image_,dab_roi_mask_f,[1,0,0]);

            fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(121),imshow(dab_image_mask),title('Mask')
            subplot(122),imshow(dab_image_mask_f),title('Mask + filter')
            
            fname3 = strcat(fname,'_3_DAB_masks.tif');
            saveas(fig3,fullfile(results_folder,fname3));
            close(fig3);
            
            % alternatively, plot mask and overlay image
%             fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
%             subplot(121),imshow(dab_roi_mask_f),title('Mask')
%             subplot(122),imshow(dab_image_mask_f),title('Image + Filtered Mask')
% 
%             fname3 = strcat(fname,'_3_H_masks.tif');
%             saveas(fig3,fullfile(results_folder,fname3));
%             close(fig3);

            
            mask_name1 = strcat(fname,'_thresh_mask');
            save(fullfile(results_folder,mask_name1),'dab_roi_mask');
            fprintf('Mask %s saved\n',mask_name1);

            
            dab_roi_mask = dab_roi_mask_f;

            mask_name2 = strcat(fname,'_thresh_filt_mask');
            save(fullfile(results_folder,mask_name2),'dab_roi_mask');
            fprintf('Mask %s saved\n',mask_name2);
            
            %% Save ROI images
            roi_image = cat(3,h_image_,dab_image_);
%             magnification = 10; % load from inside roi_file later
%             img_fname = fullfile(roi_img_folder,strcat(fname,'_',num2str(magnification),'x'));
            
            img_fname = fullfile(roi_img_folder,strcat(fname));
            save(img_fname,'roi_image');

            %% Density (% area covered)
            positive_pixels = sum(dab_roi_mask,[1,2]);
            total_pixels = size(dab_roi_mask,1)*size(dab_roi_mask,2);

            density = round(positive_pixels/total_pixels*100,1);
            fprintf('Area covered for %s ROI = %1.1f %% \n', fname, density)
            
            %% Save results
            results.fname = fname;
            results.thresh_pixel = pixel_thresh;
            results.density = density;
            results.positive_pixels = positive_pixels;
            results.total_pixels = total_pixels;
            
            fname_r = strcat(fname,'_results.mat');
            save(fullfile(results_folder,fname_r),'results');
            fprintf('Results %s saved\n',fname_r);
            
            fprintf('\n')
            
        end
    end
end


