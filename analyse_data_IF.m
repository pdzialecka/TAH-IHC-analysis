function [] = analyse_data_IF(files,close_figs)
    %% Analyse IF data across selected ROIs: simple colocalisation
    % @author: pdzialecka
    
    %% Default input options
%     if ~exist('load_rois','var')
%         load_rois = 1;
%     end
    
    if ~exist('close_figs','var')
        close_figs = 1;
    end
    
    %% G&B colormaps
    [r_colormap,g_colormap,~] = create_rgb_colormaps();
    
    %% Plot setts
    pixel_size = 0.2840; % from image info; else,1um = 3.5271 pixels
%     color_map = [0,0,1];
    transparency = 0.5;
    
    %% ROIs
%     [~,roi_fnames,roi_no] = get_roi_list();
%     roi_order_no = 1:roi_no;

    %%
    for idx = 1:length(files)
        
        file = files(idx).name;
        folder = files(idx).folder;
        fname = file(1:end-4);
        
        %% ROI folder
%         data_folder = (folder);
%         [roi_folder,mouse_name] = find_roi_folder(data_folder);

        %% ROI images folder
%         roi_img_folder = fullfile(fileparts(fileparts(folder)),'ROI_images',mouse_name);
% 
%         if ~exist(roi_img_folder)
%             mkdir(roi_img_folder);
%         end

        %% Image type: threshold + min particle size
%         idxs_ = strfind(file,'_');
%         image_type = file(idxs_(1)+1:idxs_(2)-1);
%         img_type = find_img_type(file);
%         
% %         if contains(file,'ab_4g8')
% %             image_type = 'ab_4g8';
% % 
% %         elseif contains(file,'cfos')
% %             image_type = 'cfos';
% %             do_watershed = 1;
% % 
% %         elseif contains(file,'GFAP')
% %             image_type = 'GFAP';
% % 
% %         elseif contains(file,'Iba1')
% %             image_type = 'Iba1';
% %         end
%         
%         [pixel_thresh,min_size,max_size,do_watershed] = get_antibody_threshold(img_type);
        
    
        %% Results folder
        results_folder = fullfile(fileparts(fileparts(fileparts(folder))),'Results');

        % create subfolder with mouse name
        [~,mouse_name] = fileparts(fileparts((folder)));
        results_folder = fullfile(results_folder,mouse_name);

        if ~exist(results_folder)
            mkdir(results_folder);
        end

        %% Load deconvolved images
%         file_path = fullfile(folder,file);
%         [h_image,dab_image] = load_deconvolved_images(file_path);
        
        %% Find ROIs
%         roi_not_found = 0;
%         rois_x = {};
%         rois_y = {};
%         rois_n_area = [];
%         rois_size_um = {};

        %% Load preselected ROI
%         for roi_idx = 1:roi_no
%     % %         fname = strcat(file(1:end-11),'_roi_',roi_fnames{roi_idx},'.mat');
%             fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'.mat');
% 
%             if load_rois
%                 roi_fname = fullfile(roi_folder,fname);
% 
%                 if exist(roi_fname)
%                     roi_file = load(roi_fname).roi;
%                     rois_x{roi_idx} = roi_file.x;
%                     rois_y{roi_idx} = roi_file.y;
% %                     rois_n_area(roi_idx) = roi_file.slice_area_norm;
%                     rois_size_um{roi_idx} = roi_file.size_um;
% %                     fprintf('Successfully loaded %s ROI\n',fname);
% 
%                 else
%                     roi_not_found = 1; % flag for at least one missing ROI
%                     fprintf('ROI file %s not found; ROI selection required\n',fname)
%                 end
%             end
% 
%         end
% 
%         %% Else select ROIs
%         if ~load_rois || roi_not_found
%             [rois_x,rois_y] = select_roi(files(idx));
%         end
%         
%         %% Load slice mask
%         sm_file = dir(fullfile(roi_folder,strcat('*',img_type,'*slice_mask.mat')));
%         slice_mask = load(fullfile(sm_file.folder,sm_file.name)).slice_mask;
        
        %% ROIs
%         if strcmp(img_type,'ki67') || strcmp(img_type,'dcx') || strcmp(img_type,'sox2')
%             roi_idxs = 1:2; % DG only
%         else
%             roi_idxs = 1:roi_no;
%         end

        %% Load reference image
        ref_file = dir(fullfile(fileparts(fileparts(fileparts(fileparts(fileparts(files(1).folder))))),'Reference_images','IF',strcat('*.tif')));
        ref_fname = fullfile(ref_file(1).folder,ref_file(1).name);
        ref_image = read_file(ref_fname,1);
        
        ref_dapi_img = ref_image(:,:,1);
        ref_iba1_img = ref_image(:,:,2);
        ref_ab_4g8_img = ref_image(:,:,3);
        
        %% Analysis per ROI
%         for roi_idx = 1:roi_no
            
            %%
%             fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'');            
            
            %% Extract ROIs from main image (slower)
%             dab_image_roi = dab_image(rois_x{roi_idx},rois_y{roi_idx});
%             h_image_roi = h_image(rois_x{roi_idx},rois_y{roi_idx});
% 
        %% Load ROI image
        file_path = fullfile(folder,file);
        image = read_file(file_path);
        
        %%
        img_size = size(image,[1,2]);
        orig_img_size = img_size;
        
        if img_size(2) == 1400
            roi_size = [2800,1400]; % better resolution
        elseif img_size(2) == 700
            roi_size = [1400,700];
        end
        
        %% Cut image to fixed size
        size_to_cut = img_size(1)-roi_size(1);
        
        if mod(size_to_cut,2) % odd
            cut_size = [ceil(size_to_cut/2) floor(size_to_cut/2)];
        else
            cut_size = [size_to_cut/2 size_to_cut/2];
        end
        
        % keep the middle pixels
        cut_idxs = [1:cut_size(1), img_size(1)-(cut_size(2)-1):img_size(1)];
        image(cut_idxs,:,:) = [];
        
        img_size = size(image,[1,2]);
        
        if ~all(img_size==roi_size)
            fprintf('WARNING: Incorrect ROI size detected (%s)\n',file);
        end
        
        %% Resize too small images for convenience
        % worse quality of images
        if img_size(2) == 700
            image = imresize(image,2);
        end

        %% Extract subimages and rotate
        dapi_img = imrotate(image(:,:,1),90);
        iba1_img = imrotate(image(:,:,2),90);
        ab_4g8_img = imrotate(image(:,:,3),90);
        
        %% Save final ROIs
        img_fname = fullfile(fileparts(folder),strcat(fname,'.tif'));
        imwrite(dapi_img,img_fname,'Compression','none','WriteMode','overwrite');
        imwrite(iba1_img,img_fname,'Compression','none','WriteMode','append');
        imwrite(ab_4g8_img,img_fname,'Compression','none','WriteMode','append');

        %% Smooth images
        spatial_avg = 1;
        if spatial_avg
            k3 = 9;

            box_kernel = ones(k3,k3)*1/(k3*k3);
            
            dapi_img_f = imfilter(dapi_img,box_kernel,'replicate');
            iba1_img_f = imfilter(iba1_img,box_kernel,'replicate');
            ab_4g8_img_f = imfilter(ab_4g8_img,box_kernel,'replicate');
            
            fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(211),imshow(iba1_img),colormap(r_colormap),title('Original')
            subplot(212),imshow(iba1_img_f),colormap(r_colormap),title('Filtered')

            fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(211),imshow(ab_4g8_img),colormap(g_colormap),title('Original')
            subplot(212),imshow(ab_4g8_img_f),colormap(g_colormap),title('Filtered')

            fname1 = strcat(fname,'_1_iba1_filtered.tif');
            saveas(fig1,fullfile(results_folder,fname1));
            if close_figs
                close(fig1);
            end

            fname2 = strcat(fname,'_2_4g8_filtered.tif');
            saveas(fig2,fullfile(results_folder,fname2));
            if close_figs
                close(fig2);
            end
% 
            dapi_img = dapi_img_f;
            iba1_img = iba1_img_f;
            ab_4g8_img = ab_4g8_img_f;
        end
        
        %% Normalise images
        normalise_imgs = 0;
        
        if normalise_imgs
            dapi_img_n = imadjust(dapi_img);
            iba1_img_n = imadjust(iba1_img);
            ab_4g8_img_n = imadjust(ab_4g8_img);

            fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(211),imshow(iba1_img),colormap(r_colormap),title('Original filtered')
            subplot(212),imshow(iba1_img_n),colormap(r_colormap),title('Normalised')

            fig4 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(211),imshow(ab_4g8_img),colormap(g_colormap),title('Original filtered')
            subplot(212),imshow(ab_4g8_img_n),colormap(g_colormap),title('Normalised')

            fname3 = strcat(fname,'_3_iba1_normalised.tif');
            saveas(fig3,fullfile(results_folder,fname3));
            if close_figs
                close(fig3);
            end

            fname4 = strcat(fname,'_4_4g8_normalised.tif');
            saveas(fig4,fullfile(results_folder,fname4));
            if close_figs
                close(fig4);
            end
        end
        
        %% Normalise image brightness
        normalise_brightness = 1;
        
        if normalise_brightness
            iba1_img_n = imhistmatch(iba1_img,ref_iba1_img);
            ab_4g8_img_n = imhistmatch(ab_4g8_img,ref_ab_4g8_img);

            roi_b_img_folder = fullfile(fileparts(fileparts(fileparts(folder))),'ROI_images_norm',mouse_name);
            if ~exist(roi_b_img_folder,'dir')
                mkdir(roi_b_img_folder);
            end

            img_fname = fullfile(roi_b_img_folder,strcat(fname,'.tif'));

            imwrite(iba1_img_n,img_fname,'Compression','none','WriteMode','overwrite');
            imwrite(ab_4g8_img_n,img_fname,'Compression','none','WriteMode','append');

            
            fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(211),imshow(iba1_img),colormap(r_colormap),title('Original filtered')
            subplot(212),imshow(iba1_img_n),colormap(r_colormap),title('Normalised')

            fig4 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(211),imshow(ab_4g8_img),colormap(g_colormap),title('Original filtered')
            subplot(212),imshow(ab_4g8_img_n),colormap(g_colormap),title('Normalised')

            fname3 = strcat(fname,'_3_iba1_bc.tif');
            saveas(fig3,fullfile(results_folder,fname3));
            if close_figs
                close(fig3);
            end

            fname4 = strcat(fname,'_4_4g8_bc.tif');
            saveas(fig4,fullfile(results_folder,fname4));
            if close_figs
                close(fig4);
            end
            
%             iba1_img = iba1_img_n;
%             ab_4g8_img = ab_4g8_img_n;
            
        end
        
        %% Create microglia and ab masks
        if normalise_imgs
            iba1_thresh = 0.8*255;
            ab_4g8_thresh = 0.8*255;
            
            microglia_mask = iba1_img_n>iba1_thresh;
            ab_mask = ab_4g8_img_n>ab_4g8_thresh;
            
        elseif normalise_brightness
            iba1_thresh = 0.4; % graythresh(iba1_img)
            microglia_mask = imbinarize(iba1_img_n,iba1_thresh);
            
            ab_4g8_thresh = 0.6;
            ab_mask = imbinarize(ab_4g8_img_n,ab_4g8_thresh);
            
        else % variable threshold
            iba1_thresh = 0.4; % graythresh(iba1_img)
            microglia_mask = imbinarize(iba1_img,iba1_thresh);
            
            ab_4g8_thresh = 0.5;
            ab_mask = imbinarize(ab_4g8_img,ab_4g8_thresh);
        end
        
        
        % remove too small and too big elements
        [~,min_size_1,max_size_1] = get_antibody_threshold('iba1');
        microglia_mask_f = bwpropfilt(microglia_mask,'EquivDiameter',[min_size_1/pixel_size,max_size_1/pixel_size]);

        [~,min_size_2,max_size_2] = get_antibody_threshold('moc23');
        ab_mask_f = bwpropfilt(ab_mask,'EquivDiameter',[min_size_2/pixel_size,max_size_2/pixel_size]);

        %% Visualise the masks
        iba1_image_mask = labeloverlay(iba1_img,microglia_mask,...
            'Colormap',[0,0,1],'Transparency',transparency);
        iba1_image_mask_f = labeloverlay(iba1_img,microglia_mask_f,...
            'Colormap',[0,0,1],'Transparency',transparency);

        ab_4g8_image_mask = labeloverlay(ab_4g8_img,ab_mask,...
            'Colormap',[0,0,1],'Transparency',transparency);
        ab_4g8_image_mask_f = labeloverlay(ab_4g8_img,ab_mask_f,...
            'Colormap',[0,0,1],'Transparency',transparency);
        
        fig5 = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(211),imshow(iba1_image_mask),title('Threshold')
        subplot(212),imshow(iba1_image_mask_f),title('Threshold + filtered')

        fig6 = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(211),imshow(ab_4g8_image_mask),title('Threshold')
        subplot(212),imshow(ab_4g8_image_mask_f),title('Threshold + filtered')

        fname5 = strcat(fname,'_5_iba1_mask.tif');
        saveas(fig5,fullfile(results_folder,fname5));
        if close_figs
            close(fig5);
        end

        fname6 = strcat(fname,'_6_4g8_mask.tif');
        saveas(fig6,fullfile(results_folder,fname6));
        if close_figs
            close(fig6);
        end
        
        microglia_mask = microglia_mask_f;
        ab_mask = ab_mask_f;
        
        %% Create colocolised mask
        colocalised_mask = microglia_mask & ab_mask;
        
        composite_img = imfuse(iba1_img,ab_4g8_img,'ColorChannels',[1 2 0]);
        composite_image_mask = labeloverlay(composite_img,colocalised_mask,...
            'Colormap',[0,0,1],'Transparency',0.2);
       
        fig7 = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(211),imshow(composite_img),title('Composite image')
        subplot(212),imshow(composite_image_mask),title('Composite image + overlap mask')
        
        
        fname7 = strcat(fname,'_colocalised_mask.tif');
        saveas(fig7,fullfile(results_folder,fname7));
        if close_figs
            close(fig7);
        end
        
        %% Save masks
        mask_name = strcat(fname,'_masks');
        save(fullfile(results_folder,mask_name),'microglia_mask','ab_mask','colocalised_mask');
        
        %% Compute % microglia that are ab+
        microglia_area = sum(microglia_mask,[1,2]);
        ab_area = sum(ab_mask,[1,2]);
        microglia_ab_pos_area = sum(colocalised_mask,[1,2]);
        
        microglia_ab_ratio = round((microglia_ab_pos_area/microglia_area)*100,2);
        ab_microglia_ratio = round((ab_area/microglia_area)*100,2);
        fprintf('Microglia ab+ = %1.2f%% for %s\n',microglia_ab_ratio,file);
        fprintf('Ab microglia+ = %1.2f%% for %s\n',ab_microglia_ratio,file);

        %% Save results
        results = [];
        results.orig_img_size = orig_img_size;
        results.img_size = img_size;
        results.total_area = img_size(1)*img_size(2);
        results.microglia_area = microglia_area;
        results.ab_area = ab_area;
        results.microglia_ab_ratio = microglia_ab_ratio;
        results.ab_microglia_ratio = ab_microglia_ratio;
        
%             results.density = density;
%             results.positive_pixels = positive_pixels;
%             results.total_pixels = total_pixels;
%             results.particles_found = particles_found;
%             results.particle_no = particle_no;
%             
%             if strcmp(img_type,'cfos')
%                 results.pos_particle_ratio = pos_particle_ratio;
%             end
%             
        fname_r = strcat(file(1:end-4),'_results.mat');
        save(fullfile(results_folder,fname_r),'results');
        fprintf('Results %s saved\n',fname_r);

        
        %%
%             %% Normalise brightness based on the reference image
%             if strcmp(img_type,'ki67') || strcmp(img_type,'cfos')
%                 normalise_brightness = 0;
%             else
%                 normalise_brightness = 1;
%             end
%             
%             if normalise_brightness
%                 h_image_roi = imhistmatch(h_image_roi,h_image_ref);
%     %             if ~contains(fname,'cfos')
%                 dab_image_roi = imhistmatch(dab_image_roi,dab_image_ref);
%     %             end
% 
%                 roi_b_img_folder = fullfile(fileparts(fileparts(folder)),'ROI_images_norm',mouse_name);
%                 if ~exist(roi_b_img_folder,'dir')
%                     mkdir(roi_b_img_folder);
%                 end
%     
% %                 roi_image = cat(2,h_image_roi,dab_image_roi);
%                 img_fname = fullfile(roi_b_img_folder,strcat(fname,'.tif'));
% 
%                 imwrite(h_image_roi,img_fname,'Compression','none','WriteMode','overwrite');
%                 imwrite(dab_image_roi,img_fname,'Compression','none','WriteMode','append');
%                
%             end
%             
%             %% Find slice mask for ROI
%             roi_x = rois_x{roi_idx};
%             roi_y = rois_y{roi_idx};
%             slice_mask_roi = slice_mask(roi_y,roi_x);
%             roi_size_um = rois_size_um{roi_idx};
%                         
%             %% Display zoomed in ROI w slice mask
%             h_image_roi_smask = labeloverlay(h_image_roi,~slice_mask_roi,...
%                 'Colormap',[0,0,1],'Transparency',0.7);
%             dab_image_roi_smask = labeloverlay(dab_image_roi,~slice_mask_roi,...
%                 'Colormap',[0,0,1],'Transparency',0.7);
% 
%             fig0 = figure('units','normalized','outerposition',[0 0 1 1]);
%             ax(1) = subplot(211); imshow(h_image_roi_smask)
%             ax(2) = subplot(212); imshow(dab_image_roi_smask)
%             colormap(ax(1),h_colormap)
%             colormap(ax(2),dab_colormap)
% 
%             
%             fname_ = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx});
%             fname0 = strcat(fname_,'_3_H_DAB_',num2str(roi_size_um(1)),'x',num2str(roi_size_um(2)),'um_slice_mask.tif');
%             saveas(fig0,fullfile(results_folder,fname0));
%             close(fig0);
%         
%             %% Calculate slice area per ROI
%             slice_area = sum(slice_mask_roi,[1,2]);
%             total_area = size(h_image_roi,1)*size(h_image_roi,2);
%             slice_area_norm = slice_area/total_area; % as % of total area
% 
% %             roi.slice_area = slice_area;
% %             roi.total_area = total_area;
% %             roi.slice_area_norm = slice_area_norm;
% 
%             rois_n_area(roi_idx) = slice_area_norm;
%             tot_n_area = rois_n_area(roi_idx);
%      
%             %% Denoise images
%             % This filter gets rid of most of the random pepper noise
%             
%             % Gaussian and box filters comparable:
%             % box blur applied many times = approximation of Gaussian blur
%             % gaussian f requires sigma estimation and potentially larger filter size
%             
% %             h_image_ = h_image_roi;
% %             dab_image_ = dab_image_roi;
% 
%             % Spatial smoothing
%             spatial_avg = 1;
%             if spatial_avg
%                 k3 = 3;
% 
%                 box_kernel = ones(k3,k3)*1/(k3*k3);
%                 h_image_f = imfilter(h_image_roi,box_kernel,'replicate');
%                 dab_image_f = imfilter(dab_image_roi,box_kernel,'replicate');
% 
% %                 sigma = 10;
% %                 h_image_f = imgaussfilt(h_image_,sigma,'FilterSize',k3,'Padding','replicate');
% %                 dab_image_f = imgaussfilt(dab_image_,sigma,'FilterSize',k3,'Padding','replicate');
% 
%             %     dab_image_f = imsharpen(image_dab_,'Radius',9,'Amount',1);
%             
%             
%                 % from Nir's scripts: contrast enhancing + adaptive filtering
%                 % doesn't improve ab_4g8 images
% %                 n_tiles = 8;
% %                 dab_image_ = adapthisteq(dab_image_, 'NumTiles',[n_tiles n_tiles],'ClipLimit',0.001);   % Default: NumTiles=[8 8], ClipLimit=0.01
% % 
% %                 flt_loc_pixs_size = k3; % 5
% %                 [dab_image_f, noise] = wiener2(dab_image_, [flt_loc_pixs_size flt_loc_pixs_size]);
% %                 disp([' estimated noise level before filtering: ', num2str(noise)]);
% 
%                
%                 fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
%                 subplot(211),imshow(h_image_roi),colormap(h_colormap),title('Original')
%                 subplot(212),imshow(h_image_f),colormap(h_colormap),title('Filtered')
% 
%                 fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
%                 subplot(211),imshow(dab_image_roi),colormap(dab_colormap),title('Original')
%                 subplot(212),imshow(dab_image_f),colormap(dab_colormap),title('Filtered')
%                 
%                 fname1 = strcat(fname,'_1_H_filtered.tif');
%                 saveas(fig1,fullfile(results_folder,fname1));
%                 if close_figs
%                     close(fig1);
%                 end
%                 
%                 fname2 = strcat(fname,'_2_DAB_filtered.tif');
%                 saveas(fig2,fullfile(results_folder,fname2));
%                 if close_figs
%                     close(fig2);
%                 end
% 
%                 dab_image_roi = dab_image_f;
%                 h_image_roi = h_image_f;
%             end
%             
%             %% Remove slice mask from images
%             h_image_ = h_image_roi;
%             dab_image_ = dab_image_roi;
%             
%             h_image_(~slice_mask_roi) = 255;
%             dab_image_(~slice_mask_roi) = 255;
%             
%             %%
%             if strcmp(img_type,'cfos') || strcmp(img_type,'ct695') % || strcmp(img_type,'ki67')
%                 
% %                 cell_thresh = 0.8;
%                 base_h_level = round(mean(mean(h_image_)));
%                 cell_thresh = 0.9*base_h_level/255;
% %                 [cell_thresh,~] = graythresh(h_image_roi); % automatic otsu thresh
%                 h_cell_mask = ~imbinarize(h_image_,cell_thresh);
%                 
% %                 res_cell_mask = ~imbinarize(res_image_roi,0.5);
% %                 h_cell_mask = h_cell_mask | res_cell_mask;
%                 
% %                 h_cell_mask = imbinarize(h_image_roi,'adaptive');
%                 
%                 k = 9;
%                 kernel = 1/(k*k)*ones([k,k]);
%                 h_cell_mask_2 = imfilter(h_cell_mask,kernel);
% %                 figure,imshow(h_cell_mask_2)
%                 
%                 % remove small components
%                 min_d = 3;
%                 min_con_pixels = round((pi*((min_d/2)/pixel_size)^2));
%                 h_cell_mask_f = bwareaopen(h_cell_mask_2,min_con_pixels,8);
%                 
%                 h_cell_mask = h_cell_mask_f;
% 
%                 h_image_mask_f = labeloverlay(h_image_roi,h_cell_mask_f,...
%                     'Colormap',color_map,'Transparency',transparency);
%                 
%                                 
%                 fig6 = figure('units','normalized','outerposition',[0 0 1 1]);
%                 subplot(211),imshow(h_image_roi),colormap(h_colormap)
%                 subplot(212),imshow(h_image_mask_f)
%                 
%                 fname6 = strcat(fname,'_0_H_cell_mask.tif');
%                 saveas(fig6,fullfile(results_folder,fname6));
%                 if close_figs
%                     close(fig6);
%                 end
%                 
%                 h_cell_mask_area = sum(h_cell_mask_f,[1,2]);
%                 avg_cell_r = 5;
%                 avg_cell_area = pi*avg_cell_r^2;
%                 avg_h_cell_no = round(h_cell_mask_area/avg_cell_area);
%             end
%             
%             %% Threshold
%             % white background: 255
%             % areas of interest: lower pixel values
%             
% %             color_map = [0,0,1];
% % %             color_map = [0,240,71]/255;
% %             transparency = 0.2;
% 
%             % method 1: manually chosen threshold
% %             dab_roi_mask = dab_image_ < pixel_thresh;
%             dab_roi_mask = ~imbinarize(dab_image_,pixel_thresh);
% %             figure,imshow(dab_roi_mask)
%         
%             % REMOVE NON SLICE AREA FROM THE THRESHOLD MASK
% %             dab_roi_mask(~slice_mask_roi) = 0;
% 
% 
%             if strcmp(img_type,'ct695') % || strcmp(img_type,'ki67')
%                 
%                 % remove plaques
%                 dab_roi_mask = bwpropfilt(dab_roi_mask,'EquivDiameter',[0,max_size/pixel_size]);
%                 
%                 % keep only overlap between h and dab channels
%                 dab_roi_mask_app = dab_roi_mask & h_cell_mask;
%                 
%                 fig7 = figure('units','normalized','outerposition',[0 0 1 1]);
%                 h_image_h_mask = labeloverlay(h_image_roi,h_cell_mask,...
%                     'Colormap',[1,0,0],'Transparency',transparency);
%                 h_image_h_dab_masks = labeloverlay(h_image_h_mask,dab_roi_mask,...
%                     'Colormap',[0,0,1],'Transparency',transparency);
%                 h_image_final_mask = labeloverlay(h_image_roi,dab_roi_mask_app,...
%                     'Colormap',[1,0,1],'Transparency',0);
%                 
%                 subplot(211),imshow(h_image_h_dab_masks),title('H (red), DAB (blue)')
%                 subplot(212),imshow(h_image_final_mask),title('Colocalised')
%                 
%                 fname7 = strcat(fname,'_3_0_H_cell_mask.tif');
%                 saveas(fig7,fullfile(results_folder,fname7));
%                 if close_figs
%                     close(fig7);
%                 end
%                 
%                 dab_roi_mask = dab_roi_mask_app;
%             end
%             
%             
%             dab_image_mask = labeloverlay(dab_image_roi,dab_roi_mask,...
%                 'Colormap',color_map,'Transparency',transparency);
% %             figure,imshow(dab_image_mask)
% 
%             % Mask filtering from Nir's scripts
% % 
%             % Morphological opening using a disc kernel (smoothes contours, breaks
%             % narrow bridges, eliminates small islands)
% %             strel_disk_size = 1;
% %             se = strel('disk',strel_disk_size);
% %             dab_roi_mask_2 = imopen(dab_roi_mask, se);
% %             disp([' Morphological opening. Disc size used = ', num2str(strel_disk_size)]);
% %             figure,imshow(dab_roi_mask_2)
% 
%             % Spatial smoothing of mask?
% %             dab_roi_mask_f = imfilter(dab_roi_mask,box_kernel,'replicate');
% %             dab_roi_mask_f = medfilt2(dab_roi_mask,[3,3]);
% % 
% %             figure,imshow(dab_roi_mask)
% %             figure,imshow(dab_roi_mask_f)
% 
%             % Remove small regions of connected pixels
%             % min diameter = 10um -> min area = ~309 pixels
%             % OR remove only very items smaller at this point + filter
%             % based on diameter later
%             
%             pixel_size = 0.504;
%             min_d = 2;
%             min_pix_area = round((pi*((min_d/2)/pixel_size)^2));
%             
%             min_con_pixels = min_pix_area; % 20: max ~10 um long
% %             min_con_pixels = min_pix_area;
%             connectivity = 8; % default: 4
%             dab_roi_mask_f = bwareaopen(dab_roi_mask,min_con_pixels,connectivity);
%             dab_image_mask_f = labeloverlay(dab_image_roi,dab_roi_mask_f,...
%                 'Colormap',color_map,'Transparency',transparency);
% 
%             fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
%             subplot(211),imshow(dab_image_mask),title('Threshold only')
%             subplot(212),imshow(dab_image_mask_f),title('Threshold + filtered')
%             
%             fname3 = strcat(fname,'_3_DAB_mask_1_full.tif');
%             saveas(fig3,fullfile(results_folder,fname3));
%             if close_figs
%                 close(fig3);
%             end
%             
%             % alternatively, plot mask and overlay image
% %             fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
% %             subplot(211),imshow(dab_roi_mask_f),title('Mask')
% %             subplot(212),imshow(dab_image_mask_f),title('Image + Filtered Mask')
% % 
% %             fname3 = strcat(fname,'_3_H_masks.tif');
% %             saveas(fig3,fullfile(results_folder,fname3));
% %             close(fig3);
% 
%             
% %             mask_name1 = strcat(fname,'_thresh_mask');
% %             save(fullfile(results_folder,mask_name1),'dab_roi_mask');
% %             fprintf('Mask %s saved\n',mask_name1);
%             
%             dab_roi_mask = dab_roi_mask_f;
% 
% %             mask_name2 = strcat(fname,'_thresh_filt_mask');
% %             save(fullfile(results_folder,mask_name2),'dab_roi_mask');
% %             fprintf('Mask %s saved\n',mask_name2);
%             
%             %% Separation of cells (watershed)
%             % https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
%             
%             if do_watershed
% 
%     %             bw = dab_image_;
%     %             bw(~dab_roi_mask) = 0;
%                 bw = dab_roi_mask;
% 
%                 D = -bwdist(~bw);
%     %             figure,imshow(D,[])
%                 Ld = watershed(D,26);
% 
%                 bw2 = bw;
%                 bw2(Ld == 0) = 0;
%     %             figure,imshow(bw2)
% 
%     %             m = imhmin(D,5); % can detect intensity
% 
%                 mask = imextendedmin(D,2);
%     %             figure,imshowpair(bw,mask,'blend')
% 
%                 D2 = imimposemin(D,mask);
%                 Ld2 = watershed(D2);
%                 bw3 = bw;
%                 bw3(Ld2 == 0) = 0;
% 
%                 fig5 = figure('units','normalized','outerposition',[0 0 1 1]);
%                 subplot(211),imshow(bw),title('Mask')
%                 subplot(212),imshow(bw3),title('Mask + watershed')
% 
%                 dab_roi_mask = bw3;
%                 
%                 fname5 = strcat(fname,'_3_DAB_mask_2_watershed.tif');
%                 saveas(fig5,fullfile(results_folder,fname5));
%                 if close_figs
%                     close(fig5);
%                 end
% 
%             end
%             
%             %%
% %             watershed_mask = watershed(dab_roi_mask,4);
% %             
% %             sep_dab_roi_mask = dab_roi_mask;
% %             sep_dab_roi_mask(watershed_mask == 0) = 0;
% % 
% %             figure,imshow(watershed_mask)
% %             figure,imshow(sep_dab_roi_mask)
% %             
% % %             dab_roi_mask = sep_dab_roi_mask;
% 
%             %%
%             % not v useful
% %             dab_roi_mask_2 = free_label((dab_roi_mask),6);
% %             
% %             figure,imshow(dab_roi_mask)
% %             figure,imshow(dab_roi_mask_2)
%             
%             %% Count number of cells
%             % remove line artefacts (eccentricity ~= 1)
%             max_eccentricity = 0.99;
%             
%             % remove cells smaller than x um diameter
%             min_length = (min_size/pixel_size);
%             max_length = (max_size/pixel_size);
%             
%             artefacts_idxs = [];
%             too_small_idxs = [];
%                        
%             % find cells in initial mask
%             particles_found = regionprops(dab_roi_mask,'Area','Centroid',...
%                 'BoundingBox','Circularity','Eccentricity',...
%                 'MajorAxisLength','MinorAxisLength','EquivDiameter');
%             
%             
%             % line ratio would work well but can't filter using bwpropfilt
% %             line_ratio = ([cells_found(:).MajorAxisLength]-[cells_found(:).MinorAxisLength])./[cells_found(:).MajorAxisLength]
% 
%             % remove likely artefacts (long particles) - disabled for now
% %             artefact_idxs = find([particles_found(:).Eccentricity]>max_eccentricity);
% %             artefact_idxs = [];
%             if strcmp(img_type,'ki67')
%                 max_eccentricity = 1; % 0.9
%                 artefact_idxs = find([particles_found.Eccentricity]>max_eccentricity);
%             else
%                 max_eccentricity = 1;
%                 artefact_idxs = [];
%             end
%             
% %             too_small_idxs = find([particles_found(:).MajorAxisLength]<min_length);
%             too_small_idxs = find([particles_found(:).EquivDiameter]<min_length);
%             too_large_idxs = find([particles_found(:).EquivDiameter]>max_length);
% %             too_diff = setdiff(too_small_idxs_2,too_small_idxs);
% 
%             % EquivDiameter basically based on Area
% %             [cells_found(:).Area]-([cells_found(:).EquivDiameter]./2).^2*pi
%             
%             % look at whether the area is acceptable?
% %             too_small_area = find([cells_found(:).Area]<min_pix_area/4);
% %             too_diff = setdiff(too_small_area,too_small_idxs);
%             
%             rm_idxs = [artefact_idxs, too_small_idxs, too_large_idxs];
%             
%             
%             fig4 = figure('units','normalized','outerposition',[0 0 1 1]);
%             subplot(211),hold on,title('All particles')
% %             imshow(dab_roi_mask)
%             imshow(dab_image_mask_f)
%             
%             for i = 1:length(particles_found) %rm_idxs
%                 if any(artefact_idxs==i)
%                     col = 'm';
%                     if any(too_small_idxs==i) || any(too_large_idxs==i) 
%                         col = 'r';
%                     end
%                 elseif any(too_small_idxs==i) || any(too_large_idxs==i) 
%                     col = 'r';
%                 else
%                     col = 'g';
%                 end
%                 h = rectangle('Position',particles_found(i).BoundingBox);
%                 set(h,'EdgeColor',col,'LineWidth',1);
%                 1;
%             end
%             
% %             cells_found(rm_idxs) = [];
% %             cell_no = length(cells_found);
% %             
% %             figure,imshow(dab_roi_mask),hold on
% %             for i = 1:cell_no
% %                 h = rectangle('Position',cells_found(i).BoundingBox);
% %                 set(h,'EdgeColor','r');
% %                 1;
% %             end
%             
% 
%             % create updated mask with accepted cells only
%             dab_roi_mask_f = bwpropfilt(dab_roi_mask,'Eccentricity',[0,max_eccentricity]);
%             dab_roi_mask_f = bwpropfilt(dab_roi_mask_f,'EquivDiameter',[min_length,max_length]);
%             
%             dab_image_f_mask = labeloverlay(dab_image_roi,dab_roi_mask_f,...
%                 'Colormap',color_map,'Transparency',transparency); % red
% 
% 
% %             subplot(212),hold on,title('Accepted particles')
% % %             imshow(dab_roi_plaques_mask)
% %             imshow(dab_image_f_mask)
% %             
% %             fname4 = strcat(fname,'_3_DAB_mask_3_final.tif');
% %             saveas(fig4,fullfile(results_folder,fname4));
% %             if close_figs
% %                 close(fig4);
% %             end
%             
%             %% Recreate accepted particles
%             particles_found = regionprops(dab_roi_mask_f,'Area','Centroid',...
%                 'BoundingBox','Circularity','Eccentricity',...
%                 'MajorAxisLength','EquivDiameter','Image');
%             
%             particle_no = length(particles_found);
%             
%             % normalise number of particles by area
%             particle_no = round(particle_no/tot_n_area);
%             
%             fprintf('Analysing ROI %s:\n',fname)
%             fprintf('No of particles = %d\n',particle_no)
%             
%             
%             %% Estimate number of cfos positive cells
%             if strcmp(img_type,'cfos')
%                 pos_particle_ratio = round(particle_no/avg_h_cell_no*100,1);
%                 fprintf('Cfos positive cells = %1.1f %% \n',pos_particle_ratio);
%             end
%             
%             %% Save final masks
%             mask_name1 = strcat(fname,'_mask_all');
%             save(fullfile(results_folder,mask_name1),'dab_roi_mask');
% %             fprintf('Mask %s saved\n',mask_name1);
% 
%             dab_roi_mask = dab_roi_mask_f;
% 
%             mask_name2 = strcat(fname,'_mask_accepted');
%             save(fullfile(results_folder,mask_name2),'dab_roi_mask');
% %             fprintf('Mask %s saved\n',mask_name2);
% 
% %             fprintf('Mask for %s saved\n',fname);
%             
%             %% Fiji approach
% %             setup_miji();
% %             Miji;
% %             IJ = ij.IJ;
% %             
% %             %%
% %             imp = copytoImagePlus(uint32(dab_roi_mask_f));
% %             imp.show();
% %             
% %             MIJ.run("Analyze Particles...");
% 
%             %% Density (% area covered)
% %             tot_n_area = rois_n_area(roi_idx);
%             positive_pixels = sum(dab_roi_mask,[1,2]);
%             
%             % normalise total pixels to slice area
%             total_pixels = round(size(dab_roi_mask,1)*size(dab_roi_mask,2)*tot_n_area);
% %             total_pixels = size(dab_roi_mask,1)*size(dab_roi_mask,2);
% 
%             density = round(positive_pixels/total_pixels*100,2);
% %             fprintf('Area covered for %s ROI = %1.1f %% \n', fname, density)
%             fprintf('Area covered = %1.2f %% \n', density)
% 
%             
%             %% Save final mask figure
%             title_txt = sprintf('Accepted particles. Cell count = %d, density = %1.2f',...
%                 particle_no, density);
%             subplot(212),hold on,title(title_txt)
% %             imshow(dab_roi_plaques_mask)
%             imshow(dab_image_f_mask)
%             
%             fname4 = strcat(fname,'_3_DAB_mask_3_final.tif');
%             saveas(fig4,fullfile(results_folder,fname4));
%             if close_figs
%                 close(fig4);
%             end
%             
%             %% Save results
%             results.fname = fname;
%             results.thresh_pixel = pixel_thresh;
%             
%             % info on slice mask
%             results.slice_area = slice_area;
%             results.total_area = total_area;
%             results.slice_area_norm = slice_area_norm;
% 
%             % main results
%             results.density = density;
%             results.positive_pixels = positive_pixels;
%             results.total_pixels = total_pixels;
%             results.particles_found = particles_found;
%             results.particle_no = particle_no;
%             
%             if strcmp(img_type,'cfos')
%                 results.pos_particle_ratio = pos_particle_ratio;
%             end
%             
%             fname_r = strcat(fname,'_results.mat');
%             save(fullfile(results_folder,fname_r),'results');
%             fprintf('Results %s saved\n',fname_r);
%             
%             fprintf('\n')
            
    end
%     end
end


