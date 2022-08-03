function [] = analyse_data(files,load_rois,close_figs)
    %% Analyse IHC data across selected ROIs
    % @author: pdzialecka
    
    % Variables to explore:
    % * pixel_thresh_factor (default: 1)
    % * min_con_pixels (default: 25 at 20x, 100 at 10x)
    % * k3 (spatial filter variable)
    % * roi_size / magnification (may affect other parameters)
    
    %% Default input options
    if ~exist('load_rois','var')
        load_rois = 1;
    end
    
    if ~exist('close_figs','var')
        close_figs = 1;
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
% 
%         if ~exist(roi_img_folder)
%             mkdir(roi_img_folder);
%         end

        %% Image type: threshold + min particle size
%         idxs_ = strfind(file,'_');
%         image_type = file(idxs_(1)+1:idxs_(2)-1);
        img_type = find_img_type(file);
        
%         if contains(file,'moc23')
%             image_type = 'moc23';
% 
%         elseif contains(file,'cfos')
%             image_type = 'cfos';
%             do_watershed = 1;
% 
%         elseif contains(file,'GFAP')
%             image_type = 'GFAP';
% 
%         elseif contains(file,'Iba1')
%             image_type = 'Iba1';
%         end
        
        [pixel_thresh,min_size,do_watershed] = get_antibody_threshold(img_type);
        
    
        %% Results folder
        results_folder = fullfile((fileparts(fileparts(folder))),'Results');

        % create subfolder with mouse name
        [~,mouse_name] = fileparts((folder));
        results_folder = fullfile(results_folder,mouse_name);

        if ~exist(results_folder)
            mkdir(results_folder);
        end

        %% Load deconvolved images
%         file_path = fullfile(folder,file);
%         [h_image,dab_image] = load_deconvolved_images(file_path);
        
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
            
            %% Extract ROIs from main image (slower)
%             dab_image_roi = dab_image(rois_x{roi_idx},rois_y{roi_idx});
%             h_image_roi = h_image(rois_x{roi_idx},rois_y{roi_idx});
% 
            %% Load ROI image
            file_path = fullfile(roi_img_folder,strcat(fname,'.tif'));
            [h_image_roi,dab_image_roi,~] = load_deconvolved_images(file_path);
            
            %% Denoise images
            % This filter gets rid of most of the random pepper noise
            
            % Gaussian and box filters comparable:
            % box blur applied many times = approximation of Gaussian blur
            % gaussian f requires sigma estimation and potentially larger filter size
            
            dab_image_ = dab_image_roi;
            h_image_ = h_image_roi;

            % Spatial smoothing
            spatial_avg = 1;
            if spatial_avg
                k3 = 3;

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
                if close_figs
                    close(fig1);
                end
                
                fname2 = strcat(fname,'_2_DAB_filtered.tif');
                saveas(fig2,fullfile(results_folder,fname2));
                if close_figs
                    close(fig2);
                end

                dab_image_ = dab_image_f;
                h_image_ = h_image_f;
            end

            %% Threshold
            % white background: 255
            % areas of interest: lower pixel values
            
            color_map = [0,0,1];
%             color_map = [0,240,71]/255;
            transparency = 0.2;

            % method 1: manually chosen threshold
%             dab_roi_mask = dab_image_ < pixel_thresh;
            dab_roi_mask = ~imbinarize(dab_image_,pixel_thresh);
            
            % method 2: automated Otsu's threshold, can use to find approx value
%             pixel_thresh = graythresh(dab_image_);
%             pixel_thresh_factor = 1;
%             pixel_thresh = pixel_thresh*pixel_thresh_factor
%             dab_roi_mask = ~imbinarize(dab_image_,pixel_thresh);
            

%             figure,imshow(dab_image_),colormap(dab_colormap)
            
            dab_image_mask = labeloverlay(dab_image_,dab_roi_mask,...
                'Colormap',color_map,'Transparency',transparency);
%             figure,imshow(dab_image_mask)

            % Mask filtering from Nir's scripts
% 
            % Morphological opening using a disc kernel (smoothes contours, breaks
            % narrow bridges, eliminates small islands)
%             strel_disk_size = 1;
%             se = strel('disk',strel_disk_size);
%             dab_roi_mask_2 = imopen(dab_roi_mask, se);
%             disp([' Morphological opening. Disc size used = ', num2str(strel_disk_size)]);
%             figure,imshow(dab_roi_mask_2)

            % Spatial smoothing of mask?
%             dab_roi_mask_f = imfilter(dab_roi_mask,box_kernel,'replicate');
%             dab_roi_mask_f = medfilt2(dab_roi_mask,[3,3]);
% 
%             figure,imshow(dab_roi_mask)
%             figure,imshow(dab_roi_mask_f)

            % Remove small regions of connected pixels
            % min diameter = 10um -> min area = ~309 pixels
            % OR remove only very items smaller at this point + filter
            % based on diameter later
            
            pixel_size = 0.504;
            min_d = 2;
            min_pix_area = round((pi*((min_d/2)/pixel_size)^2));
            
            min_con_pixels = min_pix_area; % 20: max ~10 um long
%             min_con_pixels = min_pix_area;
            connectivity = 8; % default: 4
            dab_roi_mask_f = bwareaopen(dab_roi_mask,min_con_pixels,connectivity);
            dab_image_mask_f = labeloverlay(dab_image_,dab_roi_mask_f,...
                'Colormap',color_map,'Transparency',transparency);

            fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(121),imshow(dab_image_mask),title('Threshold only')
            subplot(122),imshow(dab_image_mask_f),title('Threshold + filtered')
            
            fname3 = strcat(fname,'_3_DAB_mask_1_full.tif');
            saveas(fig3,fullfile(results_folder,fname3));
            if close_figs
                close(fig3);
            end
            
            % alternatively, plot mask and overlay image
%             fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
%             subplot(121),imshow(dab_roi_mask_f),title('Mask')
%             subplot(122),imshow(dab_image_mask_f),title('Image + Filtered Mask')
% 
%             fname3 = strcat(fname,'_3_H_masks.tif');
%             saveas(fig3,fullfile(results_folder,fname3));
%             close(fig3);

            
%             mask_name1 = strcat(fname,'_thresh_mask');
%             save(fullfile(results_folder,mask_name1),'dab_roi_mask');
%             fprintf('Mask %s saved\n',mask_name1);
            
            dab_roi_mask = dab_roi_mask_f;

%             mask_name2 = strcat(fname,'_thresh_filt_mask');
%             save(fullfile(results_folder,mask_name2),'dab_roi_mask');
%             fprintf('Mask %s saved\n',mask_name2);
            
            %% Separation of cells (watershed)
            % https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
            
            if do_watershed

    %             bw = dab_image_;
    %             bw(~dab_roi_mask) = 0;
                bw = dab_roi_mask;

                D = -bwdist(~bw);
    %             figure,imshow(D,[])
                Ld = watershed(D,26);

                bw2 = bw;
                bw2(Ld == 0) = 0;
    %             figure,imshow(bw2)

    %             m = imhmin(D,5); % can detect intensity

                mask = imextendedmin(D,2);
    %             figure,imshowpair(bw,mask,'blend')

                D2 = imimposemin(D,mask);
                Ld2 = watershed(D2);
                bw3 = bw;
                bw3(Ld2 == 0) = 0;

                fig5 = figure('units','normalized','outerposition',[0 0 1 1]);
                subplot(121),imshow(bw),title('Mask')
                subplot(122),imshow(bw3),title('Mask + watershed')

                dab_roi_mask = bw3;
                
                fname5 = strcat(fname,'_3_DAB_mask_2_watershed.tif');
                saveas(fig5,fullfile(results_folder,fname5));
                if close_figs
                    close(fig5);
                end

            end
            
            %%
%             watershed_mask = watershed(dab_roi_mask,4);
%             
%             sep_dab_roi_mask = dab_roi_mask;
%             sep_dab_roi_mask(watershed_mask == 0) = 0;
% 
%             figure,imshow(watershed_mask)
%             figure,imshow(sep_dab_roi_mask)
%             
% %             dab_roi_mask = sep_dab_roi_mask;

            %%
            % not v useful
%             dab_roi_mask_2 = free_label((dab_roi_mask),6);
%             
%             figure,imshow(dab_roi_mask)
%             figure,imshow(dab_roi_mask_2)
            
            %% Count number of cells
            % remove line artefacts (eccentricity ~= 1)
            max_eccentricity = 0.99;
            
            % remove cells smaller than x um diameter
            min_length = (min_size/pixel_size);
            
            artefacts_idxs = [];
            too_small_idxs = [];
                       
            % find cells in initial mask
            particles_found = regionprops(dab_roi_mask,'Area','Centroid',...
                'BoundingBox','Circularity','Eccentricity',...
                'MajorAxisLength','MinorAxisLength','EquivDiameter');
            
            
            % line ratio would work well but can't filter using bwpropfilt
%             line_ratio = ([cells_found(:).MajorAxisLength]-[cells_found(:).MinorAxisLength])./[cells_found(:).MajorAxisLength]
            artefact_idxs = find([particles_found(:).Eccentricity]>max_eccentricity);
            
%             too_small_idxs = find([particles_found(:).MajorAxisLength]<min_length);
            too_small_idxs = find([particles_found(:).EquivDiameter]<min_length);
%             too_diff = setdiff(too_small_idxs_2,too_small_idxs);

            % EquivDiameter basically based on Area
%             [cells_found(:).Area]-([cells_found(:).EquivDiameter]./2).^2*pi
            
            % look at whether the area is acceptable?
%             too_small_area = find([cells_found(:).Area]<min_pix_area/4);
%             too_diff = setdiff(too_small_area,too_small_idxs);
            
            rm_idxs = [artefact_idxs, too_small_idxs];
            
            
            fig4 = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(121),hold on,title('All particles')
%             imshow(dab_roi_mask)
            imshow(dab_image_mask_f)
            
            for i = 1:length(particles_found) %rm_idxs
                if any(artefact_idxs==i)
                    col = 'm';
                    if any(too_small_idxs==i)
                        col = 'r';
                    end
                elseif any(too_small_idxs==i)
                    col = 'r';
                else
                    col = 'g';
                end
                h = rectangle('Position',particles_found(i).BoundingBox);
                set(h,'EdgeColor',col,'LineWidth',1);
                1;
            end
            
%             cells_found(rm_idxs) = [];
%             cell_no = length(cells_found);
%             
%             figure,imshow(dab_roi_mask),hold on
%             for i = 1:cell_no
%                 h = rectangle('Position',cells_found(i).BoundingBox);
%                 set(h,'EdgeColor','r');
%                 1;
%             end
            

            % create updated mask with plaques only
            dab_roi_plaques_mask = bwpropfilt(dab_roi_mask,'Eccentricity',[0,max_eccentricity]);
%             dab_roi_plaques_mask = bwpropfilt(dab_roi_plaques_mask,'MajorAxisLength',[min_length,inf]);
            dab_roi_plaques_mask = bwpropfilt(dab_roi_plaques_mask,'EquivDiameter',[min_length,inf]);
            
            dab_image_plaque_mask = labeloverlay(dab_image_,dab_roi_plaques_mask,...
                'Colormap',color_map,'Transparency',transparency); % red

%             figure,
            subplot(122),hold on,title('Accepted particles')
%             imshow(dab_roi_plaques_mask)
            imshow(dab_image_plaque_mask)
            
            fname4 = strcat(fname,'_3_DAB_mask_3_final.tif');
            saveas(fig4,fullfile(results_folder,fname4));
            if close_figs
                close(fig4);
            end
            
            %% Recreate mask with accepted particles only
            particles_found = regionprops(dab_roi_plaques_mask,'Area','Centroid',...
                'BoundingBox','Circularity','Eccentricity',...
                'MajorAxisLength','EquivDiameter','Image');
            
            particle_no = length(particles_found);
            fprintf('Analysing ROI %s:\n',fname)
            fprintf('No of particles = %d\n',particle_no)
            
            %% Save final masks
            mask_name1 = strcat(fname,'_mask_all');
            save(fullfile(results_folder,mask_name1),'dab_roi_mask');
%             fprintf('Mask %s saved\n',mask_name1);

            dab_roi_mask = dab_roi_plaques_mask;

            mask_name2 = strcat(fname,'_mask_accepted');
            save(fullfile(results_folder,mask_name2),'dab_roi_mask');
%             fprintf('Mask %s saved\n',mask_name2);

%             fprintf('Mask for %s saved\n',fname);
            
            %%
%             dab_roi_mask = dab_roi_plaques_mask;
            
            %% Fiji approach
%             setup_miji();
%             Miji;
%             IJ = ij.IJ;
%             
%             %%
%             imp = copytoImagePlus(uint32(dab_roi_mask_f));
%             imp.show();
%             
%             MIJ.run("Analyze Particles...");
            
            %% Save ROI images
%             roi_image = cat(3,h_image_,dab_image_);
% %             magnification = 10; % load from inside roi_file later
% %             img_fname = fullfile(roi_img_folder,strcat(fname,'_',num2str(magnification),'x'));
%             
%             img_fname = fullfile(roi_img_folder,strcat(fname));
%             save(img_fname,'roi_image');

            %% Density (% area covered)
            positive_pixels = sum(dab_roi_mask,[1,2]);
            total_pixels = size(dab_roi_mask,1)*size(dab_roi_mask,2);

            density = round(positive_pixels/total_pixels*100,1);
%             fprintf('Area covered for %s ROI = %1.1f %% \n', fname, density)
            fprintf('Area covered = %1.1f %% \n', density)

            
            %% Save results
            results.fname = fname;
            results.thresh_pixel = pixel_thresh;
            results.density = density;
            results.positive_pixels = positive_pixels;
            results.total_pixels = total_pixels;
            results.particles_found = particles_found;
            results.particle_no = particle_no;
            
            fname_r = strcat(fname,'_results.mat');
            save(fullfile(results_folder,fname_r),'results');
            fprintf('Results %s saved\n',fname_r);
            
            fprintf('\n')
            
        end
    end
end


