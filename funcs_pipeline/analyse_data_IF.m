function [] = analyse_data_IF(files,close_figs)
    %% Analyse IF data across selected ROIs: simple colocalisation
    % @author: pdzialecka
    
    %% Default input options
    if ~exist('close_figs','var')
        close_figs = 1;
    end
    
    %% G&B colormaps
    [r_colormap,g_colormap,~] = create_rgb_colormaps();
    
    %% Plot setts
    pixel_size = 0.2840; % from image info; else,1um = 3.5271 pixels
    transparency = 0.5;
    
    %% Spatial average
    k3 = 9;
    box_kernel = ones(k3,k3)*1/(k3*k3);
    
    smooth_ref = 1;

    %%
    for idx = 1:length(files)
        file = files(idx).name;
        folder = files(idx).folder;
        fname = file(1:end-4);
        
        %%
%         if contains(file,'CA1')
%             smooth_ref = 0;
%         else
%             smooth_ref = 1;
%         end
    
        %% Results folder
        results_folder = fullfile(fileparts(fileparts(fileparts(folder))),'Results');

        % create subfolder with mouse name
        [~,mouse_name] = fileparts(fileparts((folder)));
        results_folder = fullfile(results_folder,mouse_name);

        if ~exist(results_folder)
            mkdir(results_folder);
        end

        %% Load reference image
        warning ('off','all');
        roi_name = find_img_type(file);
        
        ref_file = dir(fullfile(fileparts(fileparts(fileparts(fileparts(fileparts(files(1).folder))))),'Reference_images','IF',strcat('*',roi_name,'*.tif')));
        ref_fname = fullfile(ref_file(1).folder,ref_file(1).name);
        ref_image = read_file(ref_fname,1);
        
        ref_dapi_img = ref_image(:,:,1);
        ref_iba1_img = ref_image(:,:,2);
        ref_ab_4g8_img = ref_image(:,:,3);
        
        % remove noise from reference images
        if smooth_ref
            ref_dapi_img = imfilter(ref_dapi_img,box_kernel,'replicate');
            ref_iba1_img = imfilter(ref_iba1_img,box_kernel,'replicate');
            ref_ab_4g8_img = imfilter(ref_ab_4g8_img,box_kernel,'replicate');
        end
        
        %% Load ROI image
        file_path = fullfile(folder,file);
        image = read_file(file_path);
        warning ('on','all');
        
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
        
        %% Store original images
        iba1_img_orig = iba1_img;
        ab_4g8_img_orig = ab_4g8_img;

        %% Smooth images + subtract background
        if smooth_ref
            subtract_bkg = 1;
            spatial_avg = 1;

            if spatial_avg
    %             k3 = 9;
    %             box_kernel = ones(k3,k3)*1/(k3*k3);

    %             dapi_img_f = imfilter(dapi_img,box_kernel,'replicate');
                iba1_img_f = imfilter(iba1_img,box_kernel,'replicate');
                ab_4g8_img_f = imfilter(ab_4g8_img,box_kernel,'replicate');

                if subtract_bkg
                    se2 = strel('disk',100);

                    iba1_img_f = imsubtract(iba1_img_f,imopen(iba1_img_f,se2));
                    ab_4g8_img_f = imsubtract(ab_4g8_img_f,imopen(ab_4g8_img_f,se2));
                end

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

                fname2 = strcat(fname,'_4_4g8_filtered.tif');
                saveas(fig2,fullfile(results_folder,fname2));
                if close_figs
                    close(fig2);
                end
    % 
    %             dapi_img = dapi_img_f;
                iba1_img = iba1_img_f;
                ab_4g8_img = ab_4g8_img_f;
            end
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

            fname3 = strcat(fname,'_1_iba1_normalised.tif');
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
            
            fname3 = strcat(fname,'_2_iba1_bc.tif');
            saveas(fig3,fullfile(results_folder,fname3));
            if close_figs
                close(fig3);
            end

            fname4 = strcat(fname,'_5_4g8_bc.tif');
            saveas(fig4,fullfile(results_folder,fname4));
            if close_figs
                close(fig4);
            end
            
            iba1_img = iba1_img_n;
            ab_4g8_img = ab_4g8_img_n;
            
        end
        
        %% Smooth images + subtract background
        if ~smooth_ref
            subtract_bkg = 1;
            spatial_avg = 1;

            if spatial_avg
    %             k3 = 9;
    %             box_kernel = ones(k3,k3)*1/(k3*k3);

    %             dapi_img_f = imfilter(dapi_img,box_kernel,'replicate');
                iba1_img_f = imfilter(iba1_img,box_kernel,'replicate');
                ab_4g8_img_f = imfilter(ab_4g8_img,box_kernel,'replicate');

                if subtract_bkg
                    se2 = strel('disk',100);

                    iba1_img_f = imsubtract(iba1_img_f,imopen(iba1_img_f,se2));
                    ab_4g8_img_f = imsubtract(ab_4g8_img_f,imopen(ab_4g8_img_f,se2));
                end

                fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
                subplot(211),imshow(iba1_img),colormap(r_colormap),title('Original')
                subplot(212),imshow(iba1_img_f),colormap(r_colormap),title('Filtered')

                fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
                subplot(211),imshow(ab_4g8_img),colormap(g_colormap),title('Original')
                subplot(212),imshow(ab_4g8_img_f),colormap(g_colormap),title('Filtered')

                fname1 = strcat(fname,'_2_iba1_filtered.tif');
                saveas(fig1,fullfile(results_folder,fname1));
                if close_figs
                    close(fig1);
                end

                fname2 = strcat(fname,'_5_4g8_filtered.tif');
                saveas(fig2,fullfile(results_folder,fname2));
                if close_figs
                    close(fig2);
                end
    % 
    %             dapi_img = dapi_img_f;
                iba1_img = iba1_img_f;
                ab_4g8_img = ab_4g8_img_f;
            end
        end
        
        %% Compare the original composite image with the filtered one 
        composite_img_orig = imfuse(iba1_img_orig,ab_4g8_img_orig,'ColorChannels',[1 2 0]);
        composite_img = imfuse(iba1_img,ab_4g8_img,'ColorChannels',[1 2 0]);
       
        fig7 = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(211),imshow(composite_img_orig),title('Original composite image')
        subplot(212),imshow(composite_img),title('Filtered composite image')
        
        fname7 = strcat(fname,'_7_composite_img.tif');
        saveas(fig7,fullfile(results_folder,fname7));
        if close_figs
            close(fig7);
        end

        %% Create microglia and ab masks
        if normalise_imgs
            iba1_thresh = 0.9*255;
            ab_4g8_thresh = 0.8*255;
            
            microglia_mask = iba1_img>iba1_thresh;
            ab_mask = ab_4g8_img>ab_4g8_thresh;
            
        elseif normalise_brightness
            if smooth_ref
                if contains(file,'CA1') % less sensitive to exclude WM
                    iba1_thresh = 0.6; % 0.6
                    ab_4g8_thresh = 0.5; % 0.5
                else
                    iba1_thresh = 0.5;
                    ab_4g8_thresh = 0.4;
                end
            else
                if contains(file,'CA1')
                    iba1_thresh = 0.4;
                    ab_4g8_thresh = 0.5;
                else
                    iba1_thresh = 0.4;
                    ab_4g8_thresh = 0.2;
                end
            end
            
            microglia_mask = imbinarize(iba1_img,iba1_thresh);
            ab_mask = imbinarize(ab_4g8_img,ab_4g8_thresh);
            
        else % variable threshold
            iba1_thresh = 0.4; % graythresh(iba1_img)
            microglia_mask = imbinarize(iba1_img,iba1_thresh);
            
            ab_4g8_thresh = 0.5;
            ab_mask = imbinarize(ab_4g8_img,ab_4g8_thresh);
        end
        
        
        % connect mask elements
        se = strel('disk',5);
        microglia_mask_f = imclose(microglia_mask,se);
        ab_mask_f = imclose(ab_mask,se);
        
        % remove too small and too big elements
        [~,min_size_1,max_size_1] = get_antibody_threshold('iba1');
%         min_size_1 = min_size_1/2;
        microglia_mask_f = bwpropfilt(microglia_mask_f,'EquivDiameter',[min_size_1/pixel_size,max_size_1/pixel_size]);

        [~,min_size_2,max_size_2] = get_antibody_threshold('moc23');
        min_size_2 = min_size_1; % smaller threshold to detect small cores
        ab_mask_f = bwpropfilt(ab_mask_f,'EquivDiameter',[min_size_2/pixel_size,max_size_2/pixel_size]);

        
        % Visualise the masks
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

        fname5 = strcat(fname,'_3_iba1_mask.tif');
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
        
        %% Create colocolised mask (area)
        colocalised_mask = microglia_mask & ab_mask;
        
%         composite_img = imfuse(iba1_img,ab_4g8_img,'ColorChannels',[1 2 0]);
        composite_image_mask = labeloverlay(composite_img,colocalised_mask,...
            'Colormap',[0,0,1],'Transparency',0.2);
       
        fig8 = figure('units','normalized','outerposition',[0 0 1 1]);
%         subplot(211),imshow(composite_img),title('Composite image')
%         subplot(212),
        imshow(composite_image_mask)%,title('Composite image + overlap mask')
        
        
        fname8 = strcat(fname,'_8_colocalised_area.tif');
        saveas(fig8,fullfile(results_folder,fname8));
        if close_figs
            close(fig8);
        end
        
        %% Count number of cells
        microglia_found = regionprops(microglia_mask,'Area','Centroid',...
            'BoundingBox','Circularity','Eccentricity',...
            'MajorAxisLength','MinorAxisLength','EquivDiameter');
        
        ab_found = regionprops(ab_mask,'Area','Centroid',...
            'BoundingBox','Circularity','Eccentricity',...
            'MajorAxisLength','MinorAxisLength','EquivDiameter');
        
        microglia_no = length(microglia_found);
        ab_no = length(ab_found);
        
        %% Find cells in proximity with each other
        microglia_centroids = cat(1,microglia_found.Centroid);
        ab_centroids = cat(1,ab_found.Centroid);
        
        max_dist_um = 50;
        max_dist = round(max_dist_um/pixel_size);
        
        ab_c_count = 0;
        microglia_per_ab_count = zeros(1,ab_no);
        microglia_c_count = 0;
        
        fig9 = figure;
        imshow(composite_img); hold on

        for i = 1:ab_no
            ab_centroid_i = ab_centroids(i,:);
           
            diff_i = (ab_centroid_i-microglia_centroids)/pixel_size;
            dist_i = sqrt(diff_i(:,1).^2+diff_i(:,2).^2);
           
            microglia_idxs = find(dist_i<max_dist);
%             dist_i(microglia_idxs)*pixel_size

            if ~isempty(microglia_idxs)

                h_ab = rectangle('Position',ab_found(i).BoundingBox);
                set(h_ab,'EdgeColor','g','LineWidth',1);
                
                for j = 1:length(microglia_idxs)
                    h_microglia = rectangle('Position',microglia_found(microglia_idxs(j)).BoundingBox);
                    set(h_microglia,'EdgeColor','r','LineWidth',1);
                end
                
                
                microglia_per_ab_count(i) = length(microglia_idxs);
                ab_c_count = ab_c_count+1;
                microglia_c_count = microglia_c_count+length(microglia_idxs);
            end
           
        end
        
        fname9 = strcat(fname,'_9_colocalised_count.tif');
        saveas(fig9,fullfile(results_folder,fname9));
        if close_figs
            close(fig9);
        end

        
        %% Compute % of colocalisation: count
        microglia_ab_ratio = round((microglia_c_count/microglia_no)*100,2);
        ab_microglia_ratio = round((ab_c_count/ab_no)*100,2);
        
        fprintf('*** %s file ***\n',file)
        fprintf('Colocalisation ab count = %d; microglia count = %d\n', ab_c_count,microglia_c_count);
        fprintf('Microglia no = %d. Microglia ab+ = %1.2f%%\n',microglia_no,microglia_ab_ratio);
        fprintf('Ab no = %d. Ab microglia+ = %1.2f%%\n',ab_no,ab_microglia_ratio);
                
        %% Compute % of colocalisation: density (simplest approach)
        microglia_area = sum(microglia_mask,[1,2]);
        ab_area = sum(ab_mask,[1,2]);
        colocalised_area = sum(colocalised_mask,[1,2]);
        
        microglia_ab_area_ratio = round((colocalised_area/microglia_area)*100,2);
        ab_microglia_area_ratio = round((colocalised_area/ab_area)*100,2);
        fprintf('Microglia ab+ area = %1.2f%%\n',microglia_ab_area_ratio);
        fprintf('Ab microglia+ area = %1.2f%%\n',ab_microglia_area_ratio);
        
        %% Save masks
        mask_name = strcat(fname,'_masks');
        save(fullfile(results_folder,mask_name),'microglia_mask','ab_mask','colocalised_mask');
        
        %% Save results
        results = [];
        results.orig_img_size = orig_img_size;
        results.img_size = img_size;
        results.total_area = img_size(1)*img_size(2);
        
        results.microglia_no = microglia_no;
        results.microglia_c_count = microglia_c_count;
        results.microglia_ab_ratio = microglia_ab_ratio;
        results.microglia_area = microglia_area;
        results.colocalised_area = colocalised_area;
        results.microglia_ab_area_ratio = microglia_ab_area_ratio;
        
        results.ab_no = ab_no;
        results.ab_c_count = ab_c_count;
        results.microglia_per_ab_count = microglia_per_ab_count;
        results.ab_microglia_ratio = ab_microglia_ratio;
        results.ab_area = ab_area;
        results.ab_microglia_area_ratio = ab_microglia_area_ratio;
        
        fname_r = strcat(file(1:end-4),'_results.mat');
        save(fullfile(results_folder,fname_r),'results');
        fprintf('Results %s saved\n',fname_r);

    end
end


