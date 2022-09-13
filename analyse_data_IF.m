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

    %%
    for idx = 1:length(files)
        
        file = files(idx).name;
        folder = files(idx).folder;
        fname = file(1:end-4);
    
        %% Results folder
        results_folder = fullfile(fileparts(fileparts(fileparts(folder))),'Results');

        % create subfolder with mouse name
        [~,mouse_name] = fileparts(fileparts((folder)));
        results_folder = fullfile(results_folder,mouse_name);

        if ~exist(results_folder)
            mkdir(results_folder);
        end

        %% Load reference image
        ref_file = dir(fullfile(fileparts(fileparts(fileparts(fileparts(fileparts(files(1).folder))))),'Reference_images','IF',strcat('*.tif')));
        ref_fname = fullfile(ref_file(1).folder,ref_file(1).name);
        ref_image = read_file(ref_fname,1);
        
        ref_dapi_img = ref_image(:,:,1);
        ref_iba1_img = ref_image(:,:,2);
        ref_ab_4g8_img = ref_image(:,:,3);
       
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
           
        fname_r = strcat(file(1:end-4),'_results.mat');
        save(fullfile(results_folder,fname_r),'results');
        fprintf('Results %s saved\n',fname_r);

    end
end


