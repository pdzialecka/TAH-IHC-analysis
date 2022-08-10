function [rois_x,rois_y] = select_roi_auto(file_)
    %% Automatically find and save all required ROIs for a given image
    % @author: pdzialecka
    
    % The functions checks if a given ROI already exists for this image.
    % If it does, it skipps the condition to not override any existing ROIs
    
    % Skeleton from manual function (select_roi) kept for code continuity
    % & easier switching between the two versions

    %% Default input options
%     if ~exist('magnification','var')
%         magnification = 20; % 10 or 20x
%     end
    
%     if ~exist('roi_size','var')
%         roi_size_um = [500,500]; % in um
%     end
    
    %%
    pixel_size = 0.504; % um
    
    % magnification approach
%     roi_1x = round([2317,4063]/pixel_size)*10; % from standard 10x size
%     roi_size = roi_1x/magnification;
    
    % arbitrary roi size in um
%     roi_size = round(roi_size_um/pixel_size);
    
    %% H DAB colormaps
    % for easier visualisation
    [h_colormap,dab_colormap] = create_hdab_colormaps();
    
    %% ROIs
    [roi_names,roi_fnames,roi_no] = get_roi_list();
    roi_order_no = 1:roi_no;
    
    rois_x = {};
    rois_y = {};

    %% Directory info
    file = file_.name;
    folder = file_.folder;
    
    %% ROI folder
%     data_folder = fileparts(folder);
    [roi_folder,mouse_name] = find_roi_folder(folder);
    
    file_fnames = {};
    file_roi_fnames = {};
    file_exists = [];
    
    for roi_idx = 1:roi_no
        fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'.mat');
        file_roi_fnames{roi_idx} = fullfile(roi_folder,fname);
        file_fnames{roi_idx} = fname(1:end-4);
        file_exists(roi_idx) = exist(fullfile(roi_folder,fname));
    end
    
    all_rois_exist = all(file_exists);
    
    %% ROI images folder
    roi_img_folder = fullfile(fileparts(fileparts(folder)),'ROI_images',mouse_name);

    if ~exist(roi_img_folder)
        mkdir(roi_img_folder);
    end
    
    %%
    if ~all_rois_exist
        
        %% Load deconvolved images
        file_path = fullfile(folder,file);
        [h_image,dab_image,res_image] = load_deconvolved_images(file_path,1);
        
%         tiff_file = Tiff(file_path);

        %% Find DG regions for each image
        [dg_regions,dg_centroids,offset,theta,rois] = find_regions(h_image,file_,1);
%         auto_find_rois = any(offset~=0,[1,2]);
        auto_find_rois = 1;

        %% Rotate images
        trans = [0,0];
        rot = [cosd(theta),sind(theta);...
            -sind(theta),cosd(theta);];
        tform = rigid2d(rot,trans);
        
        h_image = imwarp(h_image,tform,'interp','cubic','FillValues',255);    
        dab_image = imwarp(dab_image,tform,'interp','cubic','FillValues',255);    
        res_image = imwarp(res_image,tform,'interp','cubic','FillValues',255);

        %% Extract ROI info
%         dg_rois = rois{1};
%         ca1_rois = rois{2};
%         ca3_rois = rois{3};
%         cortex_rois = rois{4};
        
        %% Create slice mask
        [slice_mask,slice_mask_filled] = create_slice_mask(h_image,file_);

        %%
        for roi_idx = 1:roi_no

            %% ROI folder
            fname = file_fnames{roi_idx};
            file_roi_fname = file_roi_fnames{roi_idx};
            roi_fname = roi_fnames{roi_idx};
            
            %% ROI info
            [coords,this_roi,coords_field] = extract_roi_coords(rois,roi_fname);
%             if contains(roi_fname,'R')
%                 coords_field = 'R_coords';
%             elseif contains(roi_fname,'L')
%                 coords_field = 'L_coords';
%             end
%             
%             if contains(roi_fname,'DG')
%                 this_roi = dg_rois;
%             elseif contains(roi_fname,'CA1')
%                 this_roi = ca1_rois;
%             elseif contains(roi_fname,'CA3')
%                 this_roi = ca3_rois;
%             elseif contains(roi_fname,'cortex')
%                 this_roi = cortex_rois;
%             end
%             
%             coords = getfield(this_roi,coords_field);

            %%
            if ~exist(file_roi_fname,'file')
                %%
                roi_accepted = 0;
                use_auto_roi = 1;

                while ~roi_accepted
                    
                    %% Select ROI on H image
                    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
                    imshow(h_image),colormap(h_colormap)

                    title(sprintf('Select %s ROI',roi_names{roi_idx}))
                    
%                     if auto_find_rois && use_auto_roi
%                         base_roi_file = dir(fullfile(roi_folder,strcat('*moc23*',roi_fnames{roi_idx},'*.mat')));
%                         base_roi = load(fullfile(base_roi_file(1).folder,base_roi_file(1).name)).roi;
%                         
%                         if contains(fname,'L')
%                             offset_ = offset(1,:);
%                         elseif contains(fname,'R')
%                             offset_ = offset(2,:);
%                         end
%                         
%                         roi_x_point = (base_roi.x(1)+offset_(1))+roi_size(1)/2;
%                         roi_y_point = (base_roi.y(1)+offset_(2))+roi_size(2)/2;
%                         roi_point = drawpoint('Position',[roi_x_point,roi_y_point]);
%                         fprintf('Estimating %s ROI location\n',fname)
%                         
%                     else
%                         roi_point = drawpoint();
%                     end

                    %%
%                     roi_point.Position = round(roi_point.Position);
%                     roi_x_1 = round(roi_point.Position(1)-roi_size(1)/2);
%                     roi_y_1 = round(roi_point.Position(2)-roi_size(2)/2);
                    
                    roi_x_1 = coords(1);
                    roi_y_1 = coords(2);
                    roi_size = coords(3:4);
                    roi_size_um = this_roi.dims_um;

                    roi_x = roi_x_1:roi_x_1+roi_size(1);
                    roi_y = roi_y_1:roi_y_1+roi_size(2);

                    roi_rect = drawrectangle('Position',[roi_x_1,roi_y_1,roi_size(1),roi_size(2)]);
                    title(sprintf('%s ROI',roi_names{roi_idx}));
                    
                    %% Diplay ROI on DAB image
                    fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
                    imshow(dab_image),colormap(dab_colormap)
                    roi_rect = drawrectangle('Position',[roi_x_1,roi_y_1,roi_size(1),roi_size(2)]);
                    title(sprintf('%s ROI',roi_names{roi_idx}));

                    %% Display zoomed in ROI selected
                    h_image_roi = h_image(roi_y,roi_x);
                    dab_image_roi = dab_image(roi_y,roi_x);
                    res_image_roi = res_image(roi_y,roi_x);
                    slice_mask_roi = slice_mask(roi_y,roi_x);

                    fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
                    ax(1) = subplot(121); imshow(h_image_roi)
                    ax(2) = subplot(122); imshow(dab_image_roi)
                    colormap(ax(1),h_colormap)
                    colormap(ax(2),dab_colormap)
                    
                    %% Display zoomed in ROI w slice mask
                    h_image_roi_smask = labeloverlay(h_image_roi,~slice_mask_roi,...
                        'Colormap',[0,0,1],'Transparency',0.7);
                    dab_image_roi_smask = labeloverlay(dab_image_roi,~slice_mask_roi,...
                        'Colormap',[0,0,1],'Transparency',0.7);
                    
                    fig4 = figure('units','normalized','outerposition',[0 0 1 1]);
                    ax(1) = subplot(121); imshow(h_image_roi_smask)
                    ax(2) = subplot(122); imshow(dab_image_roi_smask)
                    colormap(ax(1),h_colormap)
                    colormap(ax(2),dab_colormap)
                    
                    slice_area = sum(slice_mask_roi,[1,2]);
                    total_area = size(h_image_roi,1)*size(h_image_roi,2);
                    slice_area_norm = slice_area/total_area; % as % of total area

                    %% Accept or reject ROI
%                     answer = questdlg('Do you want to accept this ROI?',...
%                         'Accept ROI?','Yes','No','Yes');

                    % accept auto detected ROI
                    answer = 'Yes';

                    switch answer
                        case 'Yes'
                            roi_accepted = 1;
                            rois_x{roi_idx} = roi_x;
                            rois_y{roi_idx} = roi_y;
                            
                            if auto_find_rois && use_auto_roi
                                fprintf('Auto ROI accepted\n')
                            end
                            
                        case 'No'
                            roi_accepted = 0;
                            use_auto_roi = 0;
                            if auto_find_rois && use_auto_roi
                                fprintf('Auto ROI rejected; selecting ROI manually\n')
                            end
                            close(fig1);
                            close(fig2);
                            close(fig3);
                    end

                    %% Save accepted ROI
                    if roi_accepted

                        %% Save ROI coordinates
                        % figures
                        fname1 = strcat(fname,'_1_H_full.tif');
                        saveas(fig1,fullfile(roi_folder,fname1));

                        fname2 = strcat(fname,'_2_DAB_full.tif');
                        saveas(fig2,fullfile(roi_folder,fname2));

%                         fname3 = strcat(fname,'_3_H_DAB_',num2str(magnification),'x.tif');
                        fname3 = strcat(fname,'_3_H_DAB_',num2str(roi_size_um(1)),'x',num2str(roi_size_um(2)),'um.tif');
                        saveas(fig3,fullfile(roi_folder,fname3));
                        
                        fname4 = strcat(fname,'_3_H_DAB_',num2str(roi_size_um(1)),'x',num2str(roi_size_um(2)),'um_slice_mask.tif');
                        saveas(fig4,fullfile(roi_folder,fname4));


                        close(fig1);
                        close(fig2);
                        close(fig3);
                        close(fig4);


                        % roi details
                        roi = [];
                        roi.name = roi_names{roi_idx};
                        roi.fname = roi_fnames{roi_idx};
%                         roi.magnification = magnification;
                        roi.size = roi_size;
                        roi.size_um = roi_size_um;
                        roi.x = roi_x;
                        roi.y = roi_y;
                        roi.auto_roi = auto_find_rois & use_auto_roi;
                        
                        roi.slice_area = slice_area;
                        roi.total_area = total_area;
                        roi.slice_area_norm = slice_area_norm;

                        save(file_roi_fname,'roi');
                        fprintf('ROI %s saved\n',fname);

                        %% Save ROI image as tiff
                        roi_image = cat(3,h_image_roi,dab_image_roi,res_image_roi);
                        img_fname = fullfile(roi_img_folder,strcat(fname,'.tif'));
%                         save(img_fname,'roi_image');

                        for i = 1:size(roi_image,3)
                            if i==1
                                imwrite(roi_image(:,:,i),img_fname,'Compression','none','WriteMode','overwrite');
                            else
                                imwrite(roi_image(:,:,i),img_fname,'Compression','none','WriteMode','append');
                            end
                        end
                        
                    end
                end

            else
                fprintf('%s ROI already exists\n',fname);
            end
        end
    else
        fprintf('Skipping %s; all ROIs already exist\n',file);
    end
end
