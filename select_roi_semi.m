function [rois_x,rois_y] = select_roi_semi(file_)
    %% Select and save all required ROIs for a given image
    % @author: pdzialecka
    
    % The functions checks if a given ROI already exists for this image.
    % If it does, it skipps the condition to not override any existing ROIs

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
    [roi_names,roi_fnames,roi_no,roi_sizes_um] = get_roi_list();
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
    if 1 % ~all_rois_exist
        
        %% Load deconvolved images
        file_path = fullfile(folder,file);
        [h_image,dab_image,res_image] = load_deconvolved_images(file_path,1);
        
%         tiff_file = Tiff(file_path);

        %% Find DG regions for each image
        dgs_accepted = 0;
        
        dg_fname = strcat(file(1:end-11),'_dg_points.mat');

        if ~exist(dg_fname,'file')
            while ~dgs_accepted
        %         [dg_regions,dg_centroids,offset] = find_regions(h_image,file_,0);
        %         auto_find_rois = any(offset~=0,[1,2]);

                fig0 = figure('units','normalized','outerposition',[0 0 1 1]);
                imshow(h_image),colormap(h_colormap)

                title(sprintf('Select left DG point'))
                ldg_point_roi = drawpoint();

                title(sprintf('Select right DG point'))
                rdg_point_roi = drawpoint();

                ldg_point = round(ldg_point_roi.Position);
                rdg_point = round(rdg_point_roi.Position);

                %% Calculate rotation angle
                x1 = rdg_point(1);
                y1 = rdg_point(2);
                x2 = ldg_point(1);
                y2 = ldg_point(2);

                a = y2-y1;
                b = x2-x1;

                theta = -atand(a/b);

                if isnan(theta)
                    theta = 0;
                end

                %% Rotate images
                trans = [0,0];
                rot = [cosd(theta),sind(theta);...
                    -sind(theta),cosd(theta);];
                tform = rigid2d(rot,trans);

                %% Transform dg points
                [h_image_rot,ref] = imwarp(h_image,tform,'interp','cubic','FillValues',255);

                [x1_t,y1_t] = transformPointsForward(tform,x1,y1);
                [x2_t,y2_t] = transformPointsForward(tform,x2,y2);

                % account for new img size
                x1_t = x1_t-ref.XWorldLimits(1);
                y1_t = y1_t-ref.YWorldLimits(1);
                x2_t = x2_t-ref.XWorldLimits(1);
                y2_t = y2_t-ref.YWorldLimits(1);

                rdg_point_t = [x1_t,y1_t];
                ldg_point_t = [x2_t,y2_t];

                mid_point_t = find_mid_point([ldg_point_t;rdg_point_t]);
                midline_x = mid_point_t(1);

                %%
                fig00 = figure('units','normalized','outerposition',[0 0 1 1]);
                imshow(h_image_rot),colormap(h_colormap)

                drawpoint('Position',ldg_point_t);
                drawpoint('Position',rdg_point_t);
                xline(midline_x,'g')

                %% Accept or reject dg points
                answer = questdlg('Do you want to accept this ROI?',...
                    'Accept DG points?','Yes','No','Yes');

                % TODO: save dg point coords + calculate offset
                switch answer
                    case 'Yes'
                        ldg_point = ldg_point_t;
                        rdg_point = rdg_point_t;
                        dgs_accepted = 1;
                    case 'No'
                        dgs_accepted = 0;
                end

                close(fig00);
                close(fig0);
            end
        
        else
            dg_points = load(dg_fname);
            ldg_point = dg_points.ldg_point;
            rdg_point = dg_points.rdg_point;
            theta = dg_points.theta;
        end
        
        %% Rotate images
        h_image = imwarp(h_image,tform,'interp','cubic','FillValues',255);    
        dab_image = imwarp(dab_image,tform,'interp','cubic','FillValues',255);    
        res_image = imwarp(res_image,tform,'interp','cubic','FillValues',255);

        %% Create slice mask
        % TODO: into separate func, on a composite image?
        [slice_mask,slice_mask_filled] = create_slice_mask(h_image,file_);
        
        %% Calculate offset from first image
%         auto_find_rois = 0;
        
        base_dg_file = dir(fullfile(roi_folder,'*moc23*dg_points*'));
    
        if ~isempty(base_dg_file) && ~contains(file,'moc23')
            base_dg = load(fullfile(base_dg_file(1).folder,base_dg_file(1).name));
            base_rdg_point = base_dg.rdg_point;
            base_ldg_point = base_dg.ldg_point;

            auto_find_rois = 1;
            offset = [ldg_point;rdg_point] - [base_ldg_point;base_rdg_point];
        else
            auto_find_rois = 0;
            offset = [0 0; 0 0];
        end
        
        %% Save dg points
        save(fullfile(roi_folder,dg_fname),'ldg_point','rdg_point','theta','offset');
        
        %%
        for roi_idx = 1:roi_no

            %% ROI folder
            fname = file_fnames{roi_idx};
            roi_fname = file_roi_fnames{roi_idx};
    %         data_folder = fileparts(folder);
    %         roi_folder = find_roi_folder(data_folder);
            
            %% Find ROI pixel size
            roi_size_um = roi_sizes_um{roi_idx};
            roi_size = round(roi_size_um/pixel_size);
    
            %%
    %         fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'.mat');
    %         roi_fname = fullfile(roi_folder,fname);

            %%
            if ~exist(roi_fname,'file')
                %%
%                 roi_size = round(size(dab_image)/magnification);
                roi_accepted = 0;
                use_auto_roi = 1;

                while ~roi_accepted
                    %% Select ROI on H image
                    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
                    imshow(h_image),colormap(h_colormap)

                    title(sprintf('Select %s ROI',roi_names{roi_idx}))
                    

                    % roi = drawrectangle();
                    % 
                    % roi.Position = round(roi.Position);
                    % roi_x = roi.Position(2):roi.Position(2)+roi.Position(4);
                    % roi_y = roi.Position(1):roi.Position(1)+roi.Position(3);

%                     roi_point = drawpoint();
                    
                    if auto_find_rois && use_auto_roi
                        base_roi_file = dir(fullfile(roi_folder,strcat('*moc23*',roi_fnames{roi_idx},'*.mat')));
                        base_roi = load(fullfile(base_roi_file(1).folder,base_roi_file(1).name)).roi;
                        
                        if contains(fname,'L')
                            offset_ = offset(1,:);
                        elseif contains(fname,'R')
                            offset_ = offset(2,:);
                        end
                        
                        roi_x_point = (base_roi.x(1)+offset_(1))+roi_size(1)/2;
                        roi_y_point = (base_roi.y(1)+offset_(2))+roi_size(2)/2;
                        roi_point = drawpoint('Position',[roi_x_point,roi_y_point]);
                        fprintf('Estimating %s ROI location\n',fname)
                        
                    else
                        roi_point = drawpoint();
                    end

                    %%
                    roi_point.Position = round(roi_point.Position);
                    roi_x_1 = round(roi_point.Position(1)-roi_size(1)/2);
                    roi_y_1 = round(roi_point.Position(2)-roi_size(2)/2);

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

                    fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
                    ax(1) = subplot(121); imshow(h_image_roi)
                    ax(2) = subplot(122); imshow(dab_image_roi)
                    colormap(ax(1),h_colormap)
                    colormap(ax(2),dab_colormap)

                    %% Accept or reject ROI
                    answer = questdlg('Do you want to accept this ROI?',...
                        'Accept ROI?','Yes','No','Yes');

                    % Handle response
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

                        close(fig1);
                        close(fig2);
                        close(fig3);


                        % roi details
                        roi.name = roi_names{roi_idx};
                        roi.fname = roi_fnames{roi_idx};
%                         roi.magnification = magnification;
                        roi.size = roi_size;
                        roi.size_um = roi_size_um;
                        roi.x = roi_x;
                        roi.y = roi_y;
                        roi.auto_roi = auto_find_rois & use_auto_roi;

                        save(roi_fname,'roi');
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

%                         t = Tiff(img_fname,'w');
%                         tagstruct.ImageLength = size(roi_image,1);
%                         tagstruct.ImageWidth = size(roi_image,2);
%                         tagstruct.Photometric = getTag(tiff_file,'Photometric');
%                         tagstruct.Compression = getTag(tiff_file,'Compression');
%                         tagstruct.SamplesPerPixel = 2;
%                         tagstruct.PlanarConfiguration = getTag(tiff_file,'PlanarConfiguration');
%                         tagstruct.BitsPerSample = getTag(tiff_file,'BitsPerSample');
%                         
%                         t.setTag(tagstruct);
%                         write(t,roi_image);
%                         close(t);
                        
                    end
                end

            else
                fprintf('%s ROI already exists\n',fname);
            end
        end
    else
        fprintf('Skipping %s; all ROIs already exist\n',fname);
    end
end
