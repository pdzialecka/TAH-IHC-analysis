function [rois_x,rois_y] = select_roi(file_,roi_size_um)
    %% Select and save all required ROIs for a given image
    % @author: pdzialecka
    
    % The functions checks if a given ROI already exists for this image.
    % If it does, it skipps the condition to not override any existing ROIs

    %% Default input options
%     if ~exist('magnification','var')
%         magnification = 20; % 10 or 20x
%     end
    
    if ~exist('roi_size','var')
        roi_size_um = [400,400]; % in um
    end
    
    
    %%
    pixel_size = 0.504; % um
    
    % magnification approach
%     roi_1x = round([2317,4063]/pixel_size)*10; % from standard 10x size
%     roi_size = roi_1x/magnification;
    
    % arbitrary roi size in um
    roi_size = round(roi_size_um/pixel_size);
    
    
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
        [h_image,dab_image,res_image] = load_deconvolved_images(file_path);
        
%         tiff_file = Tiff(file_path);

        %%
        for roi_idx = 1:roi_no

            %% ROI folder
            fname = file_fnames{roi_idx};
            roi_fname = file_roi_fnames{roi_idx};
    %         data_folder = fileparts(folder);
    %         roi_folder = find_roi_folder(data_folder);

            %%
    %         fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'.mat');
    %         roi_fname = fullfile(roi_folder,fname);

            %%
            if ~exist(roi_fname,'file')
                %%
%                 roi_size = round(size(dab_image)/magnification);
                roi_accepted = 0;

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

                    roi_point = drawpoint();

                    %% Display ROI on DAB image
                    roi_point.Position = round(roi_point.Position);
                    roi_x_1 = round(roi_point.Position(2)-roi_size(1)/2);
                    roi_y_1 = round(roi_point.Position(1)-roi_size(2)/2);

                    roi_x = roi_x_1:roi_x_1+roi_size(1);
                    roi_y = roi_y_1:roi_y_1+roi_size(2);

                    roi_rect = drawrectangle('Position',[roi_y_1,roi_x_1,roi_size(2),roi_size(1)]);
                    title(sprintf('%s ROI',roi_names{roi_idx}));
                    
                    %% Diplay ROI on DAB image
                    fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
                    imshow(dab_image),colormap(dab_colormap)
                    roi_rect = drawrectangle('Position',[roi_y_1,roi_x_1,roi_size(2),roi_size(1)]);
                    title(sprintf('%s ROI',roi_names{roi_idx}));

                    %% Display zoomed in ROI selected
                    h_image_roi = h_image(roi_x,roi_y);
                    dab_image_roi = dab_image(roi_x,roi_y);
                    res_image_roi = res_image(roi_x,roi_y);

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

                        case 'No'
                            roi_accepted = 0;
                            close(fig1);
                            close(fig2);
                            close(fig3);
                    end

                    %% Save accepted ROI
                    if roi_accepted

                        %% ROI coordinates
                        % figures
                        fname1 = strcat(fname,'_1_H_full.tif');
                        saveas(fig1,fullfile(roi_folder,fname1));

                        fname2 = strcat(fname,'_2_DAB_full.tif');
                        saveas(fig2,fullfile(roi_folder,fname2));

%                         fname3 = strcat(fname,'_3_H_DAB_',num2str(magnification),'x.tif');
                        fname3 = strcat(fname,'_3_H_DAB_',num2str(roi_size(1)),'_x_',num2str(roi_size(2)),'.tif');
                        saveas(fig3,fullfile(roi_folder,fname3));

                        close(fig1);
                        close(fig2);
                        close(fig3);


                        % roi details
                        roi.name = roi_names{roi_idx};
                        roi.fname = roi_fnames{roi_idx};
%                         roi.magnification = magnification;
                        roi.size = roi_size;
                        roi.x = roi_x;
                        roi.y = roi_y;

                        save(roi_fname,'roi');
                        fprintf('ROI %s saved\n',fname);

                        %% ROI image as tiff
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
