function [rois_x,rois_y] = select_roi(file_,magnification)
    %%

    %% Default input options
    if ~exist('magnification','var')
        magnification = 10; % 10 or 20x
    end
    
    %% H DAB colormaps
    % for easier visualisation
    dab_col = (1-[0.26814753,0.57031375,0.77642715]);
    h_col = (1-[0.6500286,0.704031,0.2860126]);

    dab_colormap = [(linspace(dab_col(1),1,256)'),(linspace(dab_col(2),1,256)'),(linspace(dab_col(3),1,256)')];
    h_colormap = [(linspace(h_col(1),1,256)'),(linspace(h_col(2),1,256)'),(linspace(h_col(3),1,256)')];
    
    %% ROIs
    roi_names = {'Right Hipp','Left Hipp','Right Cortex','Left Cortex','Control Area'};
    roi_no = length(roi_names);
    roi_fnames = {'hipp_R','hipp_L','cortex_R','cortex_L','control'};
    roi_order_no = 1:roi_no;
    
    rois_x = {};
    rois_y = {};
    
    %% File to be processed
    file = file_.name;
    folder = file_.folder;
    
    %%
    for roi_idx = 1:roi_no

        %% ROI folder
        roi_folder = fullfile(fileparts(fileparts(folder)),'ROIs');

        % create subfolder with mouse name
        mouse_name = file(1:9);
        roi_folder = fullfile(roi_folder,mouse_name);

        if ~exist(roi_folder)
            mkdir(roi_folder);
        end
        
        %%
        fname = strcat(file(1:end-11),'_',num2str(roi_order_no(roi_idx)),'_roi_',roi_fnames{roi_idx},'.mat');
        roi_fname = fullfile(roi_folder,fname);
        
        %%
        if ~exist(roi_fname,'file')

            %% Load deconvolved images
            image = read_file(fullfile(folder,file));
            h_image = image(:,:,1);
            dab_image = image(:,:,2);
        %     image_res = image(:,:,3);

            %% Rotate images
            h_image = imrotate(h_image,-90);
            dab_image = imrotate(dab_image,-90);

            %%
            roi_size = round(size(dab_image)/magnification);
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

                %% Diplay ROI on DAB image
                roi_point.Position = round(roi_point.Position);
                roi_x_1 = round(roi_point.Position(2)-roi_size(1)/2);
                roi_y_1 = round(roi_point.Position(1)-roi_size(2)/2);

                roi_x = roi_x_1:roi_x_1+roi_size(1);
                roi_y = roi_y_1:roi_y_1+roi_size(2);

                roi_rect = drawrectangle('Position',[roi_y_1,roi_x_1,roi_size(2),roi_size(1)]);
                title(sprintf('%s ROI',roi_names{roi_idx}));

                fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
                imshow(dab_image),colormap(dab_colormap)
                roi_rect = drawrectangle('Position',[roi_y_1,roi_x_1,roi_size(2),roi_size(1)]);
                title(sprintf('%s ROI',roi_names{roi_idx}));


                %% Display ROI selected
                h_image_roi = h_image(roi_x,roi_y);
                dab_image_roi = dab_image(roi_x,roi_y);

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

                %% Save accepted ROI coordinates
                if roi_accepted
            %                     fname = strcat(file(1:end-11),'_roi_',roi_fnames{roi_idx});

                    % figures
                    fname1 = strcat(fname,'_1_H_full.tif');
                    saveas(fig1,fullfile(roi_folder,fname1));

                    fname2 = strcat(fname,'_2_DAB_full.tif');
                    saveas(fig2,fullfile(roi_folder,fname2));

                    fname3 = strcat(fname,'_3_H_DAB_',num2str(magnification),'x.tif');
                    saveas(fig3,fullfile(roi_folder,fname3));

                    close(fig1);
                    close(fig2);
                    close(fig3);


                    % roi details
                    roi.name = roi_names{roi_idx};
                    roi.fname = roi_fnames{roi_idx};
                    roi.x = roi_x;
                    roi.y = roi_y;

%                     fname = strcat(file(1:end-11),'_roi_',roi_fnames{roi_idx},'.mat');
%                     save(fullfile(roi_folder,fname),'roi');
                    save(roi_fname,'roi');
                    fprintf('ROI %s saved\n',fname);
                    
                end
            end
            
        else
            fprintf('%s ROI already exists\n',fname);
        end
    end
end
