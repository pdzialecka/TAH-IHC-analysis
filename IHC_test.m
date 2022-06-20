%%
clc;
close all;
clear all;

%% Analysis folder
analysis_folder = 'C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC';
addpath(genpath(analysis_folder));

%% Find data folder within cohort folder
% here, use example image inside Data subfolder of analysis folder
base_folder = fullfile(analysis_folder);
cohort_folder = 'Cohort_3';  % ****
data_folder = fullfile(base_folder,'Data',cohort_folder,'Images');

%% Find all files inside the data folder
files = dir(fullfile(data_folder,'**','*.svs'));

%% Add FIJI / MIJ paths
setup_miji();

%% Start fiji
Miji;
IJ = ij.IJ;

%% File info
idx = 1; % select file idx to be processed
file = files(idx).name;
folder = files(idx).folder;
fname = fullfile(folder,file);

%% Create folder for processed images
processed_folder = fullfile(fileparts(fileparts(folder)),'Images_processed');

% create subfolder with mouse name
[~,mouse_name] = fileparts(folder);
processed_folder = fullfile(processed_folder,mouse_name);

if ~exist(processed_folder)
    mkdir(processed_folder);
end

%% Create ROI folder
roi_folder = fullfile(fileparts(fileparts(folder)),'ROIs');
roi_folder = fullfile(roi_folder,mouse_name);

if ~exist(roi_folder)
    mkdir(roi_folder);
end

%% Create results folder
results_folder = fullfile(fileparts(fileparts((folder))),'Results');
results_folder = fullfile(results_folder,mouse_name);

if ~exist(results_folder)
    mkdir(results_folder);
end

%% Load RGB svs file
image_RGB = imread(fname);

%% Open file in Fiji
imp = copytoImagePlus(image_RGB);
% image_RGB_sml = image_RGB(end/4:2*end/4,end/4:2*end/4,:);
% imp = copytoImagePlus(image_RGB_sml);
imp.show();

%% Color deconvolution in Fiji
MIJ.run("RGB Color");
imp.close();
MIJ.run("Colour Deconvolution", "vectors=[H DAB]");

%% Save results
% Save combined tif file
MIJ.selectWindow("new (RGB)")
fname_1 = fullfile(processed_folder,strcat(file(1:end-4),'.tif'));
IJ.save(fname_1)
MIJ.run("Close")

% this ensures file saved as 8-bit
MIJ.selectWindow("Colour Deconvolution")
MIJ.run("Close")

% Save deconvolved images as a stack
MIJ.run("Images to Stack", "use");
MIJ.run("Grays");
fname_2 = fullfile(processed_folder,strcat(file(1:end-4),'_deconv.tif'));
IJ.save(fname_2)

%% Close Fiji
MIJ.closeAllWindows;
MIJ.exit;

%% Load deconvolved tiff file
% [h_image,dab_image] = load_deconvolved_images(fname_2);

image = read_file(fname_2);
h_image = image(:,:,1);
dab_image = image(:,:,2);

% rotate the images
h_image = imrotate(h_image,-90);
dab_image = imrotate(dab_image,-90);


%% Define ROIs
% define ROI names
roi_names = {'Right_Hipp','Left_Hipp','Right_Cortex','Left_Cortex','Control_Area'};
% roi_names = {'Right_Hipp','Left_Hipp'};
roi_no = length(roi_names);

% select magnification of choice for ROI images, e.g. 10 or 20x
magnification = 10; % e.g. 10/20x

% compute the size of ROI
roi_size = round(size(dab_image)/magnification);

% preparing roi x y structures
rois_x = {};
rois_y = {};

% colormaps for easier visualisation
[h_colormap,dab_colormap] = create_hdab_colormaps();

%% Looping through ROIs
for roi_idx = 1:roi_no
    % roi_idx = 1;
    
    roi_accepted = 0;

    while ~roi_accepted
        %% Select ROI on H image
        fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
        imshow(h_image),colormap(h_colormap)

        title(sprintf('Select %s ROI',roi_names{roi_idx}))
        roi_point = drawpoint();

        %% Create ROI and display on H image
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
            
            roi_fname = strcat(file(1:end-4),'_roi_',roi_names{roi_idx});

            % figures
            fname1 = strcat(roi_fname,'_1_H_full.tif');
            saveas(fig1,fullfile(roi_folder,fname1));

            fname2 = strcat(roi_fname,'_2_DAB_full.tif');
            saveas(fig2,fullfile(roi_folder,fname2));

            fname3 = strcat(roi_fname,'_3_H_DAB_',num2str(magnification),'x.tif');
            saveas(fig3,fullfile(roi_folder,fname3));

            close(fig1);
            close(fig2);
            close(fig3);


            % roi details
            roi_fname = strcat(file(1:end-4),'_roi_',roi_names{roi_idx},'.mat');
            roi.name = roi_names{roi_idx};
    %         roi.name = roi_names{roi_idx};
            roi.x = roi_x;
            roi.y = roi_y;

            save(fullfile(roi_folder,roi_fname),'roi');
            fprintf('ROI %s saved\n',roi_fname);

        end
    end
end


%% Identify antibody and assign threshold
if contains(fname,'moc23')
    image_type = 'moc23';
    pixel_thresh = 0.65; % 160;

elseif contains(fname,'cfos')
    image_type = 'cfos';
    pixel_thresh = 0.2; % decrease to detect only v dark cells

elseif contains(fname,'GFAP')
    image_type = 'GFAP';
    pixel_thresh = 0.5;

elseif contains(fname,'Iba1')
    image_type = 'Iba1';
    pixel_thresh = 0.35;
end

%% Calculate % of area covered by the antibody

% spatial averaging
k3 = 3;
box_kernel = ones(k3,k3)*1/(k3*k3);


for roi_idx = 1:roi_no
    %% Load ROI information
    roi_fname = fullfile(roi_folder,strcat(file(1:end-4),'_roi_',roi_names{roi_idx},'.mat'));
    roi = load(roi_fname).roi;
    
    h_image_roi = h_image(roi.x,roi.y);
    dab_image_roi = dab_image(roi.x,roi.y);

    %% Spatial filter
    h_image_f = imfilter(h_image_roi,box_kernel,'replicate');
    dab_image_f = imfilter(dab_image_roi,box_kernel,'replicate');

    % visualize filtered images
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(121),imshow(h_image_roi),colormap(h_colormap),title('Original')
    subplot(122),imshow(h_image_f),colormap(h_colormap),title('Filtered')

    fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(121),imshow(dab_image_roi),colormap(dab_colormap),title('Original')
    subplot(122),imshow(dab_image_f),colormap(dab_colormap),title('Filtered')

    fname1 = strcat(fname(1:end-4),'_1_H_filtered.tif');
    saveas(fig1,fullfile(fname1));
    close(fig1);

    fname2 = strcat(fname(1:end-4),'_2_DAB_filtered.tif');
    saveas(fig2,fullfile(fname2));
    close(fig2);
                
    %% Apply threshold to image to create a mask
    % create a maks from threshold
    dab_roi_mask = ~imbinarize(dab_image_f,pixel_thresh);
%     figure,imshow(dab_roi_mask)

    % overlay image with mask
    dab_image_mask = imoverlay(dab_image_f,dab_roi_mask,[1,0,0]); % red

    % Remove small regions of connected pixels
    min_con_pixels = 20;
    connectivity = 4; % 4
    dab_roi_mask = bwareaopen(dab_roi_mask,min_con_pixels,connectivity);
%     figure,imshow(dab_roi_mask)
    
    fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(121),imshow(dab_roi_mask),title('Mask')
    subplot(122),imshow(dab_image_mask),title('Image + Filtered Mask')
    
    fname3 = strcat(fname(1:end-4),'_mask.tif');
    saveas(fig3,fullfile(fname3));
    close(fig3);
    
    %% Calculate density (% area covered)
    positive_pixels = sum(sum(dab_roi_mask));
    total_pixels = size(dab_roi_mask,1)*size(dab_roi_mask,2);
    density = round(positive_pixels/total_pixels*100,1);

    %% Save results
    results.fname = fname;
    results.pixel_thresh = pixel_thresh;
    results.density = density;
    results.positive_pixels = positive_pixels;
    results.total_pixels = total_pixels;

    fname_r = strcat(file(1:end-4),'_results.mat');
    save(fullfile(results_folder,fname_r),'results');
    fprintf('Results %s saved\n',fname_r);

end


