%%
close all
clear all

%%
[h_colormap,dab_colormap] = create_hdab_colormaps();

%% Specify file to adjust
% location of the roi tif image
file_path = 'C:\Users\Pat\Desktop\TAH\Data\Cohort_5\IHC\ROI_images\AD-Hipp52\AD-Hipp52_moc23_7_roi_hipp_CA3_L.tif';

%% Find all directories
[roi_img_folder,roi_fname] = fileparts(file_path);
[a,mouse_name] = fileparts(roi_img_folder);
tmp = strfind(roi_fname,'_');
roi_fname_core = roi_fname(1:tmp(2)-1);
auto_rois_fname = strcat(roi_fname_core,'_rois_auto.mat');
slice_mask_fname = strcat(roi_fname_core,'_slice_mask.mat');
roi_folder = fullfile(fileparts(fileparts(roi_img_folder)),'ROIs',mouse_name);

%% Load roi image
[h_image_roi] = load_deconvolved_images(file_path);

%% Load roi info
rois = load(fullfile(roi_folder,auto_rois_fname));

% rois = load('C:\Users\Pat\Desktop\TAH\Data\Cohort_5\IHC\ROIs\AD-Hipp52\AD-Hipp52_moc23_rois_auto.mat');
% roi_fname = 'AD-Hipp52_moc23_7_roi_hipp_CA3_L';
[coords,this_roi,coords_field] = extract_roi_coords(rois,roi_fname);

%% Load overall slice mask
slice_mask = load(fullfile(roi_folder,slice_mask_fname)).slice_mask;

%% Extract roi location
roi_x_1 = coords(1);
roi_y_1 = coords(2);
roi_size = coords(3:4);

roi_x = roi_x_1:roi_x_1+roi_size(1);
roi_y = roi_y_1:roi_y_1+roi_size(2);

%% Select area to remove from the slice mask
figure,imshow(h_image_roi),colormap(h_colormap)
title('Select area to remove')
remove_area = drawpolygon();
title('Area to remove selected')

%% Create new remove mask
remove_mask_roi = poly2mask(remove_area.Position(:,1),remove_area.Position(:,2),size(h_image_roi,1),size(h_image_roi,2));

%% Visualise roi image with slice and remove masks
slice_mask_roi = slice_mask(roi_y,roi_x);
% figure,imshow(~slice_mask_roi)

slice_mask_roi(remove_mask_roi) = 0;
% figure,imshow(~slice_mask_roi)

h_image_slice_mask = labeloverlay(h_image_roi,~slice_mask_roi,...
    'Colormap',[0,0,1],'Transparency',0.7);
h_image_slice_mask_2 = labeloverlay(h_image_slice_mask,remove_mask_roi,...
    'Colormap',[1,0,0],'Transparency',0.7);
fig = figure;
imshow(h_image_slice_mask_2)

fig_name = strcat(roi_fname,'_remove_mask.tif');
saveas(fig,fullfile(roi_folder,fig_name));
close(fig);

%%
% slice_mask(roi_y,roi_x) = slice_mask_roi;
% figure,imshow(slice_mask)

%% Create remove mask
remove_mask = zeros(size(slice_mask));
remove_mask(roi_y,roi_x) = remove_mask_roi; 
remove_mask = logical(remove_mask);

figure,imshow(remove_mask)

%% Save remove mask
remove_mask_fname = strcat(roi_fname_core,'_remove_mask.mat');
save(fullfile(roi_folder,remove_mask_fname),'remove_mask');
