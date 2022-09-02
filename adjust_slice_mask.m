%%
close all
clear all

%% Analysis folder
analysis_folder = 'C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC';
addpath(genpath(analysis_folder));

%%
[h_colormap,dab_colormap] = create_hdab_colormaps();

%% Specify file to adjust
% location of the roi tif image
file_path = 'C:\Users\Pat\Desktop\TAH\Data\Cohort_3\IHC\ROI_images\AD-Hipp39\AD-Hipp39_cfos_8_roi_cortex_R.tif';

%% Specify number of areas to remove
area_no = 1;

%% Find all directories
[roi_img_folder,roi_fname] = fileparts(file_path);
[a,mouse_name] = fileparts(roi_img_folder);
tmp = strfind(roi_fname,'_');
roi_fname_core = roi_fname(1:tmp(2)-1);
slice_mask_fname = strcat(roi_fname_core,'_slice_mask.mat');
roi_folder = fullfile(fileparts(fileparts(roi_img_folder)),'ROIs',mouse_name);

%% Load roi image
[h_image_roi,dab_image_roi] = load_deconvolved_images(file_path);

%% Load roi info
roi = load(fullfile(roi_folder,roi_fname)).roi;

%% Load overall slice mask
slice_mask = load(fullfile(roi_folder,slice_mask_fname)).slice_mask;

%% Extract roi location
roi_x_1 = roi.x(1);
roi_y_1 = roi.y(1);
roi_size = roi.size;

roi_x = roi_x_1:roi_x_1+roi_size(1);
roi_y = roi_y_1:roi_y_1+roi_size(2);


%% Select areas to remove
slice_mask_roi = slice_mask(roi_y,roi_x);
remove_mask = zeros(size(slice_mask));
remove_mask_roi = zeros(size(slice_mask_roi));

for idx = 1:area_no

    %% Select area to remove from the slice mask
    figure,imshow(h_image_roi),colormap(h_colormap)
    title('Select area to remove')
    remove_area = drawpolygon();
    title('Area to remove selected')

    %% Create new remove mask
    remove_mask_roi = remove_mask_roi | poly2mask(remove_area.Position(:,1),remove_area.Position(:,2),size(dab_image_roi,1),size(dab_image_roi,2));
    slice_mask_roi(remove_mask_roi) = 0;

    %% Create remove mask
    remove_mask(roi_y,roi_x) = remove_mask_roi; 
    remove_mask = logical(remove_mask);

end

%% Visualise roi image with slice and remove masks
% figure,imshow(~slice_mask_roi)

h_image_slice_mask = labeloverlay(h_image_roi,~slice_mask_roi,...
    'Colormap',[0,0,1],'Transparency',0.7);
h_image_slice_mask_2 = labeloverlay(h_image_slice_mask,remove_mask_roi,...
    'Colormap',[1,0,0],'Transparency',0.7);
fig = figure;
imshow(h_image_slice_mask_2)

dab_image_slice_mask = labeloverlay(dab_image_roi,~slice_mask_roi,...
    'Colormap',[0,0,1],'Transparency',0.7);
dab_image_slice_mask_2 = labeloverlay(dab_image_slice_mask,remove_mask_roi,...
    'Colormap',[1,0,0],'Transparency',0.7);
fig2 = figure;
imshow(dab_image_slice_mask_2)

fig_name = strcat(roi_fname,'_remove_mask.tif');
saveas(fig,fullfile(roi_folder,fig_name));

%%
close(fig);
% close(fig2);

%%
% figure,imshow(remove_mask)

%% Save remove mask
remove_mask_fname = strcat(roi_fname,'_remove_mask.mat');
save(fullfile(roi_folder,remove_mask_fname),'remove_mask');
