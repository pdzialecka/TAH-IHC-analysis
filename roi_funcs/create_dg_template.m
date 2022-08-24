%% Create DG template
%% Load img
img = imread('C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC\Allen_atlas\dg_img.png');
img = rgb2gray(img);

%% Create a smooth mask
dg_template = img<100;

% spatial filter
k = 10;
kernel = 1/(k*k)*ones([k,k]);
dg_template = imfilter(dg_template,kernel,'replicate');

% clean mask
dg_template = bwareaopen(dg_template,1000,4);

% remove background
regions = regionprops(dg_template,'Area','Centroid','BoundingBox');
bb = round(regions.BoundingBox);
dg_template = dg_template(bb(2):bb(2)+bb(4),bb(1):bb(1)+bb(3));

%% Rescale to approx img dg width
% w = 2440; % ~2500 pixels
rfactor = 10;
dg_template = imresize(dg_template,rfactor);

%% Save
save('C:\Users\Pat\OneDrive - Imperial College London\AD_TI_hipp\Analysis\IHC\Allen_atlas\dg_template.mat','dg_template');

%%
template = 0.2*ones(11);
template(6,3:9) = 0.6;
template(3:9,6) = 0.6;
offsetTemplate = 0.2*ones(22);
offset = [8 6];
offsetTemplate((1:size(template,1))+offset(1), ...
    (1:size(template,2))+offset(2)) = template;

cc = xcorr2(offsetTemplate,template);
[max_cc, imax] = max(abs(cc(:)));
[ypeak, xpeak] = ind2sub(size(cc),imax(1));
corr_offset = [(ypeak-size(template,1)) (xpeak-size(template,2))];

isequal(corr_offset,offset)
