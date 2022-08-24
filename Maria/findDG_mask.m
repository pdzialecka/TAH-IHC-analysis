
function [dg_regions,dg_centroids] = findDG_mask(roi_mask)
%smooth
windowSize = 51;
kernel = ones(windowSize) / windowSize ^ 2;
blurryImage = conv2(single(roi_mask), kernel, 'same');
roi_mask = blurryImage > 0.5; % Rethreshold

%Remove too small and too large regions
roi_mask = bwareafilt(roi_mask,[5e+05 4e+06]);

roi_mask = imfill(roi_mask,'holes');
maskBoundary = bwboundaries(roi_mask);

[rows, columns, ~] = size(roi_mask);

%Load DG template outlines
load('DG_outline.mat', 'rDG_x', 'rDG_y', 'lDG_x', 'lDG_y')

%Loop through regions
for i = 1:length(maskBoundary)
    x = maskBoundary{i, 1}(:,2);
    y = maskBoundary{i, 1}(:,1);

    %remove vertical lines
    bw_region = poly2mask(x,y,rows,columns);

    [rws,clms] = size(bw_region);
    
    l = zeros(rws,clms);
    for xi = 1:1
        se = strel('line',100,xi);
        l = l+imopen(bw_region,se);
    end
      
    %Get new mask boundaries
    maskBoundary2 = bwboundaries(l);

    x = maskBoundary2{1, 1}(:,2);
    y = maskBoundary2{1, 1}(:,1);

    %Get center of each region
    centers(i,:) = [((max(x) + min(x))/2) ((max(y) + min(y))/2)];

   if length(rDG_x) < length(x) %Make same length for left DG
        rDG_x2 = interp1((1:length(rDG_x)),rDG_x,(linspace(1,length(rDG_x),length(x))))';
        rx2 = x;
        
        rDG_y2 = interp1((1:length(rDG_y)),rDG_x,(linspace(1,length(rDG_y),length(y))))';
        ry2 = y;
    else
        rx2 = interp1((1:length(x)),x,(linspace(1,length(x),length(rDG_x))))';
        rDG_x2 = rDG_x;
        
        ry2 = interp1((1:length(y)),y,(linspace(1,length(y),length(rDG_y))))';
        rDG_y2 = rDG_y;
    end
    if length(lDG_x) < length(x) %Make same length for left DG
        lDG_x2 = interp1((1:length(lDG_x)),lDG_x,(linspace(1,length(lDG_x),length(x))))';
        lx2 = x;
        
        lDG_y2 = interp1((1:length(lDG_y)),lDG_y,(linspace(1,length(lDG_y),length(y))))';
        ly2 = y;
    else
        lx2 = interp1((1:length(x)),x,(linspace(1,length(x),length(lDG_x))))';
        lDG_x2 = lDG_x;
        
        ly2 = interp1((1:length(y)),y,(linspace(1,length(y),length(lDG_y))))';
        lDG_y2 = lDG_y;
    end

    %Correlate x axis
    c_rDGx(i,1) = max(corr((rDG_x2),(rx2)));
    c_lDGx(i,1) = max(corr((lDG_x2),(lx2)));

    %Correlate y axis
    c_rDGy(i,1) = max(corr((rDG_y2),(ry2)));
    c_lDGy(i,1) = max(corr((lDG_y2),(ly2)));

    c_rDG(i,1) = (c_rDGx(i,1) + c_rDGy(i,1))/2;
    c_lDG(i,1) = (c_lDGx(i,1) + c_lDGy(i,1))/2;
end

 %remove regions on the right/left
 c_rDG(centers(:,1)<columns/2,1) = 0;
 c_lDG(centers(:,1)>columns/2,1) = 0;
 %Get right and left DG regions
 rightDG_mask = maskBoundary{c_rDG == max(c_rDG)};
 leftDG_mask = maskBoundary{c_lDG == max(c_lDG)};

 %Get properties
 bw_lDG = poly2mask(leftDG_mask(:,2),leftDG_mask(:,1),rows,columns);

 regions_l = regionprops(bw_lDG,'Area','Centroid',...
        'BoundingBox','Circularity','Eccentricity',...
        'MajorAxisLength','MinorAxisLength','ConvexArea');
  dg_regions(1) = regions_l([regions_l.Area] == max([regions_l.Area]));

 bw_rDG = poly2mask(rightDG_mask(:,2),rightDG_mask(:,1),rows,columns);

 regions_r = regionprops(bw_rDG,'Area','Centroid',...
        'BoundingBox','Circularity','Eccentricity',...
        'MajorAxisLength','MinorAxisLength','ConvexArea');
 dg_regions(2) = regions_r([regions_r.Area] == max([regions_r.Area]));

 %Get centroids
 dg_centroids = cat(1,dg_regions.Centroid);

