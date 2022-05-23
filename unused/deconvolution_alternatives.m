%% %%%%%%%% Deconvolution Matlab solutions found %%%%%%%%


%% Method 1
% Program: ColorDeconvolutionMatlab-master
% Source: https://github.com/jnkather/ColorDeconvolutionMatlab
% Issue: images look different than in fiji: higher contrast + h channel
% has darker background

%%
% alternative set of standard values (HDAB from Fiji)
He = [ 0.6500286;  0.704031;    0.2860126 ];
DAB = [ 0.26814753;  0.57031375;  0.77642715];
Res = [ 0.7110272;   0.42318153; 0.5615672 ]; % residual

% combine stain vectors to deconvolution matrix
HDABtoRGB = [He/norm(He) DAB/norm(DAB) Res/norm(Res)]';
RGBtoHDAB = inv(HDABtoRGB);
    
% separate stains = perform color deconvolution
tic
imageHDAB = SeparateStains(imageRGB,RGBtoHDAB);
toc

imageHDAB = uint8(round(imageHDAB*255));
figure,imshow(imageHDAB(:,:,2)),colormap(dab_colormap)
figure,imshow(imageHDAB(:,:,1)),colormap(h_colormap)

% show images
fig1 = figure();
set(gcf,'color','w');
subplot(2,4,1); imshow(imageRGB); title('Original');
subplot(2,4,2); imshow(imageHDAB(:,:,1),[]); title('Hematoxylin');
subplot(2,4,3); imshow(imageHDAB(:,:,2),[]); title('DAB');
subplot(2,4,4); imshow(imageHDAB(:,:,3),[]); title('Residual');

subplot(2,4,5); imhist(rgb2gray(imageRGB)); title('Original');
subplot(2,4,6); imhist(imageHDAB(:,:,1)); title('Hematoxylin');
subplot(2,4,7); imhist(imageHDAB(:,:,2)); title('DAB');
subplot(2,4,8); imhist(imageHDAB(:,:,3)); title('Residual');

%% Method 2
% Program: ColorDeconvolutionMatlab-master
% Source: https://github.com/landinig/IJ-Colour_Deconvolution2
% Info: https://blog.bham.ac.uk/intellimic/g-landini-software/colour-deconvolution-2/
% Paper: https://doi.org/10.1093/bioinformatics/btaa847

% Issue: v slow for large images (for loops). Not sure about quality of
% results

%%
% INPUT PARAMETERS:
StainingVectorID = 4;
DyeToBeRemovedID = 0;
doIcross = 0;

% IMAGE RESHAPING:
figure, imshow(uint8(imageRGB), [], 'Border', 'Tight')
ImgR = imageRGB(:,:,1);
ImgG = imageRGB(:,:,2);
ImgB = imageRGB(:,:,3);

% CALL TO THE MAIN FUNCTION:
[ImgR_back, ImgG_back, ImgB_back, Dye01_transmittance, Dye02_transmittance, Dye03_transmittance, LUTdye01, LUTdye02, LUTdye03, Q3x3Mat] = Colour_Deconvolution2(ImgR,ImgG,ImgB, StainingVectorID, DyeToBeRemovedID, doIcross);

% OUTPUT VISUALIZATION: RGB IMAGE RECONSTRUCTED BY REMOVING A DYE
ImgRGB_back(:,:,1) = ImgR_back;
ImgRGB_back(:,:,2) = ImgG_back;
ImgRGB_back(:,:,3) = ImgB_back;
figure,imshow(uint8(ImgRGB_back), 'Border', 'Tight')

% OUTPUT VISUALIZATION: SINGLE DYES: TRANSMITTANCE CHANNELS
figure,imshow(uint8(Dye01_transmittance), 'Border', 'Tight')
figure,imshow(uint8(Dye02_transmittance), 'Border', 'Tight')
figure,imshow(uint8(Dye03_transmittance), 'Border', 'Tight')
