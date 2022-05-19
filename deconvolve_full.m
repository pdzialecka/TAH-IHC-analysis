function deconvolve_full(files)
    %% Deconvolve IHC images to separate H and DAB channels
    % @author: pdzialecka
    
    % Deconvolve H DAB images using standard function in Fiji
    % Save results as tiff for easy access to entire FOV for further ROI
    % selection
    
    % The script establishes communication between Matlab and Fiji using
    % MIJI and IJ which can be used in further scripts to access plugins
    % available in Fiji
    
    %%
    close all
    clear all
    
    %% Input variables

    %% Load libraries
    addpath(genpath('Libraries'));

    %% Add FIJI / MIJ paths
    % Change these depending on local directories
    javaaddpath 'C:\Program Files\MATLAB\R2021a\java\jar\mij.jar'
    javaaddpath 'C:\Users\Pat\Desktop\Fiji.app\jars\ij-1.53q.jar'
    addpath(genpath('C:\Users\Pat\Desktop\Fiji.app\scripts'));

    %% Start fiji
    Miji;
    IJ = ij.IJ;
    
    % Alternatives tried:
    % ij.ImageJ()
    % ImageJ;
    
    % MIJ.start;
    % MIJ.setupExt('C:\Program Files\MATLAB\R2021a\java\jar');

    %%
    for idx = 1:length(files)
        %%
        file = files(idx).name;
        folder = files(idx).folder;
%         cd(folder)

        %% Load RGB svs file
        fname = fullfile(folder,file);
        image_RGB = imread(fname);

        % Save as tiff
        % imwrite(a,strcat(fname(1:end-4),'.tif'));

        %% Open file in Fiji
        imp = copytoImagePlus(image_RGB);
        imp.show();

        % MIJ.createColor('img_small') % doesn't work for an unknown reason
        % MIJ.createImage(img_small); % grayscale only

        %% Color deconvolution in Fiji
        MIJ.run("RGB Color");
        imp.close();

        MIJ.run("Colour Deconvolution", "vectors=[H DAB]");

        %% Save results
        % Save combined tif file
        MIJ.selectWindow("new (RGB)")
        fname_tif = strcat(fname(1:end-4),'.tif');
        if save_files
            IJ.save(fname_tif)
        end
        MIJ.run("Close")


        % this ensures file saved as 8-bit
        MIJ.selectWindow("Colour Deconvolution")
        MIJ.run("Close")


        % Save deconvolved images as a stack
        MIJ.run("Images to Stack", "use");
        MIJ.run("Grays");
        % MIJ.run("Delete Slice"); % remove color info slice -> this will save img as RGB

        fname_deconv_tif = strcat(fname(1:end-4),'_deconv.tif');
        IJ.save(fname_deconv_tif)

        %% Alternatively, send images back to Matlab and save from here
        % very slow for large images (full frame)
        % image_h = MIJ.getImage("new (RGB)-(Colour_1)");
        % image_dab = MIJ.getImage("new (RGB)-(Colour_2)");
        % image_res = MIJ.getImage("new (RGB)-(Colour_3)");

        % image_h = uint8(image_h);
        % image_dab = uint8(image_dab);
        % image_res = uint8(image_res);

        %% Close all windows
        MIJ.closeAllWindows;
    end

    %% Close Fiji
    MIJ.exit;

end
