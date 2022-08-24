function [] = deconvolve_full(files,save_images)
    %% Deconvolve IHC images to separate H and DAB channels
    % @author: pdzialecka
    
    % Deconvolve H DAB images using standard function in Fiji
    % Save results as tiff for easy access to entire FOV for further ROI
    % selection
    
    % The script establishes communication between Matlab and Fiji using
    % MIJI and IJ which can be used in further scripts to access plugins
    % available in Fiji
    
    %% Input variables
    if ~exist('save_images','var')
        save_images = 1;
    end

    %% Load libraries
    addpath(genpath('Libraries'));

    %% Add FIJI / MIJ paths
    setup_miji();
%     % Change these depending on local directories
%     javaaddpath 'C:\Program Files\MATLAB\R2021a\java\jar\mij.jar'
%     javaaddpath 'C:\Users\Pat\Desktop\Fiji.app\jars\ij-1.53q.jar'
%     addpath(genpath('C:\Users\Pat\Desktop\Fiji.app\scripts'));
    
    %% Process all files in the folder
    
    %% Start fiji
    Miji;
    IJ = ij.IJ;
    
    % Alternatives tried:
    % ij.ImageJ()
    % ImageJ;
    
    % MIJ.start;
    % MIJ.setupExt('C:\Program Files\MATLAB\R2021a\java\jar');

    %% Log command window messages for any errors
%     diary_file = fullfile(fileparts(fileparts(files(1).folder)),'deconvolution_log.txt');
%     diary diary_file
    
    %%
%     try
    for idx = 1:length(files)
        
        try
            %% File info
            file = files(idx).name;
            folder = files(idx).folder;
    %         cd(folder)

            %% Prepare save folder
    %             save_folder = fullfile(folder,'Processed');
            save_folder = fullfile(fileparts(fileparts(folder)),'Images_processed');

            % create subfolder with mouse name
            [~,mouse_name] = fileparts(folder);
            save_folder = fullfile(save_folder,mouse_name);

            if ~exist(save_folder)
                mkdir(save_folder);
            end
            
            %% Check if the file exists
            fname_2 = fullfile(save_folder,strcat(file(1:end-4),'_deconv.tif'));
            
            %%
            if ~exist(fname_2,'file')

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
                if save_images
                    fname_1 = fullfile(save_folder,strcat(file(1:end-4),'.tif'));
                    IJ.save(fname_1)
                end
                MIJ.run("Close")


                % this ensures file saved as 8-bit
                MIJ.selectWindow("Colour Deconvolution")
                MIJ.run("Close")


                % Save deconvolved images as a stack
                MIJ.run("Images to Stack", "use");
                MIJ.run("Grays");
                % MIJ.run("Delete Slice"); % remove color info slice -> this will save img as RGB

                if save_images
                    fname_2 = fullfile(save_folder,strcat(file(1:end-4),'_deconv.tif'));
                    IJ.save(fname_2)
                end

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
            
        catch
            fprintf('ERROR: File %s was not deconvolved\n',file);
        end
    end

    %% Close Fiji
    MIJ.exit;
    
%     catch
%         MIJ.closeAllWindows;
%         MIJ.exit;
%     end

    %% Save diary log
%     diary off
    
end

