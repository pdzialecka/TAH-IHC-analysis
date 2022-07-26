function [h_image,dab_image,res_image] = load_deconvolved_images(file_path,rotate_img)
    %% Load deconvolved channel images
    % @author: pdzialecka
    
    %%
    if ~exist('rotate_img','var')
        rotate_img = 0;
    end
    
    %%
    image = read_file(file_path);
    h_image = image(:,:,1);
    dab_image = image(:,:,2);
    res_image = image(:,:,3);

    %% Rotate images
    if rotate_img
        h_image = imrotate(h_image,-90);
        dab_image = imrotate(dab_image,-90);
        res_image = imrotate(res_image,-90);
    end
end
