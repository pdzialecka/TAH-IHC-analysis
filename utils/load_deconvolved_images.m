function [h_image,dab_image,res_image] = load_deconvolved_images(file_path)
    %% Load deconvolved channel images
    % @author: pdzialecka
    
    %%
    image = read_file(file_path);
    h_image = image(:,:,1);
    dab_image = image(:,:,2);
    res_image = image(:,:,3);

    %% Rotate images
    h_image = imrotate(h_image,-90);
    dab_image = imrotate(dab_image,-90);
    res_image = imrotate(res_image,-90);
    
end
