function [size_in_pixels] = um_to_pixel(size_in_um)
    %%
    % @author: pdzialecka
    
    %%
    pixel_size = 0.504;
    
    size_in_pixels = round(size_in_um/pixel_size);
    
end
