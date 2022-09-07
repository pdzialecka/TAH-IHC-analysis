function [r_colormap,g_colormap,b_colormap] = create_rgb_colormaps()
    %% Create R, G, B colormaps for easier visualisation of IF images
    % @author: pdzialecka

    %%
    r_col = [1,0,0];
    g_col = [0,1,0];
    b_col = [0,0,1];
    
    r_colormap = [(linspace(0,r_col(1),256)'),(linspace(0,r_col(2),256)'),(linspace(0,r_col(3),256)')];
    g_colormap = [(linspace(0,g_col(1),256)'),(linspace(0,g_col(2),256)'),(linspace(0,g_col(3),256)')];
    b_colormap = [(linspace(0,b_col(1),256)'),(linspace(0,b_col(2),256)'),(linspace(0,b_col(3),256)')];

end
