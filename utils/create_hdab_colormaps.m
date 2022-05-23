function [h_colormap,dab_colormap] = create_hdab_colormaps()
    %% Create H DAB colormaps for easier visualisation
    % @author: pdzialecka

    %%
    dab_col = (1-[0.26814753,0.57031375,0.77642715]);
    h_col = (1-[0.6500286,0.704031,0.2860126]);

    dab_colormap = [(linspace(dab_col(1),1,256)'),(linspace(dab_col(2),1,256)'),(linspace(dab_col(3),1,256)')];
    h_colormap = [(linspace(h_col(1),1,256)'),(linspace(h_col(2),1,256)'),(linspace(h_col(3),1,256)')];

end
