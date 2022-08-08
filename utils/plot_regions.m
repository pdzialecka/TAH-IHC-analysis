function [fig] = plot_regions(image,regions)
    %%
    % @author: pdzialecka
    
    %%
    centroids = round(cat(1,regions.Centroid));
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    imshow(image),hold on

    for idx = 1:length(regions)
        h = rectangle('Position',regions(idx).BoundingBox);
        set(h,'EdgeColor','r','LineWidth',1.5);
        plot(centroids(idx,1),centroids(idx,2),'r.',...
            'MarkerSize',10)
        1;
    end
    
end
