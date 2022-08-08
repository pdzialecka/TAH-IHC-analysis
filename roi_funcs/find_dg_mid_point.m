function [mid_point] = find_dg_mid_point(dg_centroids)
    %%
    % @author: pdzialecka
    
    %%
    dg_dist = dg_centroids(2,:)-dg_centroids(1,:);
    mid_point = dg_centroids(1,:)+dg_dist/2;
    
    mid_point = round(mid_point);
    
end
