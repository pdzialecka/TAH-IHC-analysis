function [mid_point] = find_mid_point(points)
    %%
    % @author: pdzialecka
    
    % Input: points = [point_1_x, point_1_y; point_2_x, point_2_y]
    
    %%
    dg_dist = points(2,:)-points(1,:);
    mid_point = points(1,:)+dg_dist/2;
    
    mid_point = round(mid_point);
    
end
