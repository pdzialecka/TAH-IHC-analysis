function[tform] = create_tform(theta)
    % @author: pdzialecka

    trans = [0,0];
    rot = [cosd(theta),sind(theta);...
        -sind(theta),cosd(theta);];
    tform = rigid2d(rot,trans);
    
end

