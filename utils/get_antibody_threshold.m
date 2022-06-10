function [pixel_thresh] = get_antibody_threshold(antibody_type)
    %% Define detection threshold for each antibody
    % @author: pdzialecka
    
    % Threshold can be between 0 and 1, with smaller values corresponding
    % to higher intensity of antibody labelling required
    
    %%
    if strcmp(antibody_type,'moc23')
        pixel_thresh = 0.65; % 160;

    elseif strcmp(antibody_type,'cfos')
        pixel_thresh = 0.2; % decrease to detect only v dark cells

    elseif strcmp(antibody_type,'GFAP')
        pixel_thresh = 0.5;

    elseif strcmp(antibody_type,'Iba1')
        pixel_thresh = 0.35;
    end
    
end
