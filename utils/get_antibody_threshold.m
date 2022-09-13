function [pixel_thresh,min_size,max_size,do_watershed,correct_brightness] = get_antibody_threshold(antibody_type)
    %% Define detection threshold for each antibody
    % @author: pdzialecka
    
    % Threshold can be between 0 and 1, with smaller values corresponding
    % to higher intensity of antibody labelling required
    
    %%
    correct_brightness = 1;
    
    if strcmp(antibody_type,'moc23')
        pixel_thresh = 0.65;
        min_size = 6; % um
        max_size = 100;
        do_watershed = 0;

    elseif strcmp(antibody_type,'12f4')
        pixel_thresh = 0.8; % 0.7 without brigthness correction
        min_size = 6;
        max_size = 100;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'ct695')
        pixel_thresh = 0.7;
        min_size = 1;
        max_size = 15;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'iba1')
        pixel_thresh = 0.5;
        min_size = 4;
        max_size = 100;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'gfap')
        pixel_thresh = 0.45; % 0.3
        min_size = 4;
        max_size = 100;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'cfos')
        pixel_thresh = 0.3; % decrease to detect only v dark cells
        min_size = 4;
        max_size = 30;
        do_watershed = 1;

    elseif strcmp(antibody_type,'ki67')
        pixel_thresh = 0.6; % 0.8 with brightness correction
        min_size = 3;
        max_size = 30;
        do_watershed = 1;
        correct_brightness = 0;

    elseif strcmp(antibody_type,'dcx')
        pixel_thresh = 0.5; % 0.3 without brightness correction
        min_size = 3;
        max_size = 30;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'sox2')
        pixel_thresh = 0.6;
        min_size = 3;
        max_size = 30;
        do_watershed = 1;

    end
end
