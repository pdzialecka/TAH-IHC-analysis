function [pixel_thresh,min_size,do_watershed] = get_antibody_threshold(antibody_type)
    %% Define detection threshold for each antibody
    % @author: pdzialecka
    
    % Threshold can be between 0 and 1, with smaller values corresponding
    % to higher intensity of antibody labelling required
    
    %%
    if strcmp(antibody_type,'moc23')
        pixel_thresh = 0.8; % 160;
        min_size = 10; % um
        do_watershed = 0;

    elseif strcmp(antibody_type,'12f4')
        pixel_thresh = 0.7;
        min_size = 10;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'ct695')
        pixel_thresh = 0.65;
        min_size = 1;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'iba1')
        pixel_thresh = 0.65;
        min_size = 5;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'gfap')
        pixel_thresh = 0.3;
        min_size = 5; % TODO: adjust
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'cfos')
        pixel_thresh = 0.6; % decrease to detect only v dark cells
        min_size = 5;
        do_watershed = 1;

    elseif strcmp(antibody_type,'ki67')
        pixel_thresh = 0.2;
        min_size = 5; % TODO: adjust
        do_watershed = 1;

    elseif strcmp(antibody_type,'dcx')
        pixel_thresh = 0.3;
        min_size = 5;
        do_watershed = 0;
        
    elseif strcmp(antibody_type,'sox2')
        pixel_thresh = 0.5;
        min_size = 5;
        do_watershed = 0;

    end
end
