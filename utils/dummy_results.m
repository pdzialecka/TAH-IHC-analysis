%% Create dummy result files for missing files
mouse_names = {'AD-Hipp41','AD-Hipp42','AD-Hipp43','AD-Hipp44','AD-Hipp45'};

for midx = 1:length(mouse_names)
    mouse_name = mouse_names{midx};
    
    fnames = {sprintf('%s_dcx_1_roi_hipp_DG_L_results.mat',mouse_name),...
              sprintf('%s_dcx_2_roi_hipp_DG_R_results.mat',mouse_name),...
              sprintf('%s_dcx_3_roi_hipp_CA1_L_results.mat',mouse_name),...
              sprintf('%s_dcx_4_roi_hipp_CA1_R_results.mat',mouse_name),...
              sprintf('%s_dcx_5_roi_hipp_CA3_L_results.mat',mouse_name),...
              sprintf('%s_dcx_6_roi_hipp_CA3_R_results.mat',mouse_name),...
              sprintf('%s_dcx_7_roi_cortex_L_results.mat',mouse_name),...
              sprintf('%s_dcx_8_roi_cortex_R_results.mat',mouse_name)};

    for i = 1:length(fnames)
        fname = fnames{i};

        results = [];

        results.fname = fname;
        results.thresh_pixel = nan;

        % info on slice mask
        results.slice_area = nan;
        results.total_area = nan;
        results.slice_area_norm = nan;

        % main results
        results.density = nan;
        results.positive_pixels = nan;
        results.total_pixels = nan;
        results.particles_found = nan;
        results.particle_no = nan;

        fname_r = strcat(fname,'_results.mat');
        save(fullfile(fname_r),'results');
        fprintf('Results %s saved\n',fname_r);

    end

end
