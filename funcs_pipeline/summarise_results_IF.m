function [] = summarise_results_IF(base_folder,cohort_case,close_figs)
    %% Create summary plot of results per image type
	% @author: pdzialecka

    %%
    if ~exist('close_figs','var')
        close_figs = 1;
    end
    
    %% H DAB colormaps
    % for easier visualisation
%     [h_colormap,dab_colormap] = create_hdab_colormaps();

    %% Data folders
%     base_folder = 'C:\Users\Pat\Desktop\TAH';
%     cohort_case = 1; % 1 = cohort 1 (13 mo). 2 = cohorts 2-5 (6 mo)

    if cohort_case == 1
        cohort_idxs = [1];
    elseif cohort_case == 2
        cohort_idxs = [2,3,4,5];
    end
    
    all_cohort_folders = dir(fullfile(base_folder,'Data','Cohort*'));
%     data_folders = all_cohort_folders(cohort_idxs);

    include_idxs = [];
    for folder_idx = 1:length(all_cohort_folders)
        cohort_no = str2num(all_cohort_folders(folder_idx).name(end));
        if any(cohort_idxs==cohort_no)
            include_idxs = [include_idxs folder_idx];
        end
    end

    data_folders = all_cohort_folders(include_idxs);

    %% Results folder
    results_folder = fullfile(base_folder,'IF_results');

    % if ~exist(results_folder)
    %     mkdir(results_folder);
    % end

    % Cohort results folder
    if cohort_case == 1
        cohort_results_folder = fullfile(results_folder,'Cohorts_1_13mo');

    elseif cohort_case == 2
        cohort_results_folder = fullfile(results_folder,'Cohorts_2-5_6mo');
    end
    
    img_type = 'iba1_4g8';
    cohort_results_folder = fullfile(cohort_results_folder,img_type);
    if ~exist(cohort_results_folder)
        mkdir(cohort_results_folder)
    end
    
    %% Extra subfolders
    % roi comparison
    comparison_folder = fullfile(cohort_results_folder,'ROI_comparison');
    if ~exist(comparison_folder)
        mkdir(comparison_folder);
    end
    
    % stats
    stats_folder = fullfile(cohort_results_folder,'Stats');
    if ~exist(stats_folder)
        mkdir(stats_folder)
    end    

    %% Load cohort info
    cohort_info_folder = fullfile(results_folder,'Cohort_info');
    cohort_infos = {};

    mouse_ids = [];
    mouse_cond_idxs = [];

    for i = 1:length(cohort_idxs)
        cohort_infos{i} = load(fullfile(cohort_info_folder,...
            sprintf('Cohort_%d_info.mat',cohort_idxs(i)))).cohort;
        mouse_ids = [mouse_ids cohort_infos{i}.mouse_ids];
        mouse_cond_idxs = [mouse_cond_idxs cohort_infos{i}.mouse_cond_idxs];
    end
    
    mouse_names = mouse_ids_to_names(mouse_ids);
    [~,m_idxs] = sort(mouse_names);
    mouse_names = mouse_names(m_idxs);
    mouse_cond_idxs = mouse_cond_idxs(m_idxs);
    mouse_no = length(mouse_ids);
    
    %% Find mouse names per condition
    sham_names = mouse_names(mouse_cond_idxs==1)';
    gamma_names = mouse_names(mouse_cond_idxs==2)';
    theta_names = mouse_names(mouse_cond_idxs==3)';
    ltd_names = mouse_names(mouse_cond_idxs==4)';
    condition_mouse_names = {sham_names,gamma_names,theta_names,ltd_names};

    %% Find result files
    result_files = [];
    roi_img_files = [];
%     roi_img_norm_files = [];
    roi_mask_files = [];
    roi_count_imgs = [];
    roi_area_imgs = [];

    for i = 1:length(data_folders)
        data_folder = fullfile(data_folders(i).folder,data_folders(i).name);
        m_names = cohort_infos{i}.mouse_names;

        for j = 1:length(m_names)
            idx_files = dir(fullfile(data_folder,'IF','Results',m_names{j},strcat('*results.mat')));
            result_files = [result_files; idx_files];
            
            idx_2_files = dir(fullfile(data_folder,'IF','ROI_images',m_names{j},'*',strcat('*.tif')));
            idx_2_files(contains({idx_2_files.name}','Merged')) = [];
            roi_img_files = [roi_img_files; idx_2_files];
%             
%             idx_4_files = dir(fullfile(data_folder,'IHC','ROI_images_norm',m_names{j},strcat('*',img_type,'*.tif')));
%             roi_img_norm_files = [roi_img_norm_files; idx_4_files];
%             
            idx_3_files = dir(fullfile(data_folder,'IF','Results',m_names{j},strcat('*masks.mat')));
            roi_mask_files = [roi_mask_files; idx_3_files];
            
            idx_4_files = dir(fullfile(data_folder,'IF','Results',m_names{j},strcat('*colocalised_count.tif')));
            roi_count_imgs = [roi_count_imgs; idx_4_files];
            
            idx_5_files = dir(fullfile(data_folder,'IF','Results',m_names{j},strcat('*colocalised_area.tif')));
            roi_area_imgs = [roi_area_imgs; idx_5_files];


        end

        % faster but wrong order
    %     idx_files = dir(fullfile(data_folder,'IHC','Results','**',strcat('*',img_type,'*results.mat')));
    %     result_files = [result_files; idx_files];
    end
    
    %% Sort result files
    [~,s_idxs1] = sort({result_files.name});
    [~,s_idxs2] = sort({roi_img_files.name});
%     [~,s_idxs4] = sort({roi_img_norm_files.name});
    [~,s_idxs3] = sort({roi_mask_files.name});
    [~,s_idxs4] = sort({roi_count_imgs.name});
    [~,s_idxs5] = sort({roi_area_imgs.name});


    result_files = result_files(s_idxs1);
    roi_img_files = roi_img_files(s_idxs2);
%     roi_img_norm_files = roi_img_norm_files(s_idxs4);
    roi_mask_files = roi_mask_files(s_idxs3);
    roi_count_imgs = roi_count_imgs(s_idxs4);
    roi_area_imgs = roi_area_imgs(s_idxs5);

%     {result_files.name}'
%     {roi_img_files.name}'
%     {roi_mask_files.name}'

    %% Conditions
    cond_names = {'Sham','40 Hz','8 Hz','LTD'};
    cond_no = length(cond_names);
    variable_names_t = {'Sham','Gamma (40 Hz)','Theta (8 Hz)','LTD (1 Hz)'};
    
    [roi_names,roi_fnames,roi_no] = get_roi_list_IF();
    
    %% Load results
    results = {};
    roi_microglia_ab_ratio = nan(roi_no,mouse_no);
    roi_ab_microglia_ratio = nan(roi_no,mouse_no);
    
    roi_microglia_ab_area_ratio = nan(roi_no,mouse_no);
    roi_ab_microglia_area_ratio = nan(roi_no,mouse_no);

    roi_microglia_per_ab_count = {};
    
    for roi_idx = 1:roi_no
        roi_fname = roi_fnames{roi_idx};
        file_idxs = find(contains({result_files.name}',roi_fname));

        for i = 1:length(file_idxs)
            file_idx = file_idxs(i);

            results{roi_idx,i} = load(fullfile(result_files(file_idx).folder,result_files(file_idx).name)).results;
            
            % count ratio
            roi_microglia_ab_ratio(roi_idx,i) = results{roi_idx,i}.microglia_ab_ratio;
            roi_ab_microglia_ratio(roi_idx,i) = results{roi_idx,i}.ab_microglia_ratio;
            
            % microglia per ab count
            roi_microglia_per_ab_count{roi_idx,i} = results{roi_idx,i}.microglia_per_ab_count;
            
            % area ratio
            roi_microglia_ab_area_ratio(roi_idx,i) = results{roi_idx,i}.microglia_ab_area_ratio;
            roi_ab_microglia_area_ratio(roi_idx,i) = results{roi_idx,i}.ab_microglia_area_ratio;
            
        end
    end

    %% Results summary: percentage of ab+ microglia
    results_ratio_1_all = plot_results(roi_microglia_ab_ratio,'microglia_ratio',...
        mouse_cond_idxs,img_type,cohort_results_folder);

    %% Results summary: percentage of microglia+ ab
    results_ratio_2_all = plot_results(roi_ab_microglia_ratio,'ab_ratio',...
        mouse_cond_idxs,img_type,cohort_results_folder);
    
    %% Results summary: percentage of ab+ microglia (AREA)
    results_ratio_1_area_all = plot_results(roi_microglia_ab_area_ratio,'microglia_area_ratio',...
        mouse_cond_idxs,img_type,cohort_results_folder);

    %% Results summary: percentage of microglia+ ab (AREA)
    results_ratio_2_area_all = plot_results(roi_ab_microglia_area_ratio,'ab_area_ratio',...
        mouse_cond_idxs,img_type,cohort_results_folder);
    
    %% Results summary: number of microglia per ab    
    results_microglia_per_ab_all = plot_results(roi_microglia_per_ab_count,'microglia_per_ab',...
        mouse_cond_idxs,img_type,cohort_results_folder);
    
    %% Statistics
    % count
    [p_ratio_1,h_ratio_1] = compute_stats(results_ratio_1_all,'microglia_ratio',...
                            img_type,stats_folder);
                        
    [p_ratio_2,h_ratio_2] = compute_stats(results_ratio_2_all,'ab_ratio',...
                            img_type,stats_folder);
                        
    % overlap area
    [p_ratio_1_area,h_ratio_1_area] = compute_stats(results_ratio_1_area_all,'microglia_area_ratio',...
                            img_type,stats_folder);
                        
    [p_ratio_2_area,h_ratio_2_area] = compute_stats(results_ratio_2_area_all,'ab_area_ratio',...
                            img_type,stats_folder);
                      
                        
    % microglia per ab
    [p_ratio_3,h_ratio_3] = compute_stats(results_microglia_per_ab_all,'microglia_per_ab',...
                            img_type,stats_folder);

    %% Plot result images for comparison
    max_n = 7;
    
    for roi_idx = 1%:roi_no
        
        roi_fname = roi_fnames{roi_idx};
        img_idxs = find(contains({roi_img_files.name}',roi_fname));
        
        roi_results_ratio_1 = results_ratio_1_all{roi_idx};
        roi_results_ratio_2 = results_ratio_2_all{roi_idx};
        roi_results_ratio_1_area = results_ratio_1_area_all{roi_idx};
        roi_results_ratio_2_area = results_ratio_2_area_all{roi_idx};
        
        comparison_subfolder = fullfile(comparison_folder,roi_fname);
        if ~exist(comparison_subfolder,'dir')
            mkdir(comparison_subfolder);
        end


        
%         if isempty(results_count_all)
%             roi_results_count = nan(size(roi_results_ratio));
%         else
%             roi_results_count = results_count_all{roi_idx};
%         end
        
 
        for j = 1:cond_no
%             fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
%             fig2 = figure('units','normalized','outerposition',[0 0 1 1]);

%             comparison_subfolder = fullfile(comparison_folder,cond_names{j});
%             if ~exist(comparison_subfolder,'dir')
%                 mkdir(comparison_subfolder);
%             end

            cond_j_names = condition_mouse_names{j};
            cond_img_idxs = img_idxs((contains({roi_img_files(img_idxs).name}',cond_j_names)));
%             {roi_img_files(cond_img_idxs).name}'

            % remove nan values
%             roi_results_density_j = roi_results_ratio(:,j);
%             roi_results_count_j = roi_results_count(:,j);
            
            roi_results_ratio_1_j = roi_results_ratio_1(:,j);
            roi_results_ratio_2_j = roi_results_ratio_2(:,j);
            roi_results_ratio_1_area_j = roi_results_ratio_1_area(:,j);
            roi_results_ratio_2_area_j = roi_results_ratio_2_area(:,j);

%             
%             missing_idxs = isnan(roi_results_density_j);
%             roi_results_density_j(missing_idxs) = [];
%             roi_results_count_j(missing_idxs) = [];
        
            
            for i = 1:length(cond_img_idxs)
                
                % composite mask + cells detected
                fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
                roi_count_img = imread(fullfile(roi_count_imgs(cond_img_idxs(i)).folder,roi_count_imgs(cond_img_idxs(i)).name));                
                set(0,'CurrentFigure',fig1)
%                 subplot(round(max_n/2),2,i),
                imshow(roi_count_img)
%                 title_txt = {roi_count_imgs(cond_img_idxs(i)).name,...
%                     sprintf('Microglia ab+ = %1.2f%%. Ab microglia+ = %1.2f%%\n',roi_results_ratio_1_j(i),roi_results_ratio_2_j(i))};
                title_txt = sprintf('Microglia ab+ = %1.2f%%. Ab microglia+ = %1.2f%%',roi_results_ratio_1_j(i),roi_results_ratio_2_j(i));
                title(title_txt,'Interpreter','none')
            
                
                % composite mask + overlap
                fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
                roi_area_img = imread(fullfile(roi_area_imgs(cond_img_idxs(i)).folder,roi_area_imgs(cond_img_idxs(i)).name));
                set(0,'CurrentFigure',fig2)
%                 subplot(round(max_n/2),2,i),
                imshow(roi_area_img)
%                 title_txt = {roi_area_imgs(cond_img_idxs(i)).name,...
%                     sprintf('Microglia ab+ area = %1.2f%%. Ab microglia+ area = %1.2f%%\n',roi_results_ratio_1_area_j(i),roi_results_ratio_2_area_j(i))};
                title_txt = sprintf('Microglia ab+ area = %1.2f%%. Ab microglia+ area = %1.2f%%',roi_results_ratio_1_area_j(i),roi_results_ratio_2_area_j(i));
                title(title_txt,'Interpreter','none')

                
                fig_1_name = sprintf('%s_%d_%s.tif',roi_count_imgs(cond_img_idxs(i)).name(1:end-4),j,cond_names{j});
                saveas(fig1,fullfile(comparison_subfolder,fig_1_name));

                fig_2_name = sprintf('%s_%d_%s.tif',roi_area_imgs(cond_img_idxs(i)).name(1:end-4),j,cond_names{j});
                saveas(fig2,fullfile(comparison_subfolder,fig_2_name));
            
                if close_figs
                    close(fig1);
                    close(fig2);
                end
            end
            
            
%             fig_1_name = sprintf('roi_images_%s_%d_%s_%d_%s_colocalised_count.tif',img_type,roi_idx,roi_fnames{roi_idx},j,cond_names{j});
%             saveas(fig1,fullfile(comparison_folder,fig_1_name));
% 
%             fig_2_name = sprintf('roi_images_%s_%d_%s_%d_%s_colocalised_area.tif',img_type,roi_idx,roi_fnames{roi_idx},j,cond_names{j});
%             saveas(fig2,fullfile(comparison_folder,fig_2_name));
%             
%             
%             if close_figs
%                 close(fig1);
%                 close(fig2);
%             end

        end
        
    end
    
end
