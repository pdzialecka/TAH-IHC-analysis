function [] = summarise_results(base_folder,cohort_case,img_type,close_figs)
    %% Create summary plot of results per image type
	% @author: pdzialecka

    %%
    if ~exist('close_figs','var')
        close_figs = 1;
    end
    
    %% H DAB colormaps
    % for easier visualisation
    [h_colormap,dab_colormap] = create_hdab_colormaps();

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
    results_folder = fullfile(base_folder,'IHC_results');

    % if ~exist(results_folder)
    %     mkdir(results_folder);
    % end

    % Cohort results folder
    if cohort_case == 1
        cohort_results_folder = fullfile(results_folder,'Cohorts_1_13mo');

    elseif cohort_case == 2
        cohort_results_folder = fullfile(results_folder,'Cohorts_2-5_6mo');
    end
    
    % add image_type subfolder
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
    roi_img_norm_files = [];
    roi_mask_files = [];

    for i = 1:length(data_folders)
        data_folder = fullfile(data_folders(i).folder,data_folders(i).name);
        m_names = cohort_infos{i}.mouse_names;

        for j = 1:length(m_names)
            idx_files = dir(fullfile(data_folder,'IHC','Results',m_names{j},strcat('*',img_type,'*results.mat')));
            result_files = [result_files; idx_files];
            
            idx_2_files = dir(fullfile(data_folder,'IHC','ROI_images',m_names{j},strcat('*',img_type,'*.tif')));
            roi_img_files = [roi_img_files; idx_2_files];
            
            idx_4_files = dir(fullfile(data_folder,'IHC','ROI_images_norm',m_names{j},strcat('*',img_type,'*.tif')));
            roi_img_norm_files = [roi_img_norm_files; idx_4_files];
            
            idx_3_files = dir(fullfile(data_folder,'IHC','Results',m_names{j},strcat('*',img_type,'*mask_accepted.mat')));
            roi_mask_files = [roi_mask_files; idx_3_files];
            
        end

        % faster but wrong order
    %     idx_files = dir(fullfile(data_folder,'IHC','Results','**',strcat('*',img_type,'*results.mat')));
    %     result_files = [result_files; idx_files];
    end
    
    %% Sort result files
    [~,s_idxs1] = sort({result_files.name});
    [~,s_idxs2] = sort({roi_img_files.name});
    [~,s_idxs4] = sort({roi_img_norm_files.name});
    [~,s_idxs3] = sort({roi_mask_files.name});

    result_files = result_files(s_idxs1);
    roi_img_files = roi_img_files(s_idxs2);
    roi_img_norm_files = roi_img_norm_files(s_idxs4);
    roi_mask_files = roi_mask_files(s_idxs3);

%     {result_files.name}'
%     {roi_img_files.name}'
%     {roi_mask_files.name}'

    %% Find mouse names per condition (with data present)
%     mouse_names_present = {};
%     for i = 1:length(result_files)
%         mouse_names_present{i} = result_files(i).name(1:9);
%     end
%     
%     mouse_names_present = unique(mouse_names_present);
%     mouse_present_idxs = contains(mouse_names,mouse_names_present);
%     
    % global idxs
%     sham_gidxs = find(mouse_cond_idxs==1 & mouse_present_idxs);
%     gamma_gidxs = find(mouse_cond_idxs==2 & mouse_present_idxs);
%     theta_gidxs = find(mouse_cond_idxs==3 & mouse_present_idxs);
%     ltd_gidxs = find(mouse_cond_idxs==4 & mouse_present_idxs);
%     
%     % present idxs
%     [~,sham_pidxs,~] = intersect(find(mouse_cond_idxs==1),sham_gidxs);
%     [~,gamma_pidxs,~] = intersect(find(mouse_cond_idxs==2),gamma_gidxs);
%     [~,theta_pidxs,~] = intersect(find(mouse_cond_idxs==3),theta_gidxs);
%     [~,ltd_pidxs,~] = intersect(find(mouse_cond_idxs==4),ltd_gidxs);
    
%     
%     sham_names = mouse_names(sham_gidxs)';
%     gamma_names = mouse_names(gamma_gidxs)';
%     theta_names = mouse_names(theta_gidxs)';
%     ltd_names = mouse_names(ltd_gidxs)';
%     
%     condition_mouse_names = {sham_names,gamma_names,theta_names,ltd_names};

    %% Conditions
    cond_names = {'Sham','40 Hz','8 Hz','LTD'};
    cond_no = length(cond_names);
    variable_names_t = {'Sham','Gamma (40 Hz)','Theta (8 Hz)','LTD (1 Hz)'};
    max_n = 7;

    [roi_names,roi_fnames,roi_no] = get_roi_list();
    
%     if strcmp(img_type,'ki67') || strcmp(img_type,'dcx') || strcmp(img_type,'sox2')
%         roi_idxs = [1:2]; % DG only
%     else
%         roi_idxs = 1:roi_no;
%     end
    
    roi_idxs = 1:roi_no;
    
    %% Results to summarise
    if strcmp(img_type,'dcx')
        check_count = 0;
    else
        check_count = 1;
    end
    
    if strcmp(img_type,'cfos')
        check_perc_positive = 1;
    else
        check_perc_positive = 0;
    end
    
    if strcmp(img_type,'iba1')
        check_size = 1;
    else
        check_size = 0;
    end
    

    %% Load results
    results = {};
    roi_density = nan(roi_no,mouse_no);
    roi_count = nan(roi_no,mouse_no);
    roi_cfos_ratio = nan(roi_no,mouse_no);
    roi_size = cell(roi_no,mouse_no);

    for roi_idx = roi_idxs % 1:roi_no
        roi_fname = roi_fnames{roi_idx};
        file_idxs = find(contains({result_files.name}',roi_fname));

        for i = 1:length(file_idxs)
            file_idx = file_idxs(i);

            results{roi_idx,i} = load(fullfile(result_files(file_idx).folder,result_files(file_idx).name)).results;
            roi_density(roi_idx,i) = results{roi_idx,i}.density;
            roi_count(roi_idx,i) = results{roi_idx,i}.particle_no;
            
            if check_perc_positive
                roi_cfos_ratio(roi_idx,i) = results{roi_idx,i}.pos_particle_ratio;
            end
            
            if check_size
                roi_size{roi_idx,i} = results{roi_idx,i}.particle_radii;
            end
        end
    end

    %% Results summary: density
    results_density_all = plot_results(roi_density,'density',...
        mouse_cond_idxs,img_type,cohort_results_folder,roi_idxs,1);
    
%     % sham
%     sham_density = roi_density(roi_idxs,mouse_cond_idxs==1);
%     sham_n = size(sham_density,2);
% 
%     % 40 Hz
%     gamma_density = roi_density(roi_idxs,mouse_cond_idxs==2);
%     gamma_n = size(gamma_density,2);
% 
%     % 8 Hz
%     theta_density = roi_density(roi_idxs,mouse_cond_idxs==3);
%     theta_n = size(theta_density,2);
% 
%     % LTD
%     ltd_density = roi_density(roi_idxs,mouse_cond_idxs==4);
%     ltd_n = size(ltd_density,2);
% 
% 
%     [nanmean(sham_density,2),nanmean(gamma_density,2),nanmean(theta_density,2),nanmean(ltd_density,2)]
%     % [nanstd(sham_density,'',2),nanstd(gamma_density,'',2),nanstd(theta_density,'',2),nanstd(ltd_density,'',2)]
% 
%     max_n = 7;
%     normalise_to_sham = 0;
%     normalise_to_control = 0;
%     
% %     if normalise_to_control
% %         control_idx = 9;
% %         
% %         roi_results = nan(max_n,4);
% %         roi_results(1:sham_n,1) = sham_density(control_idx,:);
% %         roi_results(1:gamma_n,2) = gamma_density(control_idx,:);
% %         roi_results(1:theta_n,3) = theta_density(control_idx,:);
% %         roi_results(1:ltd_n,4) = ltd_density(control_idx,:);
% %         
% %         control_results = roi_results;
% %     end
%     
%     results_density_all = {};
% 
%     for roi_idx = roi_idxs
%     %     figure,boxplot([sham_density(roi_idx,:);gamma_density(roi_idx,:)]',...
%     %         'Labels',conds(1:2));
%         roi_name = roi_names{roi_idx};
% 
%         roi_results_density = nan(max_n,4);
%         roi_results_density(1:sham_n,1) = sham_density(roi_idx,:);
%         roi_results_density(1:gamma_n,2) = gamma_density(roi_idx,:);
%         roi_results_density(1:theta_n,3) = theta_density(roi_idx,:);
%         roi_results_density(1:ltd_n,4) = ltd_density(roi_idx,:);
%         
%         
%         if normalise_to_sham
%             mean_results = mean(roi_results_density,'omitnan');
%             results_to_plot = roi_results_density./mean_results(1)*100;
%         elseif normalise_to_control
%             results_to_plot = roi_results_density./control_results*100;
%         else
%             results_to_plot = roi_results_density;
%         end
% 
%         figure,hold on
%         x = [ones(max_n,1),2*ones(max_n,1),3*ones(max_n,1),4*ones(max_n,1)];
%         b = boxchart(results_to_plot);
%         b.BoxFaceColor = [0,0,0];
%         swarmchart(x,results_to_plot,'k','filled','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
%         xticklabels(cond_names)
% %         boxplot(roi_results,'Labels',conds);
% 
%         title(roi_names{roi_idx});
%         ylabel('Area covered (%)');
% 
%         fig_name = sprintf('%s_density_%d_roi_%s.tif',img_type,roi_idx,roi_fnames{roi_idx});
%         saveas(gcf,fullfile(cohort_results_folder,fig_name));
%         if close_figs; close(gcf); end
%         
%         file_name = sprintf('%s_density_%d_roi_%s_results',img_type,roi_idx,roi_fnames{roi_idx});
%         save(fullfile(stats_folder,strcat(file_name,'.mat')),'roi_results_density');
%         results_density_all{roi_idx} = roi_results_density;
%         
%         
%         % save results as an excel table
%         row_names = {};
%         for i = 1:max_n; row_names{i} = [roi_name,' ',num2str(i)]; end
% 
%         roi_results_density_T = array2table(roi_results_density,'VariableNames',...
%             cond_names,'RowNames',row_names);
%         
%         % add mean and std
%         extra_T = array2table([mean(roi_results_density,'omitnan');std(roi_results_density,'omitnan')],...
%             'VariableNames',cond_names,'RowNames',{'Mean','Std'});
%         roi_results_density_T = [roi_results_density_T; extra_T];
%         
% %         roi_results_density_T = array2table(roi_results_density,'VariableNames',cond_names);
%         table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
%         writetable(roi_results_density_T,table_name,'WriteRowNames',true);
%         
%     end
    
    %% Results summary: cell count
%     results_count_all = {};

    if check_count
        results_count_all = plot_results(roi_count,'count',...
            mouse_cond_idxs,img_type,cohort_results_folder,roi_idxs,1);

    
%         sham_particle_no = roi_particle_no(roi_idxs,mouse_cond_idxs==1);
%         gamma_particle_no = roi_particle_no(roi_idxs,mouse_cond_idxs==2);
%         theta_particle_no = roi_particle_no(roi_idxs,mouse_cond_idxs==3);
%         ltd_particle_no = roi_particle_no(roi_idxs,mouse_cond_idxs==4);
% 
%         [nanmean(sham_particle_no,2),nanmean(gamma_particle_no,2),nanmean(theta_particle_no,2),nanmean(ltd_particle_no,2)]
%         % [nanstd(sham_density,'',2),nanstd(gamma_density,'',2),nanstd(theta_density,'',2),nanstd(ltd_density,'',2)]
% 
%         for roi_idx = roi_idxs
%             roi_name = roi_names{roi_idx};
% 
%             roi_results_count = nan(max_n,4);
%             roi_results_count(1:sham_n,1) = sham_particle_no(roi_idx,:);
%             roi_results_count(1:gamma_n,2) = gamma_particle_no(roi_idx,:);
%             roi_results_count(1:theta_n,3) = theta_particle_no(roi_idx,:);
%             roi_results_count(1:ltd_n,4) = ltd_particle_no(roi_idx,:);
% 
%             if normalise_to_sham
%                 mean_results = mean(roi_results_count,'omitnan');
%                 results_to_plot = roi_results_count./mean_results(1)*100;
%             elseif normalise_to_control
%                 results_to_plot = roi_results_count./control_results*100;
%             else
%                 results_to_plot = roi_results_count;
%             end
% 
%             figure,hold on
%             x = [ones(max_n,1),2*ones(max_n,1),3*ones(max_n,1),4*ones(max_n,1)];
%             b = boxchart(results_to_plot);
%             b.BoxFaceColor = [0,0,0];
%             swarmchart(x,results_to_plot,'k','filled','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
%             xticklabels(cond_names)
% 
%             title(roi_names{roi_idx});
%             ylabel('Cell count');
% 
%             fig_name = sprintf('%s_count_%d_roi_%s.tif',img_type,roi_idx,roi_fnames{roi_idx});
%             saveas(gcf,fullfile(cohort_results_folder,fig_name));
%             if close_figs; close(gcf); end
% 
%             file_name = sprintf('%s_count_%d_roi_%s_results',img_type,roi_idx,roi_fnames{roi_idx});
%             save(fullfile(stats_folder,strcat(file_name,'.mat')),'roi_results_count');
% 
%             results_count_all{roi_idx} = roi_results_count;
%             
%             % save results as an excel table
% %             roi_results_count_T = array2table(roi_results_count,'VariableNames',cond_names); % cond_names
%             
%             row_names = {};
%             for i = 1:max_n; row_names{i} = [roi_name,' ',num2str(i)]; end
% 
%             roi_results_count_T = array2table(roi_results_count,'VariableNames',...
%                 cond_names,'RowNames',row_names);
% 
%             % add mean and std
%             extra_T = array2table([mean(roi_results_count,'omitnan');std(roi_results_count,'omitnan')],...
%                 'VariableNames',cond_names,'RowNames',{'Mean','Std'});
%             roi_results_count_T = [roi_results_count_T; extra_T];
%             
%             table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
%             writetable(roi_results_count_T,table_name,'WriteRowNames',true);

%         end
    end
    
    %% Results summary: % of positive cells (cfos only)
%     results_ratio_all = {};
    
    if check_perc_positive
        
        results_ratio_all = plot_results(roi_cfos_ratio,'cfos_ratio',...
            mouse_cond_idxs,img_type,cohort_results_folder,roi_idxs,1);

%         sham_pos_ratio = roi_pos_particle_ratio(roi_idxs,mouse_cond_idxs==1);
%         gamma_pos_ratio = roi_pos_particle_ratio(roi_idxs,mouse_cond_idxs==2);
%         theta_pos_ratio = roi_pos_particle_ratio(roi_idxs,mouse_cond_idxs==3);
%         ltd_pos_ratio = roi_pos_particle_ratio(roi_idxs,mouse_cond_idxs==4);
% 
%         [nanmean(sham_pos_ratio,2),nanmean(gamma_pos_ratio,2),nanmean(theta_pos_ratio,2),nanmean(ltd_pos_ratio,2)];
%         % [nanstd(sham_density,'',2),nanstd(gamma_density,'',2),nanstd(theta_density,'',2),nanstd(ltd_density,'',2)]
% 
%         for roi_idx = roi_idxs
%             roi_name = roi_names{roi_idx};
% 
%             roi_results_ratio = nan(max_n,4);
%             roi_results_ratio(1:sham_n,1) = sham_pos_ratio(roi_idx,:);
%             roi_results_ratio(1:gamma_n,2) = gamma_pos_ratio(roi_idx,:);
%             roi_results_ratio(1:theta_n,3) = theta_pos_ratio(roi_idx,:);
%             roi_results_ratio(1:ltd_n,4) = ltd_pos_ratio(roi_idx,:);
% 
%             if normalise_to_sham
%                 mean_results = mean(roi_results_ratio,'omitnan');
%                 results_to_plot = (roi_results_ratio./mean_results(1)-1)*100;
%             elseif normalise_to_control
%                 results_to_plot = roi_results_ratio./control_results*100;
%             else
%                 results_to_plot = roi_results_ratio;
%             end
% 
%             figure,hold on
%             x = [ones(max_n,1),2*ones(max_n,1),3*ones(max_n,1),4*ones(max_n,1)];
%             b = boxchart(results_to_plot);
%             b.BoxFaceColor = [0,0,0];
%             swarmchart(x,results_to_plot,'k','filled','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
%             xticklabels(cond_names)
% 
%             title(roi_names{roi_idx});
%             ylabel('% cfos positive');
% 
%             fig_name = sprintf('%s_pos_ratio_%d_roi_%s.tif',img_type,roi_idx,roi_fnames{roi_idx});
%             saveas(gcf,fullfile(cohort_results_folder,fig_name));
%             if close_figs; close(gcf); end
%             
%             results_ratio_all{roi_idx} = roi_results_ratio;
% 
%             file_name = sprintf('%s_pos_ratio_%d_roi_%s_results',img_type,roi_idx,roi_fnames{roi_idx});
%             save(fullfile(stats_folder,strcat(file_name,'.mat')),'roi_results_ratio');
% 
%             % save results as an excel table
% %             roi_results_ratio_T = array2table(roi_results_ratio,'VariableNames',cond_names); % cond_names
%             row_names = {};
%             for i = 1:max_n; row_names{i} = [roi_name,' ',num2str(i)]; end
% 
%             roi_results_ratio_T = array2table(roi_results_ratio,'VariableNames',...
%                 cond_names,'RowNames',row_names);
% 
%             % add mean and std
%             extra_T = array2table([mean(roi_results_ratio,'omitnan');std(roi_results_ratio,'omitnan')],...
%                 'VariableNames',cond_names,'RowNames',{'Mean','Std'});
%             roi_results_ratio_T = [roi_results_ratio_T; extra_T];
%             
%             table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
%             writetable(roi_results_ratio_T,table_name,'WriteRowNames',true);
% 
%         end
    end
    
    %% Results summary: size of microglia
%     results_radius_all = {};
    
    if check_size
        
        results_size_all = plot_results(roi_radius,'size',...
            mouse_cond_idxs,img_type,cohort_results_folder,roi_idxs,1);

%         sham_radius = roi_radius(roi_idxs,mouse_cond_idxs==1);
%         gamma_radius = roi_radius(roi_idxs,mouse_cond_idxs==2);
%         theta_radius = roi_radius(roi_idxs,mouse_cond_idxs==3);
%         ltd_radius = roi_radius(roi_idxs,mouse_cond_idxs==4);
% 
%         for roi_idx = roi_idxs
%             roi_name = roi_names{roi_idx};
% 
% %             roi_results_ratio = nan(max_n,4);
%             roi_results_diameter_ = {};
%             roi_results_diameter_{1} = [sham_radius{roi_idx,:}]*2;
%             roi_results_diameter_{2} = [gamma_radius{roi_idx,:}]*2;
%             roi_results_diameter_{3} = [theta_radius{roi_idx,:}]*2;
%             roi_results_diameter_{4} = [ltd_radius{roi_idx,:}]*2;
%             
%             max_n_ = max(cellfun(@length,roi_results_ratio_));
%             roi_results_diameter = nan(max_n_,4);
%             roi_results_ratio(1:length(roi_results_ratio_{1}),1) = roi_results_ratio_{1};
%             roi_results_diameter(1:length(roi_results_ratio_{2}),2) = roi_results_ratio_{2};
%             roi_results_diameter(1:length(roi_results_ratio_{3}),3) = roi_results_ratio_{3};
%             roi_results_diameter(1:length(roi_results_ratio_{4}),4) = roi_results_ratio_{4};
% 
% 
% %             if normalise_to_sham
% %                 mean_results = mean(roi_results_ratio,'omitnan');
% %                 results_to_plot = (roi_results_ratio./mean_results(1)-1)*100;
% %             elseif normalise_to_control
% %                 results_to_plot = roi_results_ratio./control_results*100;
% %             else
%             results_to_plot = roi_results_diameter; % plot diameter
% %             end
% 
%             figure,hold on
%             x = [ones(max_n_,1),2*ones(max_n_,1),3*ones(max_n_,1),4*ones(max_n_,1)];
%             b = boxchart(results_to_plot);
%             b.BoxFaceColor = [0,0,0]; b.MarkerColor = [0,0,0];
%             swarmchart(x,results_to_plot,'k','filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.1)
%             xticklabels(cond_names)
% 
%             title(roi_names{roi_idx});
%             ylabel('Cell diameter (μm)');
% 
%             fig_name = sprintf('%s_diameter_%d_roi_%s.tif',img_type,roi_idx,roi_fnames{roi_idx});
%             saveas(gcf,fullfile(cohort_results_folder,fig_name));
%             if close_figs; close(gcf); end
%             
%             results_radius_all{roi_idx} = roi_results_ratio;
% 
%             file_name = sprintf('%s_diameter_%d_roi_%s_results',img_type,roi_idx,roi_fnames{roi_idx});
%             save(fullfile(stats_folder,strcat(file_name,'.mat')),'roi_results_ratio');
% 
%             % save results as an excel table
% %             row_names = {};
% %             for i = 1:max_n; row_names{i} = [roi_name,' ',num2str(i)]; end
% % 
% %             roi_results_ratio_T = array2table(roi_results_ratio,'VariableNames',...
% %                 cond_names,'RowNames',row_names);
% % 
% %             % add mean and std
% %             extra_T = array2table([mean(roi_results_ratio,'omitnan');std(roi_results_ratio,'omitnan')],...
% %                 'VariableNames',cond_names,'RowNames',{'Mean','Std'});
% %             roi_results_ratio_T = [roi_results_ratio_T; extra_T];
% %             
% %             table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
% %             writetable(roi_results_ratio_T,table_name,'WriteRowNames',true);
% 
%         end
    end
    
    %% Statistics
    % not all data normally distributed; use wilcoxon test instead of ttest
    
    p_density = nan(length(roi_idxs),3);
    h_density = nan(length(roi_idxs),3);
    
    p_count = nan(length(roi_idxs),3);
    h_count = nan(length(roi_idxs),3);
    
    p_cfos_ratio = nan(length(roi_idxs),3);
    h_cfos_ratio = nan(length(roi_idxs),3);
    
    p_size = nan(length(roi_idxs),3);
    h_size = nan(length(roi_idxs),3);
        
    
    for roi_idx = roi_idxs
        
        % DENSITY
        results_density = results_density_all{roi_idx};
    %     [h,p] = kstest(results_density); % h = 0 if normally distributed

        % compare all stims to sham
        [p_density(roi_idx,1),h_density(roi_idx,1)] = ranksum(results_density(:,1),results_density(:,2));
        [p_density(roi_idx,2),h_density(roi_idx,2)] = ranksum(results_density(:,1),results_density(:,3));
        [p_density(roi_idx,3),h_density(roi_idx,3)] = ranksum(results_density(:,1),results_density(:,4));
        
        
        % CELL COUNT
        if check_count
            results_count = results_count_all{roi_idx};

            [p_count(roi_idx,1),h_count(roi_idx,1)] = ranksum(results_count(:,1),results_count(:,2));
            [p_count(roi_idx,2),h_count(roi_idx,2)] = ranksum(results_count(:,1),results_count(:,3));
            [p_count(roi_idx,3),h_count(roi_idx,3)] = ranksum(results_count(:,1),results_count(:,4));     
        end
        
        % CFOS CFOS POSITIVE RATIO
        if check_perc_positive
            results_ratio = results_ratio_all{roi_idx};
            
            [p_cfos_ratio(roi_idx,1),h_cfos_ratio(roi_idx,1)] = ranksum(results_ratio(:,1),results_ratio(:,2));
            [p_cfos_ratio(roi_idx,2),h_cfos_ratio(roi_idx,2)] = ranksum(results_ratio(:,1),results_ratio(:,3));
            [p_cfos_ratio(roi_idx,3),h_cfos_ratio(roi_idx,3)] = ranksum(results_ratio(:,1),results_ratio(:,4));     
        end
        
        % IBA1 SIZE
        if check_size
            results_size = results_size_all{roi_idx};
            
            [p_size(roi_idx,1),h_size(roi_idx,1)] = ranksum(results_size(:,1),results_size(:,2));
            [p_size(roi_idx,2),h_size(roi_idx,2)] = ranksum(results_size(:,1),results_size(:,3));
            [p_size(roi_idx,3),h_size(roi_idx,3)] = ranksum(results_size(:,1),results_size(:,4));     
        end
    end
    
    
    %% Save stats
    var_names = {'Gamma vs Sham','Theta vs Sham','LTD vs Sham'};

    p = p_density; h = h_density;
    file_name = sprintf('%s_density_stats',img_type);
    save(fullfile(stats_folder,strcat(file_name,'.mat')),'p','h');
    
    stats_T = array2table(p,'VariableNames',var_names,'RowNames',roi_names(roi_idxs));
    table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
    writetable(stats_T,table_name,'WriteRowNames',true);

    if check_count
        p = p_count; h = h_count;
        file_name = sprintf('%s_count_stats',img_type);
        save(fullfile(stats_folder,strcat(file_name,'.mat')),'p','h');

        stats_T = array2table(p,'VariableNames',var_names,'RowNames',roi_names(roi_idxs));
        table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
        writetable(stats_T,table_name,'WriteRowNames',true);
    end
    
    if check_perc_positive
        p = p_cfos_ratio; h = h_cfos_ratio;
        file_name = sprintf('%s_cfos_ratio_stats',img_type);
        save(fullfile(stats_folder,strcat(file_name,'.mat')),'p','h');

        stats_T = array2table(p,'VariableNames',var_names,'RowNames',roi_names(roi_idxs));
        table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
        writetable(stats_T,table_name,'WriteRowNames',true);
    end
    
    if check_size
        p = p_size; h = h_size;
        file_name = sprintf('%s_size_stats',img_type);
        save(fullfile(stats_folder,strcat(file_name,'.mat')),'p','h');

        stats_T = array2table(p,'VariableNames',var_names,'RowNames',roi_names(roi_idxs));
        table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
        writetable(stats_T,table_name,'WriteRowNames',true);

    end
    
    %% Plot DAB images for comparison
    [~,~,~,~,correct_brightness] = get_antibody_threshold(img_type);
    
    for roi_idx = roi_idxs
        
        roi_fname = roi_fnames{roi_idx};
        img_idxs = find(contains({roi_img_files.name}',roi_fname));
        
        roi_results_density = results_density_all{roi_idx};
        if isempty(results_count_all)
            roi_results_count = nan(size(roi_results_density));
        else
            roi_results_count = results_count_all{roi_idx};
        end
        
 
        for j = 1:cond_no
            fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
            fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
            
            if correct_brightness
                fig3 = figure('units','normalized','outerposition',[0 0 1 1]);
            end
            
            cond_j_names = condition_mouse_names{j};
            cond_img_idxs = img_idxs((contains({roi_img_files(img_idxs).name}',cond_j_names)));
%             {roi_img_files(cond_img_idxs).name}'

            % remove nan values
            roi_results_density_j = roi_results_density(:,j);
            roi_results_count_j = roi_results_count(:,j);
            
            missing_idxs = isnan(roi_results_density_j);
            if length(cond_img_idxs) == sum(~missing_idxs)
                roi_results_density_j(missing_idxs) = [];
                roi_results_count_j(missing_idxs) = [];
            end
            
            for i = 1:length(cond_img_idxs)
                % dab images
                [~,roi_img_dab] = load_deconvolved_images(fullfile(roi_img_files(cond_img_idxs(i)).folder,roi_img_files(cond_img_idxs(i)).name));
    %             roi_img = load(fullfile(roi_img_files(img_idxs(i)).folder,roi_img_files(img_idxs(i)).name)).roi_image;
                set(0,'CurrentFigure',fig1)
                subplot(round(max_n/2),2,i),imshow(roi_img_dab)
                colormap(dab_colormap)
                title(roi_img_files(cond_img_idxs(i)).name,'Interpreter','none')
                
                % dab images + antibody masks
                dab_roi_mask = load(fullfile(roi_mask_files(cond_img_idxs(i)).folder,roi_mask_files(cond_img_idxs(i)).name)).dab_roi_mask;
                roi_img_dab_mask = labeloverlay(roi_img_dab,dab_roi_mask,...
                    'Colormap',[0,0,1],'Transparency',0.2);
                
                set(0,'CurrentFigure',fig2)
                subplot(round(max_n/2),2,i),imshow(roi_img_dab_mask)
                title_txt = {roi_mask_files(cond_img_idxs(i)).name,...
                    sprintf('Cell count = %d; Density = %1.2f %%',roi_results_count_j(i),roi_results_density_j(i))};
                title(title_txt,'Interpreter','none')
                
                 % norm dab images
                if correct_brightness
                    [~,roi_img_norm_dab] = load_deconvolved_images(fullfile(roi_img_norm_files(cond_img_idxs(i)).folder,roi_img_norm_files(cond_img_idxs(i)).name));
                    set(0,'CurrentFigure',fig3)
                    subplot(round(max_n/2),2,i),imshow(roi_img_norm_dab)
                    colormap(dab_colormap)
                    title(roi_img_norm_files(cond_img_idxs(i)).name,'Interpreter','none')
                end
            end
            
            fig_1_name = sprintf('roi_images_%s_%d_%s_%d_%s.tif',img_type,roi_idx,roi_fnames{roi_idx},j,cond_names{j});
            saveas(fig1,fullfile(comparison_folder,fig_1_name));
            
            fig_2_name = sprintf('roi_images_%s_%d_%s_%d_%s_masks.tif',img_type,roi_idx,roi_fnames{roi_idx},j,cond_names{j});
            saveas(fig2,fullfile(comparison_folder,fig_2_name));
           
            if correct_brightness
                fig_3_name = sprintf('roi_images_%s_%d_%s_%d_%s_0_norm.tif',img_type,roi_idx,roi_fnames{roi_idx},j,cond_names{j});
                saveas(fig3,fullfile(comparison_folder,fig_3_name));
            end
            
            if close_figs
                close(fig1);
                close(fig2);
                if correct_brightness; close(fig3); end
            end

        end
        
    end
    
end
