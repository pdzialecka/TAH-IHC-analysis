function [] = summarise_results(base_folder,cohort_case,image_type)
    %% Create summary plot of results per image type
	% @author: pdzialecka

    %% Data folders
%     base_folder = 'C:\Users\Pat\Desktop\TAH';
%     cohort_case = 1; % 1 = cohort 1 (13 mo). 2 = cohorts 2-5 (6 mo)

    if cohort_case == 1
        cohort_idxs = [1];
    elseif cohort_case == 2
        cohort_idxs = [2,3]; % [2,3,4,5];
    end

    all_cohort_folders = dir(fullfile(base_folder,'Data','Cohort*'));
    data_folders = all_cohort_folders(cohort_idxs);

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
    cohort_results_folder = fullfile(cohort_results_folder,image_type);

    if ~exist(cohort_results_folder)
        mkdir(cohort_results_folder)
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
    mouse_no = length(mouse_ids);

    %% Find result files
    result_files = [];

    for i = 1:length(data_folders)
        data_folder = fullfile(data_folders(i).folder,data_folders(i).name);

        cohort_names = cohort_infos{i}.mouse_names;

        for j = 1:length(cohort_names)
            idx_files = dir(fullfile(data_folder,'IHC','Results',cohort_names{j},strcat('*',image_type,'*results.mat')));

            result_files = [result_files; idx_files];
        end

        % faster but wrong order
    %     idx_files = dir(fullfile(data_folder,'IHC','Results','**',strcat('*',image_type,'*results.mat')));
    %     result_files = [result_files; idx_files];
    end


    %% Load results
    [roi_names,roi_fnames,roi_no] = get_roi_list();

    results = {};
    roi_density = nan(roi_no,mouse_no);

    for roi_idx = 1:roi_no
        roi_fname = roi_fnames{roi_idx};
        file_idxs = find(contains({result_files.name}',roi_fname));

        for i = 1:length(file_idxs)
            file_idx = file_idxs(i);

            results{roi_idx,i} = load(fullfile(result_files(file_idx).folder,result_files(file_idx).name)).results;
            roi_density(roi_idx,i) = results{roi_idx,i}.density;
        end
    end

    %% Plot summary results
    conds = {'Sham','40 Hz','8 Hz','LTD'};

    % sham
    sham_density = roi_density(:,mouse_cond_idxs==1);
    sham_n = size(sham_density,2);

    % 40 Hz
    gamma_density = roi_density(:,mouse_cond_idxs==2);
    gamma_n = size(gamma_density,2);

    % 8 Hz
    theta_density = roi_density(:,mouse_cond_idxs==3);
    theta_n = size(theta_density,2);

    % LTD
    ltd_density = roi_density(:,mouse_cond_idxs==4);
    ltd_n = size(ltd_density,2);


    [nanmean(sham_density,2),nanmean(gamma_density,2),nanmean(theta_density,2),nanmean(ltd_density,2)]
    % [nanstd(sham_density,'',2),nanstd(gamma_density,'',2),nanstd(theta_density,'',2),nanstd(ltd_density,'',2)]


    max_group_n = 6;

    for roi_idx = 1:roi_no
    %     figure,boxplot([sham_density(roi_idx,:);gamma_density(roi_idx,:)]',...
    %         'Labels',conds(1:2));

        results_to_plot = nan(max_group_n,4);
        results_to_plot(1:sham_n,1) = sham_density(roi_idx,:);
        results_to_plot(1:gamma_n,2) = gamma_density(roi_idx,:);
        results_to_plot(1:theta_n,3) = theta_density(roi_idx,:);
        results_to_plot(1:ltd_n,4) = ltd_density(roi_idx,:);

        figure,boxplot(results_to_plot,'Labels',conds);

        title(roi_names{roi_idx});
        ylabel('Area covered (%)');

        fig_name = sprintf('%s_density_roi_%s.tif',image_type,roi_fnames{roi_idx});
        saveas(gcf,fullfile(cohort_results_folder,fig_name));
        close(gcf);
    end
    
end
