function [results_all] = plot_results(quantity_to_plot,quantity_name,...
                        mouse_cond_idxs,img_type,results_folder,roi_idxs,save_results)
    %% Plot and save results
    % @author: pdzialecka
    
    %%
    if size(quantity_to_plot,1) == 6
        [roi_names,roi_fnames,roi_no] = get_roi_list_IF();
    else
        [roi_names,roi_fnames,roi_no] = get_roi_list();
    end
    
    cond_names = {'Sham','Delta','Theta','Gamma'}; % {'Sham','LTD','Theta','Gamma'};
    close_figs = 1;
    fontsize = 20;
    
    %%
    if ~exist('roi_idxs','var')
        roi_idxs = 1:roi_no;
    end
    
    if ~exist('save_results','var')
        save_results = 1;
    end
    
    %% Stats folder
    stats_folder = fullfile(results_folder,'Stats');
    if ~exist(stats_folder)
        mkdir(stats_folder)
    end  
    
    %%
    % IHC
    if strcmp(quantity_name,'density')
        ylabel_ = 'Area (%)'; %'Area covered (%)';
        y_round = 0.5;
        if strcmp(img_type,'sox2')
            y_round = 1;
        end
    elseif strcmp(quantity_name,'count')
        ylabel_ = 'Cell count';
        y_round = 10;
    elseif strcmp(quantity_name,'cfos_ratio')
        ylabel_ = '% of cfos positive cells';
        y_round = 0.25;
    elseif strcmp(quantity_name,'size')
        ylabel_ = 'Cell diameter (μm)';
        y_round = 10;
        
    % IF
    elseif strcmp(quantity_name,'microglia_ratio') || strcmp(quantity_name,'microglia_area_ratio')
        ylabel_ = '% of Aβ+ microglia';
        y_round = 10;
    elseif strcmp(quantity_name,'ab_ratio') || strcmp(quantity_name,'ab_area_ratio')
        ylabel_ = '% of microglia+ Aβ';
        y_round = 10;
    elseif strcmp(quantity_name,'microglia_per_ab')
        ylabel_ = 'Microglia no per Aβ';
        y_round = 1;
        
    % Behaviour
    elseif contains(quantity_name,'DI')
        ylabel_ = 'DI'; % 'Discrimination index';
        y_round = 1;
        roi_names = {''};
        roi_fnames = {''};
    elseif contains(quantity_name,'Time')
        ylabel_ = 'Exploration time (s)';
        y_round = 10;
        roi_names = {''};
        roi_fnames = {''};
    elseif contains(quantity_name,'index')
        ylabel_ = strcat(quantity_name(1:3),' (%)');
        y_round = 100;
        if contains(quantity_name,'SAR')
            y_round = 10;
        end
        roi_names = {''};
        roi_fnames = {''};
    end
    
    %%
    if iscell(quantity_to_plot)
        y_max = max(cellfun(@max,quantity_to_plot),[],[1,2],'omitnan')*2;
    else
        y_max = max(quantity_to_plot,[],[1,2],'omitnan');
    end
    
    y_max = ceil(y_max/y_round)*y_round; % round
    ylims_ = [0,y_max];
    
    %% Extract results per group
    % sham
    sham_results = quantity_to_plot(roi_idxs,mouse_cond_idxs==1);
    sham_n = size(sham_results,2);

    % LTD
    ltd_results = quantity_to_plot(roi_idxs,mouse_cond_idxs==4);
    ltd_n = size(ltd_results,2);

    % 8 Hz
    theta_results = quantity_to_plot(roi_idxs,mouse_cond_idxs==3);
    theta_n = size(theta_results,2);
    
    % 40 Hz
    gamma_results = quantity_to_plot(roi_idxs,mouse_cond_idxs==2);
    gamma_n = size(gamma_results,2);

%     [nanmean(sham_density,2),nanmean(gamma_density,2),nanmean(theta_density,2),nanmean(ltd_density,2)]
%     [nanstd(sham_density,'',2),nanstd(gamma_density,'',2),nanstd(theta_density,'',2),nanstd(ltd_density,'',2)]


    max_n = 7;
    normalise_to_sham = 0;
    normalise_to_control = 0;
    
%     if normalise_to_control
%         control_idx = 9;
%         
%         roi_results = nan(max_n,4);
%         roi_results(1:sham_n,1) = sham_density(control_idx,:);
%         roi_results(1:gamma_n,2) = gamma_density(control_idx,:);
%         roi_results(1:theta_n,3) = theta_density(control_idx,:);
%         roi_results(1:ltd_n,4) = ltd_density(control_idx,:);
%         
%         control_results = roi_results;
%     end
    
    results_all = {};

    %%
    for roi_idx = roi_idxs
        %%
        roi_name = roi_names{roi_idx};

        
        if iscell(quantity_to_plot)
            roi_results_ = {};
            roi_results_{1} = [sham_results{roi_idx,:}];
            roi_results_{2} = [ltd_results{roi_idx,:}];
            roi_results_{3} = [theta_results{roi_idx,:}];
            roi_results_{4} = [gamma_results{roi_idx,:}];


            max_n = max(cellfun(@length,roi_results_));
            roi_results = nan(max_n,4);
            roi_results(1:length(roi_results_{1}),1) = roi_results_{1};
            roi_results(1:length(roi_results_{2}),2) = roi_results_{2};
            roi_results(1:length(roi_results_{3}),3) = roi_results_{3};
            roi_results(1:length(roi_results_{4}),4) = roi_results_{4};
        
        else
            roi_results = nan(max_n,4);
            roi_results(1:sham_n,1) = sham_results(roi_idx,:);
            roi_results(1:ltd_n,2) = ltd_results(roi_idx,:);
            roi_results(1:theta_n,3) = theta_results(roi_idx,:);
            roi_results(1:gamma_n,4) = gamma_results(roi_idx,:);
        end
        
        
        if normalise_to_sham
            mean_results = mean(roi_results,'omitnan');
            results_to_plot = roi_results./mean_results(1)*100;
        elseif normalise_to_control
            results_to_plot = roi_results./control_results*100;
        else
            results_to_plot = roi_results;
        end
        
        
        results_all{roi_idx} = roi_results;

        %% Plot results
        figure,hold on
        x = [ones(max_n,1),2*ones(max_n,1),3*ones(max_n,1),4*ones(max_n,1)];
        b = boxchart(results_to_plot);
        b.BoxFaceColor = [0,0,0]; b.MarkerColor = [0,0,0];
        
        if max_n < 30
            s = swarmchart(x,results_to_plot,'k','filled','MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8);
        else
            s = swarmchart(x,results_to_plot,'k','filled','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0.1);
        end
        s(1).XJitterWidth = 0.025*length(cond_names); s(2).XJitterWidth = 0.025*length(cond_names);
        s(3).XJitterWidth = 0.025*length(cond_names); s(4).XJitterWidth = 0.025*length(cond_names);

        xticklabels(cond_names)
%         boxplot(roi_results,'Labels',conds);

        title(roi_names{roi_idx});
        ylabel(ylabel_); ylim(ylims_)
        set(gca,'box','off','Fontsize',fontsize)
        
        if contains(quantity_name,'DI')
            yline(0.5,'--k','LineWidth',1);
        elseif contains(quantity_name,'SAP')
            yline(6/27*100,'--k','LineWidth',1);
        end

        
        if save_results
            roi_str = sprintf('_%d_roi_%s',roi_idx,roi_fnames{roi_idx});
            if isempty(roi_fnames{roi_idx})
                roi_str = '';
            end
            
            fig_name = sprintf('%s_%s%s',img_type,quantity_name,roi_str); % roi_idx,roi_fnames{roi_idx});
            saveas(gcf,fullfile(results_folder,strcat(fig_name,'.tif')));
            saveas(gcf,fullfile(results_folder,strcat(fig_name,'.fig')));

            if close_figs; close(gcf); end

            file_name = sprintf('%s_%s%s_results',img_type,quantity_name,roi_str); % roi_idx,roi_fnames{roi_idx});
            save(fullfile(stats_folder,strcat(file_name,'.mat')),'roi_results');
        end
        
        %% Save results as an excel table
        row_names = {};
        for i = 1:max_n; row_names{i} = [roi_name,' ',num2str(i)]; end

        roi_results_T = array2table(roi_results,'VariableNames',...
            cond_names,'RowNames',row_names);
        
        % add mean and std
        extra_T = array2table([mean(roi_results,'omitnan');std(roi_results,'omitnan');sum(~isnan(roi_results))],...
            'VariableNames',cond_names,'RowNames',{'Mean','Std','n'});
        roi_results_T = [roi_results_T; extra_T];
        
        if save_results
            table_name = fullfile(stats_folder,strcat(file_name,'.xlsx'));
            writetable(roi_results_T,table_name,'WriteRowNames',true);
        end
        
    end
end
