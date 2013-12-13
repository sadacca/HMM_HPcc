function star = Data(star)
%% function 'star = Data'
% Reads-in neural spike train data for N trials, computes the peri-stimulus
% histogram for each tastant, then computes the correlation between
% selected PC and behavioral palatability data at each time
% bin, selecting the best correlating PC

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Read in data
count = 1;

for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        resultfilename = getfield(star.resultfiles(count), 'name'); % get name
        load(resultfilename, 'par'); % load structure
        
        data_matrices{i, j} = par.observations;
        data_matrices_norm{i, j} = par.norm_observations;
        n_trials(i, j) = par.n_trials;
        n_symbols(i, j) = par.n_symbols;
        n_bins(i, j) = par.n_bins;
        count = count + 1;
        
        clear par
        
    end
    
end

%% Analyze neural data
for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        if star.select_trial_flag % use only selected trials (in which target state appears)
            
            data_to_use = data_matrices{i, j};
            n_trials(i, j) = star.n_trials(i, j);
            
            for n = 1:n_trials(i, j)
                
                temp(n, :) = data_to_use(star.new_trial_indices{i, j}(n), :);
                target_trial_spikes(i, j, n, :) = temp(n, :);
                
            end
            
            data_matrices{i, j} = squeeze(target_trial_spikes(i, j, :, :)); % select only targeted trials
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         temp_indices = [1 7 8 9];
        %         data_to_use_test = data_matrices{i, j};
        %             for n = 1:length(temp_indices)
        %                 target_trial_spikes_test(i, j, n, :) = data_to_use_test(temp_indices(n), :);
        %             end
        %         spikes_test = squeeze(target_trial_spikes_test(i, j, :, :))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for k = 1:star.max_n_symbols
            
            psth_data{i, j, k} = mean(data_matrices{i, j} == k); % average spikes over trials
            
        end
        
    end
    
end


%% Smooth neural data
for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        psth_data_double(i, j, :, :) = squeeze(cell2mat(psth_data(i, j, :))); % convert to double
        
        for k = 1:star.max_n_symbols
            
            psth_data_double(i, j, :, k) = smooth(squeeze(psth_data_double(i, j, :, k)), star.smooth_factor); % average (smooth) over time bins
            
        end
        
    end
    
end

%% Correlate PCA and palatability for time series

switch star.PCA_L2_toggle
    
    case 'PCA'
        
        for i = 1:star.n_datasets
            
            for bin_loop = star.starttime:star.endtime
                
                [a(:, :, bin_loop) b(:, :, bin_loop) latent] = princomp(squeeze(psth_data_double(i, :, bin_loop, 2:end))); % do the pca
                c(:, bin_loop) = latent ./ sum(latent);
                
            end
            
            coeff{i} = a;
            score{i} = b;
            eigenvalue{i} = c;
            
            for pc_loop = 1:3
                
                [r(:, pc_loop) p(:, pc_loop)] = corr(score{i}(:, pc_loop, :), star.palatability_vals, 'type', star.correlation_toggle); % calculate correlation
                
            end
            
            rho{i} = r;
            pvalue{i} = p;
            [dummy star.best_index(i)] = max(rho{i}(:, star.pc_to_show).^2);
            star.best_index(i) = min(star.best_index(i) + star.mag_perm_toggle, star.endtime);
            
            switch star.plot_indiv_datasets_toggle
                case('on')
                    
                    figure(star.fig + 100 + i), clf
                    
                    subplot(2, 2, 1)
                    plot(squeeze(score{i}(:, star.pc_to_show, :))', 'LineWidth', 2)
                    hold on
                    minmax = ylim;
                    line([star.best_index(i) star.best_index(i)], minmax, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                    axis([0 star.endtime minmax])
                    xlabel('Time bin', 'FontSize', 14)
                    ylabel('PC value', 'FontSize', 14)
                    title(['PC projection: ' star.model_toggle ' ' star.include_states_toggle], 'FontSize', 16)
                    
                    subplot(2, 2, 2)
                    plot(repmat(star.palatability_vals, 1, n_bins(1, 1))', 'LineWidth', 2)
                    minmax = ylim;
                    axis([0 star.endtime minmax])
                    xlabel('Time bin', 'FontSize', 14)
                    ylabel('Palatability value', 'FontSize', 14)
                    title('Palatability', 'FontSize', 16)
                    legend(star.tastant_names)
                    
                    subplot(2, 2, 3)
                    plot(eigenvalue{i}(star.pc_to_show, :), 'LineWidth', 2)
                    hold on
                    line([star.best_index(i) star.best_index(i)], [0 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                    axis([0 star.endtime 0 1])
                    xlabel('Time bin', 'FontSize', 14)
                    ylabel(texlabel('R^{2}'), 'FontSize', 14)
                    title(['Var. explained PC' num2str(star.pc_to_show)], 'FontSize', 16)
                    
                    subplot(2, 2, 4)
                    plot(rho{i}(:, star.pc_to_show), 'LineWidth', 2)
                    hold on
                    line([star.best_index(i) star.best_index(i)], [-1 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                    axis([0 star.endtime -1 1])
                    xlabel('Time bin', 'FontSize', 14)
                    ylabel(texlabel('r'), 'FontSize', 14)
                    title(['Corr. PC' num2str(star.pc_to_show) ' and Pal.'], 'FontSize', 16)
                    
                    addtitle('Data', 24, [0.5 0.1]);
                    
            end
            
            clear a b c r p
            
        end
        
        switch star.plot_indiv_datasets_toggle
            
            case('off')
                
                mean_score = zeros(size(squeeze(score{1}(:, 1, :))'));
                mean_eigenvalue = zeros(size(eigenvalue{1}(1, :)));
                mean_rho = zeros(size(rho{1}(:, 1)));
                
                for i = 1:star.n_datasets
                    
                    mean_score = squeeze(score{i}(:, star.pc_to_show, :))' + mean_score;
                    mean_eigenvalue = eigenvalue{i}(star.pc_to_show, :) + mean_eigenvalue;
                    mean_rho = rho{i}(:, star.pc_to_show).^2 + mean_rho;
                    rho_vector(i, :) = rho{i}(:, star.pc_to_show).^star.r_coeff_toggle(1);
                    eigenvalue_vector(i, :) = eigenvalue{i}(star.pc_to_show, :);
                    
                end
                
                mean_score = mean_score / i;
                mean_eigenvalue = mean_eigenvalue / i;
                mean_rho = mean_rho / i;
                ste_rho = std(rho_vector, 0, 1) / i;
                ste_eigenvalue = std(eigenvalue_vector, 0, 1) / i;
                
                figure(star.fig + 1000), clf
                
                subplot(2, 2, 1)
                plot(mean_score, 'LineWidth', 2)
                hold on
                minmax = ylim;
                axis([0 star.endtime minmax])
                %         xlabel('Time bin', 'FontSize', 14)
                ylabel('PC value', 'FontSize', 14)
                title(['PC projection: ' star.model_toggle ' ' star.include_states_toggle], 'FontSize', 16)
                
                subplot(2, 2, 2)
                plot(repmat(star.palatability_vals, 1, n_bins(1, 1))', 'LineWidth', 2)
                minmax = ylim;
                axis([0 star.endtime minmax])
                %         xlabel('Time bin', 'FontSize', 14)
                ylabel('Palatability value', 'FontSize', 14)
                title('Palatability', 'FontSize', 16)
                legend(star.tastant_names)
                
                subplot(2, 2, 3)
                plot(mean_eigenvalue, 'LineWidth', 2)
                hold on
                plot(mean_eigenvalue + ste_eigenvalue, 'Color', 'b', 'LineWidth', 1)
                plot(mean_eigenvalue - ste_eigenvalue, 'Color', 'b', 'LineWidth', 1)
                axis([0 star.endtime 0 1])
                xlabel('Time bin', 'FontSize', 14)
                ylabel(texlabel('R^{2}'), 'FontSize', 14)
                title(['Var. explained PC' num2str(star.pc_to_show)], 'FontSize', 16)
                
                subplot(2, 2, 4)
                plot(mean_rho, 'LineWidth', 2)
                hold on
                plot((mean_rho + ste_rho').^star.r_coeff_toggle(2), 'Color', 'b', 'LineWidth', 1)
                plot((mean_rho - ste_rho').^star.r_coeff_toggle(2), 'Color', 'b', 'LineWidth', 1)
                
                %         plot(abs(diff(mean_rho)), 'r', 'LineWidth', 2)
                axis([0 star.endtime 0 1])
                xlabel('Time bin', 'FontSize', 14)
                ylabel(texlabel('r^{2}'), 'FontSize', 14)
                title(['Corr. PC' num2str(star.pc_to_show) ' and Pal.'], 'FontSize', 16)
                
                %         addtitle('Data', 24, [0.5 0.1]);
                
        end
        
    case 'L2'
        %% Correlate L2 and palatability for time series
        
        % Neural data
        count1 = 2;
        count2 = 1;
        count3 = 1;
        count4 = star.n_tastants;
        
        for i = 1:star.n_datasets
            
            for q = 1:star.n_tastants - 1
                
                for j = count1:count4
                    
                    psth_data_new = psth_data_double(i, j, :, 2:end);
                    psth_data_offset_new = psth_data_double(i, j - count2, :, 2:end);
                    
                    tastantdistance_data(i, count3, :) = sqrt(sum((squeeze(psth_data_new - psth_data_offset_new)).^2, 2));
                    
                    pairs{count3} = [star.tastant_names{j} '-' star.tastant_names{j - count2}];
                    
                    count3 = count3 + 1;
                    
                end
                
                count1 = count1 + 1;
                count2 = count2 + 1;
                
            end
            
            count1 = 2;
            count2 = 1;
            count3 = 1;
            
        end
        
        % Behavioral data
        count1 = 2;
        count2 = 1;
        count3 = 1;
        count4 = star.n_tastants;
        
        for i = 1:star.n_datasets
            
            for q = 1:star.n_tastants - 1
                
                for j = count1:count4
                    
                    palatability_vals = star.palatability_vals(j);
                    palatability_vals_offset = star.palatability_vals(j - count2);
                    
                    tastantdistance_behav(i, count3) = sqrt(sum((squeeze(palatability_vals - palatability_vals_offset)).^2, 2));
                    
                    pairs{count3} = [star.tastant_names{j} '-' star.tastant_names{j - count2}];
                    
                    count3 = count3 + 1;
                    
                end
                
                count1 = count1 + 1;
                count2 = count2 + 1;
                
            end
            
            count1 = 2;
            count2 = 1;
            count3 = 1;
            
        end
        
        % Correlate neural and behavioral measures
        for i = 1:star.n_datasets

            [r p] = corr(squeeze(tastantdistance_data(i, :, :)), squeeze(tastantdistance_behav(i, :))', 'type', star.correlation_toggle); % calculate correlation
            
            rho(i, :) = r;
            pvalue(i, :) = p;
            [dummy star.best_index(i)] = max(rho(i, :).^2);
            star.best_index(i) = min(star.best_index(i) + star.mag_perm_toggle, star.endtime);
            
        end

        mean_rho = mean(rho.^2, 1);
        ste_rho = std(rho.^2, 0, 1) / i;
        
        % Plot
        tastantdistance_data = squeeze(mean(tastantdistance_data, 1))';
        tastantdistance_data = tastantdistance_data(star.starttime:star.endtime, :);
        maxdist_data = max(tastantdistance_data(:));
        
        switch star.plot_normalized_distance
            case('on')
                tastantdistance_data = tastantdistance_data/maxdist_data;
        end
        
        figure(star.fig + 1000)
        
        subplot(2, 2, 1)
        % imagesc(tastantdistance_data, [0 1]), colorbar
        imagesc(tastantdistance_data), colorbar
        xlabel('Tastant pairs', 'FontSize', 14);
        ylabel('Time bin', 'FontSize', 14);
        set(gca, 'XTick', 1:star.n_tastantcombos)
        set(gca, 'XTickLabel', pairs)
        title(['Model reconstructions: ' star.model_toggle ' ' star.include_states_toggle], 'FontSize', 16)
        axis([0 star.n_tastantcombos+1 0 star.endtime - star.starttime + 1])
        
        subplot(2, 2, 2)
        imagesc(repmat(tastantdistance_behav, length(star.starttime:star.endtime), 1)), colorbar
        xlabel('Tastant pairs', 'FontSize', 14);
        ylabel('Time bin', 'FontSize', 14);
        set(gca, 'XTick', 1:star.n_tastantcombos)
        set(gca, 'XTickLabel', pairs)
        title(['Behavioral data: ' star.model_toggle ' ' star.include_states_toggle], 'FontSize', 16)
        axis([0 star.n_tastantcombos+1 0 star.endtime - star.starttime + 1])
        
        subplot(2, 2, 3)
        plot(mean_rho, 'LineWidth', 2)
        hold on
        plot((mean_rho + ste_rho).^star.r_coeff_toggle(2), 'Color', 'b', 'LineWidth', 1)
        plot((mean_rho - ste_rho).^star.r_coeff_toggle(2), 'Color', 'b', 'LineWidth', 1)
        axis([0 star.endtime 0 1])
        xlabel('Time bin', 'FontSize', 14)
        ylabel(texlabel('r^{2}'), 'FontSize', 14)
        title(['Corr. PC' num2str(star.pc_to_show) ' and Pal.'], 'FontSize', 16)
        
        addtitle('Data', 24, [0.5 0.1]);
        
end