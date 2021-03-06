function star = ModelRealignMTrials(star)
%% function 'star = ModelRealignMTrials'
% Reads-in parameters for each HM model type, then computes the
% transition-realigned putative peri-stimulus histogram that gave rise to
% the model (Mtrials = trials in which the target state actually came on)

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Read in data
count = 1;

for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        resultfilename = getfield(star.resultfiles(count), 'name'); % get name
        load(resultfilename, 'par'); % load structure
        
        LL_transition_matrices{i, j} = par.LL_transition_matrix;
        LL_emission_matrices{i, j} = par.LL_emission_matrix;
        LL_probability{i, j} = par.hmm_LL_sequence;
        I_LL_max(i, j) = par.I_LL_max;
        
        AIC_transition_matrices{i, j} = par.AIC_transition_matrix;
        AIC_emission_matrices{i, j} = par.AIC_emission_matrix;
        AIC_probability{i, j} = par.hmm_AIC_sequence;
        I_AIC_min(i, j) = par.I_AIC_min;
        
        BIC_transition_matrices{i, j} = par.BIC_transition_matrix;
        BIC_emission_matrices{i, j} = par.BIC_emission_matrix;
        BIC_probability{i, j}= par.hmm_BIC_sequence;
        I_BIC_min(i, j) = par.I_BIC_min;
        
        data_matrices{i, j} = par.observations;
        data_matrices_norm{i, j} = par.norm_observations;
        n_trials(i, j) = par.n_trials;
        n_symbols(i, j) = par.n_symbols;
        n_bins(i, j) = par.n_bins;
        count = count + 1;
        
        clear par
        
    end
    
end

switch star.state_combo_toggle
    case('on')
        star = FindBestIndices(star);
end

%% Compute LL, AIC- or BIC-based PSTHs
for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        switch star.model_toggle
            case('LL')
                emission = LL_emission_matrices{i, j};
                probability = LL_probability{i, j};
                state_index = I_LL_max(i, j);
            case('AIC')
                emission = AIC_emission_matrices{i, j};
                probability = AIC_probability{i, j};
                state_index = I_AIC_min(i, j);
            case('BIC')
                emission = BIC_emission_matrices{i, j};
                probability = BIC_probability{i, j};
                state_index = I_BIC_min(i, j);
        end
        
        switch star.include_states_toggle
            case('primary')
                [modifier_matrix target_state_sequence most_likely_states] = FindPrimaryStateSequence(n_trials(i, j), n_bins(i, j), state_index, probability);
                state_sequence(i, j, :) = target_state_sequence;
            case('secondary')
                [modifier_matrix target_state_sequence most_likely_states] = FindSecondaryStateSequence(n_trials(i, j), n_bins(i, j), state_index, probability);
                state_sequence(i, j, :) = target_state_sequence;
        end
        
        switch star.state_combo_toggle
            case('on')
                star.best_index = ones(1, star.n_datasets);
            case('off')
                star.beststate(i, j) = squeeze(state_sequence(i, j, star.best_index(i)));
        end
        
        data_to_use = squeeze(probability(star.beststate(i, j), :, :));
        star.probability_thresholded{i, j} = (data_to_use > star.threshold) .* data_to_use;
        
        state_onset_times = zeros(1, n_trials(i, j));
        
        for n = 1:n_trials(i, j)
            
            [vals indices] = find(star.probability_thresholded{i, j}(:, n)');
            if isempty(indices)
                state_onset_times(n) = 0;
            else
                state_onset_times(n) = min(indices);
            end
            
        end
        
        [dummy star.new_trial_indices{i, j}] = find(state_onset_times);
        star.state_onset_times{i, j} = nonzeros(state_onset_times)';
        star.orig_n_trials(i, j) = n_trials(i, j);
        star.n_trials(i, j) = nnz(state_onset_times);
        
        for n = 1:star.n_trials(i, j)
            
            for m = 1:state_index

                temp(m, :, n) = zeros(1, 2*n_bins(i, j));
                realigned_trial_start = n_bins(i, j) - star.state_onset_times{i, j}(n) + 1;
%                 bb=[(realigned_trial_start) i j n state_index]
                temp(m, realigned_trial_start:realigned_trial_start + n_bins(i, j) - 1, n) = probability(m, :, star.new_trial_indices{i, j}(n));
                realigned_probability(m, :, n) = temp(m, n_bins(i, j) - star.startpoint_realignment + (1:n_bins(i, j)), n);
%                 cc = n_bins(i, j) - star.startpoint_realignment + (1:n_bins(i, j)) - 1
                emission_probability_product(n, m, :, :) = repmat(emission(m, :)', 1, n_bins(i, j)) .* repmat(realigned_probability(m, :, n) , n_symbols(i, j), 1);
                
            end
            
            emission_probability_innerproduct(n, :, :) = squeeze(sum(emission_probability_product(n, :, :, :), 2)); % sum over all states
            
        end
        
        star.realigned_probability{i, j} = realigned_probability;
        
%                 data_to_check = squeeze(realigned_probability(star.beststate(i, j), 100:250, :));
%                 data_to_check = (data_to_check > star.threshold) .* data_to_check
        
        
        spike_rate_vector_trial{i, j} = emission_probability_innerproduct;
        star.psth_model{i, j} = squeeze(mean(emission_probability_innerproduct, 1))';
        
        clear emission_probability_product emission_probability_innerproduct emission probability state_index realigned_probability
        
    end
    
end

%% Plot most likely states at each time bin
switch star.plot_mostlikelystates_toggle
    case('on')
        for dataset_loop = 1:star.n_datasets
            
            figure(star.fig + 1)
            imagesc(squeeze(state_sequence(dataset_loop, :, :))')
            xlabel('Tastant', 'FontSize', 14)
            ylabel('Time', 'FontSize', 14)
            set(gca, 'XTickLabel', star.tastant_names)
            title(['Dataset ' num2str(dataset_loop)], 'FontSize', 16)
            star.fig = star.fig + 1;
            
        end
end

%% Correlate PCA and palatability for time series
switch star.PCA_L2_toggle
    
    case 'PCA'
        
        for i = 1:star.n_datasets
            
            for j = 1:star.n_tastants
                
                psth_model_double(j, :, :) = squeeze(cell2mat(star.psth_model(i, j, :)));
                
            end
            
            for bin_loop = star.starttime:star.endtime
                
                [a(:, :, bin_loop) b(:, :, bin_loop) latent] = princomp(squeeze(psth_model_double(:, bin_loop, 2:end))); % do the pca
                c(:, bin_loop) = latent ./ sum(latent);
                
            end
            
            coeff{i} = a;
            score{i} = b;
            c(find(isnan(c))) = 0; % set eigenvalue to zero if NaN
            eigenvalue{i} = c;
            
            for pc_loop = 1:3
                
                [r(:, pc_loop) p(:, pc_loop)] = corr(score{i}(:, pc_loop, :), star.palatability_vals, 'type', star.correlation_toggle); % calculate correlation
                
            end
            
            r(find(isnan(r))) = 0; % set r to zero if NaN
            rho{i} = r;
            pvalue{i} = p;
            
            switch star.plot_indiv_datasets_toggle
                case('on')
                    
                    figure(star.fig + 100 + i), clf
                    
                    subplot(2, 2, 1)
                    plot(squeeze(score{i}(:, star.pc_to_show, :))', 'LineWidth', 2)
                    hold on
                    minmax = ylim;
                    line([star.startpoint_realignment star.startpoint_realignment], minmax, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                    axis([0 star.endtime minmax])
                    %             xlabel('Time bin', 'FontSize', 14)
                    ylabel('PC value', 'FontSize', 14)
                    title(['PC projection: ' star.model_toggle ' ' star.include_states_toggle], 'FontSize', 16)
                    
                    subplot(2, 2, 2)
                    plot(repmat(star.palatability_vals, 1, n_bins(1, 1))', 'LineWidth', 2)
                    minmax = ylim;
                    axis([0 star.endtime minmax])
                    %             xlabel('Time bin', 'FontSize', 14)
                    ylabel('Palatability value', 'FontSize', 14)
                    title('Palatability', 'FontSize', 16)
                    legend(star.tastant_names)
                    
                    subplot(2, 2, 3)
                    plot(eigenvalue{i}(star.pc_to_show, :), 'LineWidth', 2)
                    hold on
                    line([star.startpoint_realignment star.startpoint_realignment], [0 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                    axis([0 star.endtime 0 1])
                    xlabel('Time bin', 'FontSize', 14)
                    ylabel(texlabel('R^{2}'), 'FontSize', 14)
                    title(['Var. explained PC' num2str(star.pc_to_show)], 'FontSize', 16)
                    
                    subplot(2, 2, 4)
                    plot(rho{i}(:, star.pc_to_show), 'LineWidth', 2)
                    hold on
                    line([star.startpoint_realignment star.startpoint_realignment], [-1 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                    axis([0 star.endtime -1 1])
                    xlabel('Time bin', 'FontSize', 14)
                    ylabel(texlabel('r'), 'FontSize', 14)
                    title(['Corr. PC' num2str(star.pc_to_show) ' and Pal.'], 'FontSize', 16)
                    
                    addtitle('ModelRealignMTrials', 24, [0.5 0.1]);
                    
            end
            
            clear psth_model_double a b c r p
            
        end
        
        switch star.plot_indiv_datasets_toggle
            
            case('off')
                
                mean_score = zeros(size(squeeze(score{1}(:, 1, :))'));
                mean_eigenvalue = zeros(size(eigenvalue{1}(1, :)));
                mean_rho = zeros(size(rho{1}(:, 1)));
                
                for i = 1:star.n_datasets
                    
                    mean_score = squeeze(score{i}(:, star.pc_to_show, :))' + mean_score;
                    mean_eigenvalue = eigenvalue{i}(star.pc_to_show, :) + mean_eigenvalue;
                    mean_rho = rho{i}(:, star.pc_to_show).^star.r_coeff_toggle(1) + mean_rho;
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
                line([star.startpoint_realignment star.startpoint_realignment], minmax, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
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
                line([star.startpoint_realignment star.startpoint_realignment], [0 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                axis([0 star.endtime 0 1])
                xlabel('Time bin', 'FontSize', 14)
                ylabel(texlabel('R^{2}'), 'FontSize', 14)
                title(['Var. explained PC' num2str(star.pc_to_show)], 'FontSize', 16)
                
                subplot(2, 2, 4)
                plot(mean_rho.^star.r_coeff_toggle(2), 'LineWidth', 2)
                hold on
                plot((mean_rho + ste_rho').^star.r_coeff_toggle(2), 'Color', 'b', 'LineWidth', 1)
                plot((mean_rho - ste_rho').^star.r_coeff_toggle(2), 'Color', 'b', 'LineWidth', 1)
                
                %         plot(abs(diff(mean_rho)), 'r', 'LineWidth', 2)
                line([star.startpoint_realignment star.startpoint_realignment], [-1 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                axis([0 star.endtime 0 1])
                xlabel('Time bin', 'FontSize', 14)
                ylabel(texlabel('r^{2}'), 'FontSize', 14)
                title(['Corr. PC' num2str(star.pc_to_show) ' and Pal.'], 'FontSize', 16)
                
                %         addtitle('ModelRealignMTrials', 24, [0.5 0.1]);
                
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
                    
                    psth_model_new = star.psth_model{i, j};
                    psth_model_offset_new = star.psth_model{i, j - count2};
                    
                    tastantdistance_model(i, count3, :) = sqrt(sum((squeeze(psth_model_new - psth_model_offset_new)).^2, 2));
                    
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
            
            [r p] = corr(squeeze(tastantdistance_model(i, :, :)), squeeze(tastantdistance_behav(i, :))', 'type', star.correlation_toggle); % calculate correlation
            
            rho(i, :) = r;
            pvalue(i, :) = p;
            
        end
        
        mean_rho = mean(rho.^2, 1);
        ste_rho = std(rho, 0, 1) / i;
        
        % Plot
        tastantdistance_model = squeeze(mean(tastantdistance_model, 1))';
        tastantdistance_model = tastantdistance_model(star.starttime:star.endtime, :);
        maxdist_model = max(tastantdistance_model(:));
        
        switch star.plot_normalized_distance
            case('on')
                tastantdistance_model = tastantdistance_model/maxdist_model;
        end
        
        figure(star.fig + 1000)
        
        subplot(2, 2, 1)
        % imagesc(tastantdistance_model, [0 1]), colorbar
        imagesc(tastantdistance_model), colorbar
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
        
        addtitle('ModelRealignMTrials', 24, [0.5 0.1]);
        
end