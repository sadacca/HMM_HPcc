function star = DataRealignMTrials(star)
%% function 'star = DataRealignMTrials'
% Reads-in neural spike train data for N trials, realigns data to
% onset/offset of specific HMM states, then computes the peri-stimulus
% histogram for each tastant (Mtrials = trials in which the target state
% actually came on)

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
        n_trials(i, j) = star.n_trials(i, j); % note use of 'star'
        n_symbols(i, j) = par.n_symbols;
        n_bins(i, j) = par.n_bins;
        count = count + 1;
        
        clear par
        
    end
    
end

%% Realign neural data to HMM state onsets/offsets
for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        %        onsettimes_test = squeeze(star.state_onset_times{i, j}(n))
        data_to_use = data_matrices{i, j};
        
        for n = 1:n_trials(i, j)
            
            temp(n, :) = zeros(2*n_bins(i, j), 1);
            realigned_trial_start = n_bins(i, j) - star.state_onset_times{i, j}(n) + 1;
            temp(n, realigned_trial_start:realigned_trial_start + n_bins(i, j) - 1) = data_to_use(star.new_trial_indices{i, j}(n), :);
            realigned_spikes(i, j, n, :) = temp(n, n_bins(i, j) - star.startpoint_realignment + (1:n_bins(i, j)));
            
        end
        
        data_matrices{i, j} = squeeze(realigned_spikes(i, j, :, :)); % realign spikes over trials
        
        %         realigned_spikes_test = squeeze(realigned_spikes(i, j, :, :))
        %         sumrealigned_spikes_test = sum(squeeze(realigned_spikes(i, j, :, :)), 2)
        
    end
    
end

%% Analyze neural data
for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        for k = 1:star.max_n_symbols
            
            psth_data_cell{i, j, k} = mean(data_matrices{i, j} == k); % average spikes over trials
            
        end
        
    end
    
end

%% Smooth neural data
for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        psth_data_double(i, j, :, :) = squeeze(cell2mat(psth_data_cell(i, j, :))); % convert to double
        
        for k = 1:star.max_n_symbols
            
            psth_data_double(i, j, :, k) = smooth(squeeze(psth_data_double(i, j, :, k)), star.smooth_factor); % average (smooth) over time bins
            
        end
        
    end
    
end

%% Plot PSTH and most likely states at each time bin
switch star.plot_psth_toggle
    case('on')
        
        for i = 1:star.n_datasets
            
            for j = 1:star.n_tastants
                
                figure(star.fig + 100 + j), % clf
                
                subplot(2, 2, 1)
                plot(squeeze(psth_data_double(i, j, :, 2:end)))
                minmax = ylim;
                line([star.startpoint_realignment star.startpoint_realignment], minmax, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                xlabel('Time bin', 'FontSize', 14)
                ylabel('Spike probability', 'FontSize', 14)
                title('Aligned data PSTH', 'FontSize', 16)
                
                subplot(2, 2, 2)
                plot(star.psth_model{i, j}(:, 2:end))
                line([star.startpoint_realignment star.startpoint_realignment], minmax, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                title('Aligned model PSTH', 'FontSize', 16)
                
                subplot(2, 2, 3)
                imagesc(squeeze(cell2mat(star.probability_thresholded(i, j))), [0 1])
                hold on
                minmax = xlim;
                line([0 minmax(2) + 1], [star.best_index(i) star.best_index(i)], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1)
                ylabel('Time bin', 'FontSize', 14)
                xlabel('Trial', 'FontSize', 14)
                title('Unaligned post. prob.', 'FontSize', 16)
                
                subplot(2, 2, 4)
                realigned_probability = squeeze(cell2mat(star.realigned_probability(i, j)));
                realigned_probability = realigned_probability(:, star.starttime:star.endtime, :);
                data_to_plot = squeeze(realigned_probability(star.beststate(i, j), :, :));
                data_to_plot = (data_to_plot > star.threshold) .* data_to_plot;
                imagesc(data_to_plot, [0 1])
                hold on
                line([0 n_trials(i, j) + 1], [star.startpoint_realignment star.startpoint_realignment], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1)
                ylabel('Time bin', 'FontSize', 14)
                xlabel('Trial', 'FontSize', 14)
                title('Aligned post. prob.', 'FontSize', 16)
                
                addtitle(['Dataset ' num2str(i) ' Tastant ' star.tastant_names(j)], 12);
                
                pause(1)
                
            end
            
        end
        
end

%% Plot ONLY most likely states at each time bin
switch star.plot_prob_per_trial_toggle
    case('on')
        
        for i = 1:star.n_datasets
            
            figure(star.fig + 10050), clf
            
            for j = 1:star.n_tastants
                
                subplot(2, 3, j)
                imagesc(squeeze(cell2mat(star.probability_thresholded(i, j))), [0 1])
                hold on
                minmax = xlim;
                line([0 minmax(2) + 1], [star.best_index(i) star.best_index(i)], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2)
                ylabel('Time bin', 'FontSize', 14)
                if j > 3, xlabel('Trial', 'FontSize', 14), end
                title(star.tastant_names{j}, 'FontSize', 16)
                
            end
            
        end
        
        for i = 1:star.n_datasets
            
            figure(star.fig + 10150), clf
            
            for j = 1:star.n_tastants
                
                subplot(2, 3, j)
                realigned_probability = squeeze(cell2mat(star.realigned_probability(i, j)));
                realigned_probability = realigned_probability(:, star.starttime:star.endtime, :);
                data_to_plot = squeeze(realigned_probability(star.beststate(i, j), :, :));
                data_to_plot = (data_to_plot > star.threshold) .* data_to_plot;
                imagesc(data_to_plot, [0 1])
                hold on
                line([0 n_trials(i, j) + 1], [star.startpoint_realignment star.startpoint_realignment], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2)
                ylabel('Time bin', 'FontSize', 14)
                if j > 3, xlabel('Trial', 'FontSize', 14), end
                title(star.tastant_names{j}, 'FontSize', 16)
                
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
            
            star.eigenvector = squeeze(a(:, 1, star.startpoint_realignment-star.include_bins_toggle:...
                star.startpoint_realignment+star.include_bins_toggle));
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
            
            % Fit sigmoid function to correlation time series
            tempdata = r.^star.r_coeff_toggle(1);
            ydata = tempdata((star.startpoint_realignment-star.include_bins_toggle):(star.startpoint_realignment+star.include_bins_toggle));
            xdata = (-star.include_bins_toggle:star.include_bins_toggle);
            resnorm_current = inf;
            B_current = zeros(1, 4);
            for k = 1:100
                alpha = rand(1, 4);
                [B_temp resnorm_temp R exitflag output lamda J] = lsqcurvefit(@Sigmoid_fit_full, alpha, xdata, ydata, [-100 -100 -5 1], [100 100 5 1+eps]);
%                 [B_temp resnorm_temp R exitflag output lamda J] = lsqcurvefit(@Sigmoid_fit_full, alpha, xdata, ydata, [0 -1 -5 1], [1 1 5 50]);
                if resnorm_temp < resnorm_current
                    B_current = B_temp;
                    resnorm_current = resnorm_temp;
                end
            end
            B = B_current
            resnorm = resnorm_current
            sigmoid_output = Sigmoid_fit_full(B, xdata);
            temp_output = zeros(1, star.endtime); temp_output(star.startpoint_realignment+1:star.endtime) = 1;
            temp_output(star.startpoint_realignment-star.include_bins_toggle:star.startpoint_realignment+star.include_bins_toggle) = sigmoid_output;
            sigmoid_output = temp_output;
            star.data_params(i, :) = B;
            
            switch star.plot_indiv_datasets_toggle
                case('on')
                    
                    figure(star.fig + 100 + i), clf
                    
                    subplot(2, 2, 1)
                    plot(squeeze(score{i}(:, star.pc_to_show, :))', 'LineWidth', 2)
                    hold on
                    minmax = ylim;
                    line([star.startpoint_realignment star.startpoint_realignment], minmax, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
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
                    
                    addtitle('DataRealignMTrials', 24, [0.5 0.1]);
                    
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
                    mean_rho = rho{i}(:, star.pc_to_show).^star.r_coeff_toggle(1) + mean_rho;
                    rho_vector(i, :) = rho{i}(:, star.pc_to_show).^star.r_coeff_toggle(1);
                    eigenvalue_vector(i, :) = eigenvalue{i}(star.pc_to_show, :);
                    
                end
                
                mean_score = mean_score / i;
                mean_eigenvalue = mean_eigenvalue / i;
                mean_rho = mean_rho / i;
                ste_rho = std(rho_vector, 0, 1) / i;
                ste_eigenvalue = std(eigenvalue_vector, 0, 1) / i;
                
                %         tic
                %         for bin_loop = star.starttime:star.endtime
                %
                %             ci(:, bin_loop) = bootci(1, @corr, score{i}(:, star.pc_to_show, bin_loop), star.palatability_vals);
                %
                %         end
                %         toc
                %
                %         ci
                
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
                plot(sigmoid_output, 'Color', 'r', 'LineWidth', 1)
                %         plot(abs(diff(mean_rho)), 'r', 'LineWidth', 2)
                line([star.startpoint_realignment star.startpoint_realignment], [-1 1], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                axis([0 star.endtime 0 1])
                xlabel('Time bin', 'FontSize', 14)
                ylabel(texlabel('r^{2}'), 'FontSize', 14)
                title(['Corr. PC' num2str(star.pc_to_show) ' and Pal.'], 'FontSize', 16)
                
                %         addtitle('DataRealignMTrials', 24, [0.5 0.1]);
                
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
                    
                    tastantdistance_data(i, count3, :) = (sum((squeeze(abs(psth_data_new - psth_data_offset_new))).^star.L1_or_L2, 2)).^(1/star.L1_or_L2);
                    
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
                    
                    tastantdistance_behav(i, count3) = (sum((squeeze(abs(palatability_vals - palatability_vals_offset))).^star.L1_or_L2, 2)).^(1/star.L1_or_L2);
                    
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
            
            % Fit sigmoid function to correlation time series
            tempdata = r.^star.r_coeff_toggle(1);
            ydata = tempdata((star.startpoint_realignment-star.include_bins_toggle):(star.startpoint_realignment+star.include_bins_toggle));
            xdata = (-star.include_bins_toggle:star.include_bins_toggle)';
            resnorm_current = inf;
            B_current = zeros(1, 4);
            for k = 1:100
                alpha = rand(1, 4);
                                [B_temp resnorm_temp R exitflag output lamda J] = lsqcurvefit(@Sigmoid_fit_full, alpha, xdata, ydata, [-100 -100 -5 -100], [100 100 5 100]);
%                 [B_temp resnorm_temp R exitflag output lamda J] = lsqcurvefit(@Sigmoid_fit_full, alpha, xdata, ydata, [0 -1 -5 1], [1 1 5 50]);
                if resnorm_temp < resnorm_current
                    B_current = B_temp;
                    resnorm_current = resnorm_temp;
                end
            end
            B = B_current
            resnorm = resnorm_current
            sigmoid_output = Sigmoid_fit_full(B, xdata);
            temp_output = zeros(1, star.endtime); temp_output(star.startpoint_realignment+1:star.endtime) = 1;
            temp_output(star.startpoint_realignment-star.include_bins_toggle:star.startpoint_realignment+star.include_bins_toggle) = sigmoid_output;
            sigmoid_output = temp_output;
            star.data_params(i, :) = B;
            
        end
        
        mean_rho = mean(rho.^2, 1);
        ste_rho = std(rho.^2, 0, 1) / i;
        
        % Plot
        tastantdistance_data = squeeze(mean(tastantdistance_data, 1))';
        tastantdistance_data= tastantdistance_data(star.starttime:star.endtime, :);
        maxdist_data = max(tastantdistance_data(:));
        
        switch star.plot_normalized_distance
            case('on')
                tastantdistance_data= tastantdistance_data/maxdist_data;
        end
        
        figure(star.fig + 10000)
        
        imagesc(tastantdistance_data), colorbar
        xlabel('Tastant pairs', 'FontSize', 14);
        ylabel('Time bin', 'FontSize', 14);
        set(gca, 'XTick', 1:star.n_tastantcombos)
        set(gca, 'XTickLabel', pairs)
        title(['Model reconstructions: ' star.model_toggle ' ' star.include_states_toggle], 'FontSize', 16)
        axis([0 star.n_tastantcombos+1 0 star.endtime - star.starttime + 1])
        
        figure(star.fig + 1000)
        
        subplot(2, 2, 1)
        % imagesc(tastantdistance_data, [0 1]), colorbar
        %         imagesc(tastantdistance_data), colorbar
        %         xlabel('Tastant pairs', 'FontSize', 14);
        %         ylabel('Time bin', 'FontSize', 14);
        %         set(gca, 'XTick', 1:star.n_tastantcombos)
        %         set(gca, 'XTickLabel', pairs)
        %         title(['Model reconstructions: ' star.model_toggle ' ' star.include_states_toggle], 'FontSize', 16)
        %         axis([0 star.n_tastantcombos+1 0 star.endtime - star.starttime + 1])
        
        plot(tastantdistance_data)
        ylabel('Tastant pair corr.', 'FontSize', 14);
        xlabel('Time bin', 'FontSize', 14);
        title(['Model reconstructions: ' star.model_toggle ' ' star.include_states_toggle], 'FontSize', 16)
        
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
        plot(sigmoid_output, 'Color', 'r', 'LineWidth', 1)
        axis([0 star.endtime 0 1])
        xlabel('Time bin', 'FontSize', 14)
        ylabel(texlabel('r^{2}'), 'FontSize', 14)
        title(['Corr. PC' num2str(star.pc_to_show) ' and Pal.'], 'FontSize', 16)
        
        addtitle('DataRealignMTrials', 24, [0.5 0.1]);
        
end

%% Histogram (best indices - state onsets)
switch star.plot_realignstats_toggle
    case('on')
        
        for i = 1:star.n_datasets
            
            for j = 1:star.n_tastants
                
                figure(star.fig + 10000 + i)
                
                subplot(2, 3, j)
                shiftdist{i, j} = star.state_onset_times{i, j};
                star.n_trials_stateON_b4_Rbest(i, j) = sum(star.state_onset_times{i, j} < star.best_index(i), 2);
                hist(shiftdist{i, j}, n_bins(i, j))
                minmax = [0 3];
                line([star.best_index(i) star.best_index(i)], minmax, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                axis([0 star.endtime 0 2])
                xlabel('Bin number', 'FontSize', 14)
                ylabel('Count', 'FontSize', 14)
                title(star.tastant_names{j}, 'FontSize', 16)    
                %                 title('State onset - best index', 'FontSize', 16)          
                %                 addtitle(['Dataset ' num2str(i)], 16, [0.5 0.1]);
                
            end
            
        end
        
end

%% Shift distance (onset - best corr.) vs delta realigned corr.
% cover = 50;
% for i = 1:star.n_datasets
%
%     for j = 1:star.n_tastants
%
%         meanshift_per_tastant(i, j) = mean(star.best_index(i) - star.state_onset_times{i, j});
%         meanonset_per_tastant(i, j) = mean(star.state_onset_times{i, j});
%
%     end
%
%     meanshift_per_dataset(i) = mean(meanshift_per_tastant(i, :));
%     meanonset_per_dataset(i) = mean(meanonset_per_tastant(i, :));
%
%     switch star.PCA_L2_toggle
%
%         case 'PCA'
%             meancorr_pre_stateON(i) = mean(rho{i}(star.startpoint_realignment-1-cover:star.startpoint_realignment-1,...
%                 star.pc_to_show).^2);
%             meancorr_post_stateON(i) = mean(rho{i}(star.startpoint_realignment:star.startpoint_realignment+cover-1,...
%                 star.pc_to_show).^2);
%
%             stdcorr_pre_stateON(i) = std(rho{i}(star.startpoint_realignment-1-cover:star.startpoint_realignment-1,...
%                 star.pc_to_show).^2) / cover;
%             stdcorr_post_stateON(i) = std(rho{i}(star.startpoint_realignment:star.startpoint_realignment+cover-1,...
%                 star.pc_to_show).^2) / cover;
%
%         case 'L2'
%             meancorr_pre_stateON(i) = mean(rho(i, star.startpoint_realignment-1-cover:star.startpoint_realignment-1,...
%                 star.pc_to_show).^2);
%             meancorr_post_stateON(i) = mean(rho(i, star.startpoint_realignment:star.startpoint_realignment+cover-1,...
%                 star.pc_to_show).^2);
%
%             stdcorr_pre_stateON(i) = std(rho(i, star.startpoint_realignment-1-cover:star.startpoint_realignment-1,...
%                 star.pc_to_show).^2) / cover;
%             stdcorr_post_stateON(i) = std(rho(i, star.startpoint_realignment:star.startpoint_realignment+cover-1,...
%                 star.pc_to_show).^2) / cover;
%
%     end
%
% end
%
% meancorr_postpre_stateON = meancorr_post_stateON - meancorr_pre_stateON;
%
% mean_meancorr_postpre_stateON = mean(meancorr_postpre_stateON);
%
% [h p ci stats] = ttest(meancorr_postpre_stateON);
%
% data_to_use = meanshift_per_dataset;
% star.best_index;
% meanshift_per_tastant;
% meanshift_per_dataset;
%
% [r_shift p_shift] = corr(data_to_use', meancorr_postpre_stateON');
%
% figure(star.fig + 100000)
% scatter(data_to_use, meancorr_postpre_stateON, 'filled')
% Xminmax = xlim;
% Yminmax = ylim;
% xlabel('Number bins (best corr. bin - state onset)', 'FontSize', 14)
% ylabel('Mean corr. diff. (post - pre align)', 'FontSize', 14)
% axis([Xminmax Yminmax])
% r_handle = [texlabel('r^{2}') ' = ' num2str(r_shift.^2)];
% p_handle = ['p = ' num2str(p_shift)];
% text(10, 0.2, r_handle, 'FontSize', 16)
% text(10, 0.1, p_handle, 'FontSize', 16)

