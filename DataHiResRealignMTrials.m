function star = DataHiResRealignMTrials(star)
%% function 'star = DataHiResRealignMTrials(star)'
% Reads-in neural spike train data for N trials, computes the peri-stimulus
% histogram for each tastant, then computes the correlation between
% selected PC and behavioral palatability data at each time bin, selecting
% the best correlating PC (HiRes = the use of all spikes, rather than the
% pre-processed inputs to HMM) (Mtrials = trials in which the target state
% actually came on)

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Read in data
count = 1;

for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        resultfilename = getfield(star.resultfiles(count), 'name'); % get name
        load(resultfilename, 'par'); % load structure
        
        data_matrices{i, j} = par.newdat;
        n_trials(i, j) = star.n_trials(i, j); % note use of 'star'
        n_symbols(i, j) = par.n_symbols;
        n_bins(i, j) = size(par.newdat, 2);
        count = count + 1;
        
        clear par
        
    end
    
end

%% Realign neural data to HMM state onsets/offsets
star.startpoint_realignment = star.startpoint_realignment * star.hires_factor;
for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
%        onsettimes_test = squeeze(star.state_onset_times{i, j}(n))
        data_to_use = data_matrices{i, j};
        
        for n = 1:n_trials(i, j)

            temp(n, :) = zeros(2*n_bins(i, j), 1);
            realigned_trial_start = n_bins(i, j) - star.state_onset_times{i, j}(n) * star.hires_factor;
            temp(n, realigned_trial_start:realigned_trial_start + n_bins(i, j) - 1) = data_to_use(star.new_trial_indices{i, j}(n), :);
            realigned_spikes(i, j, n, :) = temp(n, n_bins(i, j) - star.startpoint_realignment + (1:n_bins(i, j)) - 1);
            
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
            
            psth_data{i, j, k} = mean(data_matrices{i, j} == k); % average spikes over trials
            
        end
        
    end
    
end

%% Smooth neural data
star.smooth_factor = star.smooth_factor * star.hires_factor;
for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        psth_data_double(i, j, :, :) = squeeze(cell2mat(psth_data(i, j, :))); % convert to double
        
        for k = 1:star.max_n_symbols
            
            psth_data_double(i, j, :, k) = smooth(squeeze(psth_data_double(i, j, :, k)), star.smooth_factor); % average (smooth) over time bins
            
        end
        
    end
    
end

%% Plot most likely states at each time bin
star.best_index = star.best_index * star.hires_factor;
switch star.plot_psth_toggle
    case('on')
        
        for i = 1:star.n_datasets
            
            for j = 1:star.n_tastants
                
                figure(star.fig + 1000), clf
                
                subplot(2, 2, 1)
                plot(squeeze(psth_data_double(i, j, :, :)))
                minmax = ylim;
                line([star.startpoint_realignment star.startpoint_realignment], minmax, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1)
                xlabel('Time bin', 'FontSize', 14)
                ylabel('Spike probability', 'FontSize', 14)
                title('Aligned data PSTH', 'FontSize', 16)
                
                subplot(2, 2, 2)
                plot(star.psth_model{i, j})
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
                
                pause(3)
                
            end
            
        end
        
end

%% Histogram (best indices - state onsets)
switch star.plot_realignstats_toggle
    case('on')
        
        for i = 1:star.n_datasets
            
            for j = 1:star.n_tastants
                
                figure(star.fig + 10000 + i)
                
                subplot(2, 3, j)
                star.shiftdist{i, j} = star.best_index(i) - star.state_onset_times{i, j};
                hist(star.shiftdist{i, j}, 2*n_bins(i, j))
                axis([-star.endtime star.endtime 0 2])
                xlabel('Distance', 'FontSize', 14)
                ylabel('Frequency', 'FontSize', 14)
%                 title('State onset - best index', 'FontSize', 16)
                
                addtitle(['Dataset ' num2str(i) ' Tastant ' star.tastant_names(j)], 12);
                                
            end
            
        end
        
end

%% Correlate PCA and palatability for time series
star.starttime = 1;
star.endtime = star.endtime * star.hires_factor * 2;

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
            
        end
        
        mean_score = mean_score / i;
        mean_eigenvalue = mean_eigenvalue / i;
        mean_rho = mean_rho / i;
        
        figure(star.fig + 1000), clf
        
        subplot(2, 2, 1)
        plot(mean_score, 'LineWidth', 2)
        hold on
        minmax = ylim;
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
        plot(mean_eigenvalue, 'LineWidth', 2)
        hold on
        axis([0 star.endtime 0 1])
        xlabel('Time bin', 'FontSize', 14)
        ylabel(texlabel('R^{2}'), 'FontSize', 14)
        title(['Var. explained PC' num2str(star.pc_to_show)], 'FontSize', 16)
        
        subplot(2, 2, 4)
        plot(mean_rho, 'LineWidth', 2)
        hold on
        plot(abs(diff(mean_rho)), 'r', 'LineWidth', 2)
        axis([0 star.endtime 0 1])
        xlabel('Time bin', 'FontSize', 14)
        ylabel(texlabel('r^{2}'), 'FontSize', 14)
        title(['Corr. PC' num2str(star.pc_to_show) ' and Pal.'], 'FontSize', 16)
        
        addtitle('Data', 24, [0.5 0.1]);

end
