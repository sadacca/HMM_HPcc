function CompHiddMarkStates(star)
%% function 'CompHiddMarkStates'
% Reads-in parameters for each HM model type, and computes the distances
% between tastant pairs per state

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Read in HMM
count = 1;

for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        resultfilename = getfield(star.resultfiles(count), 'name'); % get name
        load(resultfilename, 'par'); % load structure
        
        LL_transition_matrices{i, j} = par.LL_transition_matrix;
        LL_emission_matrices{i, j} = par.LL_emission_matrix;
        LL_sequence{i, j} = par.hmm_LL_sequence;
        I_LL_max(i, j) = par.I_LL_max;
        
        AIC_transition_matrices{i, j} = par.AIC_transition_matrix;
        AIC_emission_matrices{i, j} = par.AIC_emission_matrix;
        AIC_sequence{i, j} = par.hmm_AIC_sequence;
        I_AIC_min(i, j) = par.I_AIC_min;
        
        BIC_transition_matrices{i, j} = par.BIC_transition_matrix;
        BIC_emission_matrices{i, j} = par.BIC_emission_matrix;
        BIC_sequence{i, j}= par.hmm_BIC_sequence;
        I_BIC_min(i, j) = par.I_BIC_min;
        
        data_matrices{i, j} = par.observations;
        data_matrices_norm{i, j} = par.norm_observations;
        n_trials(i, j) = par.n_trials;
        n_observations(i, j) = par.n_observations;
        n_symbols(i, j) = par.n_symbols;
        n_bins(i, j) = par.n_bins;
        count = count + 1;
        
        clear par
        
    end
    
end

%% Compute LL, AIC- or BIC-based emission matrices

for i = 1:star.n_datasets
    
    for j = 1:star.n_tastants
        
        switch star.model_toggle
            case('LL')
                transition = LL_transition_matrices{i, j};
                emission = LL_emission_matrices{i, j};
                emission_matrices{i, j} = LL_emission_matrices{i, j};
                sequence = LL_sequence{i, j};
                n_states(i) = min(I_LL_max(i, :));
                index = I_LL_max(i, j);
            case('AIC')
                transition = AIC_transition_matrices{i, j};
                emission = AIC_emission_matrices{i, j};
                emission_matrices{i, j} = AIC_emission_matrices{i, j};
                sequence = AIC_sequence{i, j};
                n_states(i) = min(I_AIC_min(i, :));
                index = I_AIC_min(i, j);
            case('BIC')
                transition = BIC_transition_matrices{i, j};
                emission = BIC_emission_matrices{i, j};
                emission_matrices{i, j} = BIC_emission_matrices{i, j};
                sequence = BIC_sequence{i, j};
                n_states(i) = min(I_BIC_min(i, :));
                index = I_BIC_min(i, j);
        end
        
%         [hist_sequence_indices max_PostProb_indices_trials order_alltrials hist_store_order] = ...
%             FindCorrectStateSequence(n_trials(i, j), n_observations(i, j), n_symbols(i, j), ...
%             star.max_seq_length, sequence, star.consecutive_bins_threshold, ...
%             star.include_repeat_states_toggle);

        [modifier_matrix target_state_sequence most_likely_states] = FindPrimaryStateSequence(n_trials(i, j), n_bins(i, j), index, sequence);
                        
        compressed_sequence = CompressPrimaryStateSequence(n_observations(i, j), star.max_seq_length, ...
            target_state_sequence, star.consecutive_bins_threshold, star.include_repeat_states_toggle);

        switch star.trial_average_toggle
            case('on')
                accumulating_sum = 0;
                for k = 1:n_trials(i, j)
                    accumulating_sum = accumulating_sum + emission(order_alltrials(k, :), :);
                end
                emission_seq_adjusted{i, j} = accumulating_sum / n_trials(i, j);
            case('off')
                emission_seq_adjusted{i, j} = emission(compressed_sequence, :);
        end
        
        store_target_state_sequence(i, j, :) = compressed_sequence;
        store_sequence_indices(i, j, :) = target_state_sequence;
        
    end
    
end


%% Compute distances between all tastant pairs
count1 = 2;
count2 = 1;
count3 = 1;
count4 = star.n_tastants;

for i = 1:star.n_datasets
    
    for q = 1:star.n_tastants - 1
        
        for j = count1:count4
            
            tastantdistance{count3} = sqrt(sum((emission_seq_adjusted{i, j} - emission_seq_adjusted{i, j - count2}).^2, 2));
            
            pairs{count3} = [star.tastant_names{j} '-' star.tastant_names{j - count2}];
            
            count3 = count3 + 1;
            
        end
        
        count1 = count1 + 1;
        count2 = count2 + 1;
        
    end
    
    store_tastantdistance{i} = cell2mat(tastantdistance);
    store_tastantdistance_matrix(i, :, :) = store_tastantdistance{i};
    
    count1 = 2;
    count2 = 1;
    count3 = 1;
    
end

%% Plot
switch star.plot_indiv_datasets_toggle
    case('on')
        for i = 1:star.n_datasets
            figure(i + 100)
            individual_datasets = store_tastantdistance{i};
            individual_datasets = individual_datasets / max(individual_datasets(:));
            imagesc(individual_datasets), colormap('jet'), colorbar
            ylabel('State', 'FontSize', 14);
            xlabel('Tastant pairs', 'FontSize', 14);
            set(gca, 'YTick', 1:n_states(i))
            set(gca, 'YTickLabel', 1:n_states(i))
            set(gca, 'XTick', 1:star.n_tastantcombos)
            set(gca, 'XTickLabel', pairs)
            title('Tastant-pair distance', 'FontSize', 16)
        end
end

tastantdistance_mean = squeeze(mean(store_tastantdistance_matrix, 1));
maxdist = max(tastantdistance_mean(:));
switch star.plot_normalized_distance
    case('on')
        tastantdistance_mean = tastantdistance_mean/maxdist;
end

figure(star.fig)
imagesc(tastantdistance_mean), colormap('jet'), colorbar
ylabel('State', 'FontSize', 14);
xlabel('Tastant pairs', 'FontSize', 14);
set(gca, 'YTick', 1:max(n_states))
set(gca, 'YTickLabel', 1:max(n_states))
set(gca, 'XTick', 1:star.n_tastantcombos)
set(gca, 'XTickLabel', pairs)
title([star.model_toggle ' tastant-pair distances'], 'FontSize', 16)

%% Plot most common state sequences for each tastant
for i = 1:star.n_datasets
    
    data_to_plot = squeeze(store_target_state_sequence(i, :, :));
    figure(star.fig + 100 + i)
    imagesc(data_to_plot')
    ylabel('State sequences', 'FontSize', 14);
    xlabel('Tastant', 'FontSize', 14);
    set(gca, 'YTick', 1:star.max_seq_length)
    set(gca, 'YTickLabel', 1:star.max_seq_length)
    set(gca, 'XTick', 1:star.n_tastants)
    set(gca, 'XTickLabel', star.tastant_names)
    title([star.model_toggle ' state sequence'], 'FontSize', 16)
    
    data_to_plot = squeeze(store_sequence_indices(i, :, :));
    figure(star.fig + 200 + i)
    imagesc(data_to_plot')
    ylabel('Time-series state sequence', 'FontSize', 14);
    xlabel('Tastant', 'FontSize', 14);
    set(gca, 'XTick', 1:star.n_tastants)
    set(gca, 'XTickLabel', star.tastant_names)
    title([star.model_toggle ' states per time bin'], 'FontSize', 16)
    
end

%% Correlate PCA and palatability for HMM states
switch star.do_pca_hmm_toggle
    case('on')
        
        for i = 1:star.n_datasets
            
            for j = 1:star.n_tastants
                
                emission = emission_seq_adjusted{i, j};

                for k = 1:star.max_seq_length
                
                target_states_emission(i, j, k, :) = squeeze(emission(k, :));
                
                end
                
            end
            
            for k = 1:star.max_seq_length
                      
            [coeff score latent] = princomp(squeeze(target_states_emission(i, :, k, :))); % do the pca
            palatability_vals = [11.69, 16.61, 8.375, -35.61, 10.48, -72.34]'; % mean # licks / 15s to tastant - mean # licks / 15s to H20
            [rho(i, k) pvalue(i, k)] = corr(squeeze(score(:, 1, :)), palatability_vals, 'type', star.correlation_toggle); % calculate correlation
            
            end
            
            r = rho(i, :)
            p = pvalue(i, :)
            
            figure(star.fig + i)
            barh(fliplr(rho(i, :)))
            set(gca, 'YTick', 1:star.max_seq_length)
            set(gca, 'YTickLabel', fliplr(1:star.max_seq_length))
            ylabel('State', 'FontSize', 14);
            xlabel('Rho', 'FontSize', 14);
            title([star.model_toggle ' correlation with behav. palatability'], 'FontSize', 16)
            
        end
        
end


