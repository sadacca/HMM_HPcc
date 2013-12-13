function star = FindBestIndices(star)
%% function 'star = FindBestIndices'
% Reads-in neural spike train data for N trials, and finds the emission
% matrix indices that best correlate with palatability.

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

%% Compute LL, AIC- or BIC-based PSTHs
switch star.model_toggle
    case('LL')
        emission = LL_emission_matrices;
        index = I_LL_max;
    case('AIC')
        emission = AIC_emission_matrices;
        index = I_AIC_min;
    case('BIC')
        emission = BIC_emission_matrices;
        index = I_BIC_min;
end

%% Compute new emission matrices permuting all tastant-state combos
count = 0;

for i = 1:star.n_datasets
    
    for m1 = 1:index(i, 1)
        
        for m2 = 1:index(i, 2)
            
            for m3 = 1:index(i, 3)
                
                for m4 = 1:index(i, 4)
                    
                    for m5 = 1:index(i, 5)
                        
                        for m6 = 1:index(i, 6)
                            
                            count = count + 1;
                            indices(:, count, i) = [m1 m2 m3 m4 m5 m6];
                            star.emission_permutation{:, :, count, i} = [emission{i, 1}(m1, 2:end); emission{i, 2}(m2, 2:end); ...
                                emission{i, 3}(m3, 2:end); emission{i, 4}(m4, 2:end); emission{i, 5}(m5, 2:end); ...
                                emission{i, 6}(m6, 2:end)];
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    total_combos(i) = count;
    count = 0;
    
end
total_combos

%% PCA and correlation with palatability
for i = 1:star.n_datasets
    
    for count = 1:total_combos(i)
        
        [a(:, :, count) b(:, :, count) latent] = princomp(star.emission_permutation{:, :, count, i}); % do the pca
        c(:, count) = latent ./ sum(latent);
        
    end
    
    coeff{i} = a;
    score{i} = b;
    eigenvalue{i} = c;
    clear a b c
    
end



for i = 1:star.n_datasets
    
    for pc_loop = 1:1
        
        [r(:, pc_loop) p(:, pc_loop)] = corr(score{i}(:, pc_loop, :), star.palatability_vals, 'type', star.correlation_toggle); % calculate correlation
        
    end
    
    rho{i} = r;
%     star.corr_states_pal(i, :) = r;
    pvalue{i} = p;
    clear r p
    
    [dummy best_count(i)] = max(rho{i});
    star.beststate(i, :) = indices(:, best_count(i), i)';
    
    figure(2)
    plot(smooth(rho{i}, 1))
    
end

