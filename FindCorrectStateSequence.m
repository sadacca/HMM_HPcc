function [hist_sequence_indices max_LL_indices_trials order_trials hist_store_order] = ...
    FindCorrectStateSequence(ntrials, nobs, nsymbols, nstates, sequence, ...
    consecutive_bins_threshold, include_repeat_states_toggle)
%% function '[hist_sequence_indices max_LL_indices_trials order_trials ...
% hist_store_order] = FindCorrectStateSequence(ntrials, nobs, nsymbols, nstates, sequence)'

% Finds the most common observed state sequence across trials, for a given
% tastant. The outputs are indices to rows of the emission matrices for
% given recording sessions (datasets)

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Main loop
for trial_loop = 1:ntrials;
    
    seq_orig = squeeze(sequence(:, :, trial_loop));
    order = zeros(1, nstates);
    state_counter = 2;
    consecutive_bins_counter = 1;
    [dummy indices] = max(seq_orig);
    order(1) = indices(1);
    
    for loop = 2:nobs
        switch include_repeat_states_toggle
            case('off')
                if ~ismember(indices(loop), order) & consecutive_bins_counter >= consecutive_bins_threshold
                    order(state_counter) = indices(loop);
                    state_counter = state_counter + 1;
                    consecutive_bins_counter = 1;
                else
                    consecutive_bins_counter = consecutive_bins_counter + 1;
                end
            case('on')
                if indices(loop)~=indices(loop - 1) & consecutive_bins_counter >= consecutive_bins_threshold
                    order(state_counter) = indices(loop);
                    state_counter = state_counter + 1;
                    consecutive_bins_counter = 1;
                else
                    consecutive_bins_counter = consecutive_bins_counter + 1;
                end
        end
        if state_counter == nstates + 1
            break
        end
    end
    
    max_LL_indices_trials(trial_loop, :) = indices;
    order_trials(trial_loop, :) = order;
    
end

%% Fill-in empty slots

for trial_loop = 1:ntrials
    
    for state_loop = 1:nstates
        
        if order_trials(trial_loop, state_loop) == 0
            order_trials(trial_loop, state_loop) = order_trials(trial_loop, state_loop - 1) ;
        end
        
    end
    
end

hist_store_order = hist(order_trials, 1:nstates);

%% Invert

for correction_loop = 1:nstates
    
    n_zero_elements(correction_loop) = ntrials - nnz(order_trials(:, correction_loop));
    
end

%% Histogram

hist_store_order(1, :) = hist_store_order(1, :) - n_zero_elements;
hist_store_order = hist_store_order / ntrials;
[dummy hist_sequence_indices] = max(hist_store_order);
