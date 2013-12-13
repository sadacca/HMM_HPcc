function compressed_sequence = CompressPrimaryStateSequence(nobs, nstates, sequence, ...
    consecutive_bins_threshold, include_repeat_states_toggle)
%% function '[hist_sequence_indices max_LL_indices_trials order_trials ...
% hist_store_order = CompressPrimaryStateSequence(ntrials, nobs, nstates, sequence, ...
%     consecutive_bins_threshold, include_repeat_states_toggle)'

% Takes as input the primary (most common) state sequence across trials
% (for every time bin), for a given tastant. Outputs a 'compressed' version
% of the input, eliminating all consecutive states (i.e. in neighbouring
% bins)

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Main loop

indices = squeeze(sequence);
order = zeros(1, 30);
order(1) = indices(1);
state_counter = 1;
consecutive_bins_counter = 1;
state_length_counter(length(order)) = 0;

for loop = 2:nobs
    switch include_repeat_states_toggle
        case('off')
            if ~ismember(indices(loop), order)
                state_counter = state_counter + 1;
                order(state_counter) = indices(loop);
            else
                state_length_counter(state_counter) = state_length_counter(state_counter) + 1;
            end
        case('on')
            if indices(loop)~=indices(loop - 1)
                state_counter = state_counter + 1;
                order(state_counter) = indices(loop);
            else
                state_length_counter(state_counter) = state_length_counter(state_counter) + 1;
            end
    end
end

%% Apply consecutive bins rule

number_repeats = 0;

for state_counter = 1:length(order) - 1
    
    if state_length_counter(state_counter) < consecutive_bins_threshold
        order(state_counter) = 0;
        number_repeats = number_repeats + 1;
    end
                
end

new_order = zeros(1, nstates);
new_order(1:length(nonzeros(order)')) = nonzeros(order)';

%% Fill-in empty slots

for state_loop = 1:nstates
    
    if new_order(state_loop) == 0
        new_order(state_loop) = new_order(state_loop - 1) ;
    end
    
end

compressed_sequence = new_order(1:nstates);
