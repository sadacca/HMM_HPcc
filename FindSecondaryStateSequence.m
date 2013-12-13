function [modifier_matrix secondary_state_sequence most_likely_states] = FindSecondaryStateSequence(ntrials, nobs, nstates, sequence)
%% function [modifier_matrix secondary_state_sequence most_likely_states] = FindSecondaryStateSequence(ntrials, nobs, sequence)

% Finds the second most common observed state sequence across trials, for a
% given tastant. The outputs are indices to rows of the emission matrices
% for given recording sessions (datasets)

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Main loop
for trial_loop = 1:ntrials;
        
    seq_orig = squeeze(sequence(:, :, trial_loop));
    [dummy indices] = max(seq_orig);
    most_likely_states(trial_loop, :) = indices;
    
end

%% Histogram
hist_sequence_indices = hist(most_likely_states, 1:nstates);

%% Second Maximum
[frequency primary_state_sequence] = max(hist_sequence_indices);
for nobs_loop = 1:nobs
    
hist_sequence_indices(primary_state_sequence(nobs_loop), nobs_loop) = zeros;

end

[frequency secondary_state_sequence] = max(hist_sequence_indices);
for nobs_loop = 1:nobs
    
    if frequency(nobs_loop) == 0
        secondary_state_sequence(nobs_loop) = primary_state_sequence(nobs_loop);
    end

end

modifier_matrix = zeros(nstates, nobs);
for nobs_loop = 1:nobs
    
modifier_matrix(secondary_state_sequence(nobs_loop), nobs_loop) = ones;

end