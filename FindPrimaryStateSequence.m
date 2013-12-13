function [modifier_matrix primary_state_sequence most_likely_states] = FindPrimaryStateSequence(ntrials, nobs, nstates, sequence)
%% function [modifier_matrix primary_state_sequence] = FindPrimaryStateSequence(ntrials, nobs, sequence)

% Finds the most common observed state sequence across trials, for a given
% tastant. The outputs are indices to rows of the emission matrices for
% given recording sessions (datasets)

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

%% Maximum
[frequency primary_state_sequence] = max(hist_sequence_indices);

modifier_matrix = zeros(nstates, nobs);
for nobs_loop = 1:nobs
    
modifier_matrix(primary_state_sequence(nobs_loop), nobs_loop) = ones;

end
