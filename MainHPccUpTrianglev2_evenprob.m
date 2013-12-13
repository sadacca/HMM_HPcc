function MainHPccUpTrianglev2_evenprob(par)
%% function 'MainHpcc(par)'
% Performs Hidden Markov Model (HMM) analysis on neural spike train data
% recorded from multiple units simultaneously: identifies putative hidden
% network state sequences from spike trains

% Code optimized to run in conjunction with Matlab Parallel Computing
% Toolbox

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Load data

warning off;
filename='currentHMMdata.mat';
load(filename)
[par temp] = GetDataParameters(par); % get data params
par = GetSortSpikeData(par, temp); % get sorted spike observations
delete(filename)
cd ..

%% Get parameter values from structure (parfor compatible)
thresh = par.thresh;
its = par.its;
max_HMM_states = par.max_HMM_states;
n_HM_models = par.n_HM_models;
n_HM_cycles = par.n_HM_cycles;
n_workers = par.n_workers;
observations = par.observations;
n_symbols = par.n_symbols;
n_observations = par.n_observations;
n_trials = par.n_trials;

%% HMM train and estimate
%set(0,'RecursionLimit',500)
tic

for n = 1:n_HM_models
    
    pALL = diag(0.99*ones(1,max_HMM_states(n)));
    
    if n == 1      
    else
        for ii=1:max_HMM_states(n)-1
            pALL(ii,ii+1:end)=(1-sum(pALL(ii,:)))/(length(pALL(ii,ii:end))-1);
        end   
    end
    
    pALL(end,end)=1;
    p = pALL;
    
      for m = 1:n_HM_cycles
            
        q = rand(max_HMM_states(n), n_symbols); % emission matrix (q) randomly initialized.
        q = normalise(q, 2); % normalise
        
        [Pp Qq LL] = hmmtrain(observations, p, q, 'Verbose', false,'Maxiterations', its, ...
            'Algorithm', 'Baumwelch', 'Tolerance', 1e-6); % parameter estimation
        
        n_transition_params(n, m) = nnz(Pp) - size(Pp, 2) + sum(Pp(:) == 1); % # free transition matrix params
        n_emission_params(n, m) = nnz(Qq) - size(Qq, 2) + sum(Qq(:) == 1); % # free emission matrix params
        
        AIC_hmmtrain(n, m) = -2*LL(end) + 2*(n_transition_params(n, m) + n_emission_params(n, m)); % store AIC values
        BIC_hmmtrain(n, m) = -2*LL(end) + ...
            (n_transition_params(n, m) + n_emission_params(n, m))*log(n_observations*n_trials); % store BIC values
        
        if size(Pp, 1) < max_HMM_states(end)
            Pp(max_HMM_states(end), :) = 0; % pad variably-sized model matrices with zeros
            Pp(:, max_HMM_states(end)) = 0;
            Qq(max_HMM_states(end), :) = 0;
        end
        
        P(:, :, n, m) = Pp; % store parameters
        Q(:, :, n, m) = Qq;
        LL_hmmtrain(n, m) = LL(end);
        
    end
    
end

timeforonetastant = toc

%% Store variables in structure
par.n_transition_params = n_transition_params;
par.n_emission_params = n_emission_params;
par.AIC_hmmtrain = AIC_hmmtrain;
par.BIC_hmmtrain = BIC_hmmtrain;
par.P = P;
par.Q = Q;
par.LL_hmmtrain = LL_hmmtrain;

%% Select best fitting model (maximizing log likelihood)
for i = 1:par.n_HM_models;
    
    for j = 1:par.n_HM_cycles
        
        if par.LL_hmmtrain(i, j) == max(max(par.LL_hmmtrain)) % find indices to max LL model
            par.I_LL_max = i; par.J_LL_max = j; % save indices
            
            par.LL_transition_matrix = par.P(1:par.max_HMM_states(i), ...
                1:par.max_HMM_states(i), i, j); % access relevant values
            par.LL_emission_matrix = par.Q(1:par.max_HMM_states(i), :, i, j);
            
        end
        
    end
    
end

%% Select most parsimonious model (minimizing AIC)
for i = 1:par.n_HM_models;
    
    for j = 1:par.n_HM_cycles
        
        if par.AIC_hmmtrain(i, j) == min(min(par.AIC_hmmtrain)) % find indices to min AIC model
            par.I_AIC_min = i; par.J_AIC_min = j; % save indices
            
            par.AIC_transition_matrix = par.P(1:par.max_HMM_states(i), ...
                1:par.max_HMM_states(i), i, j); % access relevant values
            par.AIC_emission_matrix = par.Q(1:par.max_HMM_states(i), :, i, j);
            
        end
        
    end
    
end

%% Select most parsimonious model (minimizing BIC)
for i = 1:par.n_HM_models;
    
    for j = 1:par.n_HM_cycles
        
        if par.BIC_hmmtrain(i, j) == min(min(par.BIC_hmmtrain)) % find indices to min BIC model
            par.I_BIC_min = i; par.J_BIC_min = j; % save indices
            
            par.BIC_transition_matrix = par.P(1:par.max_HMM_states(i), ...
                1:par.max_HMM_states(i), i, j); % access relevant values
            par.BIC_emission_matrix = par.Q(1:par.max_HMM_states(i), :, i, j);
            
        end
        
    end
    
end

%% Decode best HMM solution for max LL, and min AIC/BIC models
par.biggest_symbol_value = max(max(par.observations));
par.norm_observations = par.observations/par.biggest_symbol_value; % normalise obs. for display

for j = 1:par.n_trials
    
    [hmm_sequence LL] = hmmdecode(par.observations(j, :), ...
        par.LL_transition_matrix, par.LL_emission_matrix); % decode most likely state sequence and compute associated LL
    par.hmm_LL_LL(j) = LL;
    par.hmm_LL_sequence(:, :, j) = hmm_sequence;
    
end

%% Decode best HMM (AIC) solution
for j = 1:par.n_trials
    
    [hmm_sequence LL] = hmmdecode(par.observations(j, :), ...
        par.AIC_transition_matrix, par.AIC_emission_matrix); % decode most likely state sequence and compute associated LL
    par.hmm_AIC_LL(j) = LL;
    par.hmm_AIC_sequence(:, :, j) = hmm_sequence;
    
end

%% Decode best HMM (BIC) solution
for j = 1:par.n_trials
    
    [hmm_sequence LL] = hmmdecode(par.observations(j, :), ...
        par.BIC_transition_matrix, par.BIC_emission_matrix); % decode most likely state sequence and compute associated LL
    par.hmm_BIC_LL(j) = LL;
    par.hmm_BIC_sequence(:, :, j) = hmm_sequence;
    
end

%% Define filename and save par.mat files
d1 = num2str(par.n_neurons);
d2 = num2str(par.n_trials);
d3 = num2str(par.tastant);
d4 = num2str(par.trialduration);
d5 = num2str(par.bintime);
d6 = num2str(par.max_HMM_states);
d7 = num2str(par.n_HM_cycles);

par.resultfilename = [par.datafilename '.' d1 '.' d2 '.' d3 '.' d4 '.' d5 '.' d6 '.' d7 '.mat']; % define .mat filename

cd(par.resultdir) % go to results directory
save(par.resultfilename, 'par') % save

quit

%% Plot
% ReadPlotAnalFile(par.resultfilename, 'off', 'LL');
