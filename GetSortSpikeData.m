function par = GetSortSpikeData(par, temp)
%% function 'par = GetSortSpikeData(par, temp)'
% Conditions multi-unit spike-train observations for HM analyses. If
% multiple spikes (time bins = maximal temporal resolution), choose one
% randomly. Then re-sample according to desired temporal resolution,
% according to neuron that spikes most often during epoch of interest

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Sample data to ensure HMM constraints are satisfied
if isstruct(temp)
    for i = 1:par.n_neurons
        par.rawdata(i, :, :) = temp.dat(i, par.tastant, 1:par.n_trials, :); % label spikes as symbols (1=no spikes, 2=neuron(1) spike, 3=neuron(2) spike, etc)
    end
else
    par.rawdata = temp; % read data
end

for j = 1:par.n_trials
    for k = 1:par.tottime
        tempvar = find(par.rawdata(:, j, k) > 1);
        if tempvar
            par.random = randperm(length(tempvar));
            par.newdat(j, k) = par.rawdata(tempvar(par.random(1)), j, k); % if multiple spikes, choose one randomly
        else
            par.newdat(j, k) = 1;
        end
    end
end

%% Re-sample
for j = 1:par.n_trials
    par.start = par.trialstart;
    par.fin = par.start + par.bintime - 1;
    for k = 1:par.n_bins
        par.test = hist(par.newdat(j, par.start:par.fin), 1:par.n_neurons + 1); % if bins are > 1 msec, choose most common symbol in bin
        [a b] = max(par.test(2:end)); % this constraint ensures that the 'no spikes' symbol (1) does not dominate all bins
        if a == 0
            par.resample(j, k) = b;
        else
            par.resample(j, k) = b + 1;
        end
        par.start = par.start + par.bintime;
        par.fin = par.fin + par.bintime;
    end
end

%% Output
par.observations = par.resample;
par.n_observations = length(par.observations);

% if isstruct(temp)
%     par.observations = par.resample;
%     par.n_observations = length(par.observations);
% else
%    par.dummy_observations = par.resample;
% end
