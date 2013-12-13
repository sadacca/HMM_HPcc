 function labelled_spike_array = MakeDummyData(par)
%% function 'MakeDummyData(par)'
% Creates new sets of 'dummy' spike trains based on PSTH
% data calculated at X ms resolution

% Date of creation: 2/15/2011
% Author: Tony Vladusich, Brandeis University

%% Load data
[par temp] = GetDataParameters(par); % get data params
par = GetSortSpikeData(par, temp); % get sorted spike observations
clear temp

%% Average neural data
        
for k = 1:par.n_neurons
    
    symbol = k + 1;
    psth_data(:, k) = mean(par.observations == symbol); % average spikes over trials (1=no spikes, 2=neuron(1) spike, 3=neuron(2) spike, etc)
    
end

%% Smooth neural data
        
for k = 1:par.n_neurons
    
    psth_data(:, k) = smooth(squeeze(psth_data(:, k)), par.smooth_factor); % average (smooth) over time bins
    
end

%% Generate dummy data
        
spike_array = zeros(par.n_neurons, par.n_trials, par.trialduration);

for trial = 1:par.n_trials
    
    for bin = 1:par.trialduration
        
        for neuron = 1:par.n_neurons
            
            randval = rand(1);
            if ( psth_data(bin, neuron) ) > randval
                spike_array(neuron, trial, bin) = 1; % create spike array from PSTH data
            end
            
        end
        
    end
    
end

for trial = 1:par.n_trials
    
    for neuron = 1:par.n_neurons
        
        labelled_spike_array(neuron, trial, :) = neuron * spike_array(neuron, trial, :) + 1; % label spikes as symbols (1=no spikes, 2=neuron(1) spike, 3=neuron(2) spike, etc)
        
    end
        
end

%% Plot figures

switch par.PlottingToggle
    case('on')
        
        par = GetSortSpikeData(par, labelled_spike_array);

        close all
        
        for j = 1:par.n_trials
            
            figure(j), clf
            plot(par.observations(j, :), 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 3.5)
            axis([0 par.n_observations + 1 1.1  par.n_neurons + 1.1]);
            set(gca, 'YTick', [0:par.n_neurons])
            set(gca, 'YTickLabel', [0:par.n_neurons])
            addxlabel('Time bin', 16, [0.5 0.005]);
            addylabel('Neurons', 14, [0.05 0.5]);
            
            figure(j + 100), clf
            plot(par.dummy_observations(j, :), 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 3.5)
            axis([0 par.n_observations + 1 1.1  par.n_neurons + 1.1]);
            set(gca, 'YTick', [0:par.n_neurons])
            set(gca, 'YTickLabel', [0:par.n_neurons])
            addxlabel('Time bin', 16, [0.5 0.005]);
            addylabel('Neurons', 14, [0.05 0.5]);
            
        end
        
        tilefigs
        
end
