function par = ReadPlotFile(par, PlottingToggle)
%% function 'par = ReadPlotFile(par, PlottingToggle)'
% Loads and plots original spiking observations

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Get data
[par temp] = GetDataParameters(par); % data params
par = GetSortSpikeData(par, temp); % get sorted data

%% Plot figures
switch PlottingToggle
    case('on')
        
        close all
        
        for j = 1:par.n_trials
            
            figure(j), clf
            plot(par.observations(j, :), 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 3.5)
            axis([0 par.n_observations + 1 1.1  par.n_neurons + 1.1]);
            set(gca, 'YTickLabel', [1:par.n_neurons])
            addxlabel('Time bin', 16, [0.5 0.005]);
            addylabel('Neurons', 14, [0.05 0.5]);
            
        end
        
        tilefigs
        
end

