function [par temp] = GetDataParameters(par)
%% function [par temp] = GetDataParameters(par)
% Loads data and sets some default parameters for analyzing multi-unit
% spike-train observations

% example of required inputs for structure par
% par.datafilename = 'spikeGC1.mat'; 
% par.trialstart = 1; 
% par.trialduration = 5000; 
% par.bintime = 1; 
% par.tastant = 1; 

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Set parameters
temp = load(par.datafilename); % load file
par.n_neurons = size(temp.dat, 1); % # neurons recorded
par.n_symbols = par.n_neurons + 1; % # HMM symbols required
par.n_trials = size(temp.dat, 3) - 1; % # trials recorded
par.tottime = size(temp.dat, 4); % total recording time (msecs)
par.trialend = par.trialstart + par.trialduration - 1; % end trial (msecs)
par.n_bins = par.trialduration/par.bintime; % # bins