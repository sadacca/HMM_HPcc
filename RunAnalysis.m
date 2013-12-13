function [data_slope data_P0 data_delP data_T0 data_delT data_best_indices] = RunAnalysis
%% function 'RunAnalysis'
% Computes neural time series that have been re-aligned to start of the HMM
% state that correlates best with palatability measure

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Initialize parameters
close all, clear all
star.basedirname = '\HMM'; % base directory name
star.datasetdirname = '\Sadacca'; % dataset directory name (Sadacca, Fontanini)
star.resultdirname = '\Results\newresults1ms'; % results directory name
star.simtype = '\AllData'; % all data or N-fold cross-validation (AllDummyData, AllData)
addpath([star.basedirname star.datasetdirname star.resultdirname star.simtype]); % add path
star.resultfiles = dir([star.basedirname star.datasetdirname star.resultdirname star.simtype '\spikeGC4.mat*.mat']); % use single dataset
% star.resultfiles = dir([star.basedirname star.datasetdirname star.resultdirname star.simtype '/spikeGC*.mat']); % use all datasets in dir
star.n_resultfiles = length(star.resultfiles);
switch star.datasetdirname
    case '\Fontanini'
        star.n_datasets = 12;
        star.n_tastants = 4;
        star.tastant_names = {'S' 'Q' 'C' 'N'}; % note: code won't work for this dataset (no behavioral measures available)
    case '\Sadacca'
        star.n_tastants = 6;
        star.n_datasets = star.n_resultfiles / star.n_tastants;
        star.tastant_names = {'N1' 'N2' 'N3' 'N4' 'S' 'Q'};
end
star.PCA_L2_toggle = 'PCA'; % analysis type (PCA, L2)
star.L1_or_L2 = 2; % 1 = L1 norm, 2 = L2 norm
star.plot_normalized_distance = 'on'; % normalize neural difference vector
star.n_tastantcombos = (star.n_tastants - 1)*star.n_tastants/2; % # tastant pairs
star.max_n_symbols = 15; % # symbols (= # neurons + 1)
star.plot_psth_toggle = 'on'; % plot PSTH and target states per trial, pre- and post re-alignment (on, off)
star.plot_prob_per_trial_toggle= 'on'; % plot ONLY target states per trial, pre- and post re-alignment (on, off)
star.plot_indiv_datasets_toggle = 'on'; % (on, off)
star.plot_realignstats_toggle = 'on'; % plots hist. of distance between state onsets and best correlated bin (on, off)
star.state_combo_toggle = 'on'; % computes set of states across tastants that best correlates with pal. (on, off)
star.r_coeff_toggle = [2 1]; % sets coefficients for individual and mean correlations to plot (i.e r or r^2)
star.mag_perm_toggle = 0; % adds or subtracts X index values to best indices, as reality check
star.smooth_factor = 10; % PSTH smoothing filter width
star.hires_factor = 10; % multiplies parameter values for HiRes analysis
star.starttime = 1; % trial start bin
star.endtime = 250; % trial end bin
star.plot_mostlikelystates_toggle = 'on'; % plots most likely state at each time bin (for each dataset and tastant)
star.correlation_toggle = 'Pearson'; % Pearson or Spearman
star.pc_to_show = 1; % (1, 2, 3)
star.threshold = 0.5; % defines when states comes 'on'
star.startpoint_realignment = 100; % arbitrary re-alignment point
star.model_toggle = 'AIC'; % type of model (BIC, AIC, LL)
star.include_states_toggle = 'primary'; % commonality of state (primary, secondary)
temp_palatability_vals = -[11.69, 16.61, 8.375, -35.61, 10.48, -72.34]';
star.palatability_vals = temp_palatability_vals(1:star.n_tastants); % mean # licks / 15s to tastant - mean # licks / 15s to H20
star.include_bins_toggle = 50; % number of bins to include in sigmoid fits pre- and post-alignment

%% Pre-alignment

% star.fig = 1; 
% Model(star);

star.fig = 11; star.select_trial_flag = 0; % selects only trials in which target state appears (0 = all trials, 1 = target trials)
star = Data(star); 
data_best_indices = star.best_index;

%% Post-alignment, trials in which target state doesn't appear omitted

star.fig = 21;
star = ModelRealignMTrials(star);
% orig_trials_per_dataset = sum(star.orig_n_trials, 2)'
% stateON_trials_per_dataset = sum(star.n_trials, 2)'

star.fig = 31;
star = DataRealignMTrials(star);
data_P0 = star.data_params(1, 1);
data_delP = star.data_params(1, 2);
data_T0 = star.data_params(1, 3);
data_delT = star.data_params(1, 4);
data_slope = data_delP / data_delT;


% stateON_trials_per_dataset_b4_Rbest = sum(star.n_trials_stateON_b4_Rbest, 2)'

% star.fig = 41; star.select_trial_flag = 1; % selects only trials in which target state appears (0 = all trials, 1 = target trials)
% star = Data(star);

%% Post-alignment 

% star.fig = 21;
% star = ModelRealign(star);
% 
% star.fig = 31;
% DataRealign(star);

% star.fig = 31;
% star = DataHiRes(star);
% best_indices = star.best_index

%% Post-alignment, and shuffle dataset indices as a 'reality check'

% star.best_index = star.best_index(randperm(length(star.best_index))); 
% best_indices = star.best_index
% star.fig = 61;
% star = ModelRealignMTrials(star);
% 
% star.fig = 71;
% star = DataRealignMTrials(star)

%% Post-alignment, and re-compute best-correlating indices using only 'ON' trials then re-align

% star.fig = 61; star.select_trial_flag = 1;
% star = Data(star); best_indices = star.best_index
% 
% star.fig = 71;
% star = ModelRealignMTrials(star);
% 
% star.fig = 81;
% star = DataRealignMTrials(star);

% star.best_index = ones(1, star.n_datasets) * 250;


rmpath([star.basedirname star.datasetdirname star.resultdirname star.simtype]); % rm path
