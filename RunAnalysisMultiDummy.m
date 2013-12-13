function RunAnalysisMultiDummy
%% function 'RunAnalysisMultiDummy'
% Computes neural time series that have been re-aligned to start of the HMM
% state that correlates best with palatability measure

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Number of dummy sets to run
close all, clear all, %xxx = 7
number_dummy_sets = 20;
sigmoid_slopes = zeros(1, number_dummy_sets);

for dummy_set_number = 1:number_dummy_sets
    
    close all
    dummy_set_number
    %% Initialize parameters
    star.basedirname = '/Users/TheVlad/Documents/MATLAB/work/HMM'; % base directory name
    star.datasetdirname = '/Sadacca'; % dataset directory name (Sadacca, Fontanini)
    star.resultdirname = '/Results'; % results directory name
    star.simtype = '/AllDummyData'; % all data or N-fold cross-validation (AllDummyData, AllData)
    star.dummy_set_number_str = num2str(dummy_set_number); % dummy set number
    resultdir = [star.basedirname star.datasetdirname star.resultdirname star.simtype star.dummy_set_number_str]; % results dir path
    addpath(resultdir); % add path
    star.resultfiles = dir([resultdir '/spikeGC4.mat*.mat']); % use single dataset
    % star.resultfiles = dir([star.basedirname star.datasetdirname star.resultdirname star.simtype '/spikeGC*.mat']); % use all datasets in dir
    star.n_resultfiles = length(star.resultfiles);
    switch star.datasetdirname
        case '/Fontanini'
            star.n_datasets = 12;
            star.n_tastants = 4;
            star.tastant_names = {'S' 'Q' 'C' 'N'}; % note: code won't work for this dataset (no behavioral measures available)
        case '/Sadacca'
            star.n_tastants = 6;
            star.n_datasets = star.n_resultfiles / star.n_tastants;
            star.tastant_names = {'N1' 'N2' 'N3' 'N4' 'S' 'Q'};
    end
    star.PCA_L2_toggle = 'PCA'; % analysis type (PCA, L2)
    star.L1_or_L2 = 2; % 1 = L1 norm, 2 = L2 norm
    star.plot_normalized_distance = 'off'; % normalize neural difference vector
    star.n_tastantcombos = (star.n_tastants - 1)*star.n_tastants/2; % # tastant pairs
    star.max_n_symbols = 15; % # symbols (= # neurons + 1)
    star.plot_psth_toggle = 'off'; % plot PSTH and target states per trial, pre- and post re-alignment (on, off)
    star.plot_prob_per_trial_toggle= 'off'; % plot ONLY target states per trial, pre- and post re-alignment (on, off)
    star.plot_indiv_datasets_toggle = 'off'; % (on, off)
    star.plot_realignstats_toggle = 'off'; % plots hist. of distance between state onsets and best correlated bin (on, off)
    star.state_combo_toggle = 'off'; % computes set of states across tastants that best correlates with pal. (on, off)
    star.r_coeff_toggle = [2 1]; % sets coefficients for individual and mean correlations to plot (i.e r or r^2)
    star.mag_perm_toggle = 0; % adds or subtracts X index values to best indices, as reality check
    star.smooth_factor = 10; % PSTH smoothing filter width
    star.hires_factor = 10; % multiplies parameter values for HiRes analysis
    star.starttime = 1; % trial start bin
    star.endtime = 250; % trial end bin
    star.plot_mostlikelystates_toggle = 'off'; % plots most likely state at each time bin (for each dataset and tastant)
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
    dummy_best_indices(dummy_set_number) = star.best_index;
    
    %% Post-alignment, trials in which target state doesn't appear omitted
    
    star.fig = 21;
    star = ModelRealignMTrials(star);
    % orig_trials_per_dataset = sum(star.orig_n_trials, 2)'
    % stateON_trials_per_dataset = sum(star.n_trials, 2)'
    
    star.fig = 31;
    star = DataRealignMTrials(star);
    dummy_P0(dummy_set_number) = star.data_params(1, 1);
    dummy_delP(dummy_set_number) = star.data_params(1, 2);
    dummy_T0(dummy_set_number) = star.data_params(1, 3);
    dummy_delT(dummy_set_number) = star.data_params(1, 4);
    dummy_slopes(dummy_set_number) = dummy_delP(dummy_set_number) / dummy_delT(dummy_set_number);

    switch star.PCA_L2_toggle
        case 'PCA'
            eigenvectors = (abs(star.eigenvector) > 0.01);
            sum_eigenvectors(dummy_set_number) = sum(sum(eigenvectors, 1));
    end
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
    
%     tilefigs
    
    rmpath([star.basedirname star.datasetdirname star.resultdirname star.simtype star.dummy_set_number_str]); % rm path
end
% bla
[data_slope data_P0 data_delP data_T0 data_delT data_best_indices] = RunAnalysis;


[h p_slope ci stats] = ttest(dummy_slopes, data_slope); % compare dummy slopes vs original data slope
[h p_P0 ci stats] = ttest(dummy_P0, data_P0); % compare dummy P0 vs original data P0
[h p_delP ci stats] = ttest(dummy_delP, data_delP); % compare dummy delP vs original data delP
[h p_T0 ci stats] = ttest(dummy_T0, data_T0); % compare dummy T0 vs original data T0
[h p_delT ci stats] = ttest(dummy_delT, data_delT); % compare dummy delT vs original data delT

dummy_vars = [dummy_best_indices' dummy_slopes']
data_vars = [data_best_indices data_slope]

p_slope
p_P0
p_delP
p_T0
p_delT

%% Plot
close all
figure(1)
hist(dummy_slopes), hold on, line([data_slope data_slope], [0 20], 'color', 'r', 'LineStyle', '--', 'LineWidth', 2)
% axis([-0.4 1 0 20])
title('Slope')

figure(2)
subplot(2, 2, 1)
hist(dummy_delP), hold on, line([data_delP data_delP], [0 20], 'color', 'r', 'LineStyle', '--', 'LineWidth', 2)
title('Delta P')

subplot(2, 2, 2)
hist(dummy_delT), hold on, line([data_delT data_delT], [0 20], 'color', 'r', 'LineStyle', '--', 'LineWidth', 2)
title('Delta T')

subplot(2, 2, 3)
hist(dummy_P0), hold on, line([data_P0 data_P0], [0 20], 'color', 'r', 'LineStyle', '--', 'LineWidth', 2)
title('P0')

subplot(2, 2, 4)
hist(dummy_T0), hold on, line([data_T0 data_T0], [0 20], 'color', 'r', 'LineStyle', '--', 'LineWidth', 2)
title('T0')

tilefigs
% [r p] = corr(dummy_slopes', sum_eigenvectors') % correlate dummy slopes and # 'big' eigenweights