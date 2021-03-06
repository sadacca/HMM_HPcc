function par = ReadPlotAnalFile(ParFileName, PlottingToggle, ModelToggle)
%% function 'par = ReadPlotAnalFile(ParFileName, PlottingToggle, ModelToggle)'
% Loads and plots decoded HM sequences for selected model type (LL, AIC,
% BIC), along with original spiking observations

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Load data
load(ParFileName, 'par')

%% Call model
switch ModelToggle
    case('LL')
        par.hmm_sequence = par.hmm_LL_sequence;
        par.HMM_transition_matrix = par.LL_transition_matrix;
        par.HMM_emission_matrix = par.LL_emission_matrix;
    case('AIC')
        par.hmm_sequence = par.hmm_AIC_sequence;
        par.HMM_transition_matrix = par.AIC_transition_matrix;
        par.HMM_emission_matrix = par.AIC_emission_matrix;        
    case('BIC')
        par.hmm_sequence = par.hmm_BIC_sequence;
        par.HMM_transition_matrix = par.BIC_transition_matrix;
        par.HMM_emission_matrix = par.BIC_emission_matrix;
end

%% Plot figures
close all

switch PlottingToggle
    case('on')
        
        for j = 1:par.n_trials
            
            figure(j), clf
            plot(par.norm_observations(j, :) + 1.1, 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 3)
            hold on
            plot(par.hmm_sequence(:, :, j)', 'LineWidth', 2) % plot HM-estimated state sequence
            axis([0 par.n_observations + 1 0  2.2]);
            set(gca, 'YTick', [0 1])
            set(gca, 'YTickLabel', [0 1])
            addylabel('Probability', 16, [0.075 0.28]);
            addxlabel('Time bin', 16, [0.5 0.05]);
            addylabel('Spikes', 16, [0.075 0.7]);
            
        end
        
        figure(j + 1), clf
        bar(par.HMM_transition_matrix, 'stacked')
        axis([0 size(par.HMM_transition_matrix, 1) + 1 0  1]);
        title('HMM transition', 'FontSize', 16) 
        addylabel('Probability', 16, [0.05 0.5]);
        addxlabel('Parameter (grouped by row)', 16, [0.5 0.02]);
        
        figure(j + 2), clf
        bar(normalise(par.HMM_emission_matrix(:, 2:end), 2), 'stacked')
        axis([0 size(par.HMM_emission_matrix, 1) + 1 0  1]);
        title('HMM emission', 'FontSize', 16)
        
        addylabel('Probability', 16, [0.05 0.5]);
        addxlabel('Parameter (grouped by row)', 16, [0.5 0.02]);
           
        tilefigs
        
end
