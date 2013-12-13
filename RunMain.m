function RunMain
%% function 'RunMain'
% Runs Main function to create LL, BIC, AIC HM model types

% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Set parameters
par.basedirname = './HMM'; % base directory name
par.datasetdirname = '/Sadacca'; % dataset directory name
par.resultdirname = '/Results/newresults1ms'; % results directory name
par.simtype = '/AllData'; % all data or N-fold cross-validation
par.resultdir = [par.basedirname par.datasetdirname par.resultdirname par.simtype]; % results dir path
par.spikedatafiles = '/SpikeData/*.mat'; % data files
par.datafiles = dir([par.basedirname par.datasetdirname par.spikedatafiles]); % datafile names
par.n_datasets = length(par.datafiles); % # datasets
par.n_tastants = 6; % # tastants
par.bintime = 10; % choose temporal binwidth
par.trialstart = 1001; % start trial (msecs)
par.trialduration = 2500; % duration trial (msecs)
par.thresh = 0; % defines > 0 values of HM parameters
par.its = 400; % max number iterations of Baum-Welch algorithm per model cycle
par.max_HMM_states = 1:7; % max number of HM states per cycle (scalar or vector)
par.n_HM_models = length(par.max_HMM_states); % number of HM models
par.n_HM_cycles = 50; % number of cycles
par.n_workers = 1; % # parallel workers

%% Set-up Matlab parallel 'workers'
% if matlabpool('size') > 0
%     matlabpool close % close any existing workers
%     matlabpool % now open multiple workers
% else
%     matlabpool % open multiple workers
% end

%% Call Main function
for i = 1:par.n_datasets
    
    par.datafilename = getfield(par.datafiles(i), 'name'); % get name
    
        for j = 1:par.n_tastants
            
            par.tastant = j;
            Main(par);
            
        end
        
end