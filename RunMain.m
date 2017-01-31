function RunMain
%% function 'RunMain'
% Runs Main function to create LL, BIC, AIC HM model types


%% Set parameters
par.basedirname = 'C:\Users\sadaccabf\NIDA_DATA\OFC_SPC'; % base directory name
par.datasetdirname = '\HMM'; % dataset directory name
par.resultdirname = '\Results\newresults1ms'; % results directory name
par.simtype = '\AllData'; % all data or N-fold cross-validation
par.resultdir = [par.basedirname par.datasetdirname par.resultdirname par.simtype]; % results dir path
par.spikedatafiles = '\SpikeData1\*.mat'; % data files
par.datafiles = dir([par.basedirname par.datasetdirname par.spikedatafiles]); % datafile names
par.n_datasets = length(par.datafiles); % # datasets
par.n_tastants = 4; % # tastants
par.bintime = 10; % choose temporal binwidth
par.trialstart = 1; % start trial (msecs)
par.trialduration = 10000; % duration trial (msecs)
par.thresh = 0; % defines > 0 values of HM parameters
par.its = 200; % max number iterations of Baum-Welch algorithm per model cycle
par.max_HMM_states = 1:5; % max number of HM states per cycle (scalar or vector)
par.n_HM_models = length(par.max_HMM_states); % number of HM models
par.n_HM_cycles = 25; % number of cycles
par.n_workers = 1; % # parallel workers

                                   

%% Call Main function
for i = 1:par.n_datasets
    
    par.datafilename = getfield(par.datafiles(i), 'name'); % get name
    
        for j = 1:par.n_tastants
            
            par.tastant = j;
            Main(par);
            
        end
        
end