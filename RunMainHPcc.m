function RunMainHPcc
%% function 'RunMain'
% Runs Main function to create LL, BIC, AIC HM model types
% Date of creation: 5/18/2010
% Author: Tony Vladusich, Brandeis University

%% Set parameters
par.basedirname = './HMM'; % base directory name
par.datasetdirname = '/Sadacca'; % dataset directory name
par.resultdirname = '/Results/OptiDummy7'; % results directory name
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
par.its = 500; % max number iterations of Baum-Welch algorithm per model cycle
par.max_HMM_states = 1:7; % max number of HM states per cycle (scalar or vector)
par.n_HM_models = length(par.max_HMM_states); % number of HM models
par.n_HM_cycles = 25; % number of cycles
par.n_workers = 1; % # parallel workers

%% Call Main function
for i =1:par.n_datasets
    
    par.datafilename = getfield(par.datafiles(i), 'name'); % get name
    
        for j = 1:par.n_tastants
            cd ./matlab
            par.tastant = j;
            filename='currentHMMdata.mat';
            save(filename,'par')
	     try
            	system('ssh -x login-node-1-0 "qsub -l neuro /home/sadacca/matlab/HPccHMM.job"');
	     catch ME
	     ME
            exit
	     end

            while exist(['/matlab/',filename])>0
                pause(20)
            end

            cd ..
            
        end
        
end
exit