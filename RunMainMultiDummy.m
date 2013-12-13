function RunMainMultiDummy
%% function 'RunMainMultiDummy'
% Runs Main function to create LL, BIC, AIC HM model types: Generates
% multiple dummy datasets for analysis

% Date of creation of main file: 5/18/2010
% Author of main file: Tony Vladusich, Brandeis University

% Author of modified file: Brian Sadacca, Brandeis University



%% Number of dummy sets to run
number_dummy_sets = 20;

%% Set-up Matlab parallel 'workers'
if matlabpool('size') > 0
    matlabpool close % close any existing workers
    matlabpool % now open multiple workers
else
    matlabpool % open multiple workers
end

for dummy_set_number = 19:number_dummy_sets
    %% Set parameters
    par.basedirname = '/Users/TheVlad/Documents/MATLAB/work/HMM'; % base directory name
    par.datasetdirname = '/Sadacca'; % dataset directory name
    par.resultdirname = '/Results'; % results directory name
    par.simtype = '/AllDummyData'; % all data or N-fold cross-validation
    par.dummy_set_number_str = num2str(dummy_set_number); % dummy set number
    resultdir = [par.basedirname par.datasetdirname par.resultdirname par.simtype par.dummy_set_number_str]; % results dir path
    mkdir(resultdir)
    par.resultdir = resultdir; % results dir path
    par.sourcedirprefix = '/DummySpikeData'; % src dir prefix
    par.spikedatafiles = [par.sourcedirprefix par.dummy_set_number_str '/*.mat']; % data files
    par.datafiles = dir([par.basedirname par.datasetdirname par.spikedatafiles]); % datafile names
    par.sourcedir = [par.basedirname par.datasetdirname par.sourcedirprefix par.dummy_set_number_str]; % src dir
    cd(par.sourcedir); % cd to source dir
    par.n_datasets = length(par.datafiles); % # datasets
    par.n_tastants = 6; % # tastants
    par.bintime = 10; % choose temporal binwidth
    par.trialstart = 1001; % start trial (msecs)
    par.trialduration = 2500; % duration trial (msecs)
    par.thresh = 0; % defines > 0 values of HM parameters
    par.its = 200; % max number iterations of Baum-Welch algorithm per model cycle
    par.max_HMM_states = 1:7; % max number of HM states per cycle (scalar or vector)
    par.n_HM_models = length(par.max_HMM_states); % number of HM models
    par.n_HM_cycles = 50; % number of cycles
    par.n_workers = 2; % # parallel workers
    
    %% Call Main function
    for i = 1:par.n_datasets
        
        par.datafilename = getfield(par.datafiles(i), 'name'); % get name
        for j = 1:par.n_tastants
            
            par.tastant = j;
            Main(par);
            cd(par.sourcedir) % go to results directory
            
        end
        
    end
    
end