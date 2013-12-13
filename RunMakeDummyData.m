function RunMakeDummyData
%% function 'RunMakeDummyData'
% Runs function to create new sets of 'dummy' spike trains based on PSTH
% data calculated at X ms resolution

% Date of creation: 2/15/2011
% Author: Tony Vladusich, Brandeis University

%% Set parameters
par.basedirname = '/Users/TheVlad/Documents/MATLAB/work/HMM'; % base directory name
par.datasetdirname = '/Sadacca'; % dataset directory name
par.resultdirname = '/DummySpikeData'; % results directory name
par.resultdir = [par.basedirname par.datasetdirname par.resultdirname]; % results dir path
par.spikedatafiles = '/SpikeData/*.mat'; % data files
par.datafiles = dir([par.basedirname par.datasetdirname par.spikedatafiles]); % datafile names
par.n_datasets = length(par.datafiles); % # datasets
par.n_tastants = 6; % # tastants
par.bintime = 1; % choose temporal binwidth
par.trialstart = 1; % start trial (msecs)
par.trialduration = 5000; % duration trial (msecs)
par.PlottingToggle = 'off'; % toggle plot trial figs
par.max_n_symbols = 15; % max. # symbols in all datasets
par.smooth_factor = 1; % PSTH smoothing filter width (1 = no smoothing)

%% Call MakeDummyData function

for i = 1:par.n_datasets
    
    par.datafilename = getfield(par.datafiles(i), 'name'); % get name
    
    datafilename = par.datafilename
        
    for j = 1:par.n_tastants
        
        tastant = j
        
        par.tastant = j;
        labelled_spike_array = MakeDummyData(par);
        dat(:, par.tastant, :, :) = labelled_spike_array;

    end
    
    cd([par.resultdir]) % go to results directory
    save(par.datafilename, 'dat') % save
    clear dat
        
end