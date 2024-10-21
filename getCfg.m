function [Cfg, Paths, Log] = getCfg(projectPath)

% All parameters for the preprocessing, analysis and plotting of this
% project

% Init parameter structure
Paths   = struct(); 
Cfg     = struct(); 
Log     = struct();


%% Paths & Log

    Paths.projectName               = 'NonRepComparison';
    
    Paths.lw7                       = fullfile(projectPath,'lib/letswave7-master');
    Paths.lw6                       = fullfile(projectPath,'lib/letswave6-master');
    Paths.lw6ElecLoc                = fullfile(Paths.lw6,'core_functions/Standard-10-20-Cap81.locs');
    Paths.lw6OutsideMatlabPath      = fullfile(extractBefore(projectPath,"PROJECTS"),'letswave6-master');
    addpath (genpath (Paths.lw7)); 

    %Paths.toolboxes                 = fullfile(extractBefore(projectPath,"PROJECTS"),'Toolboxes');
    %Paths.plottingFunctions         = fullfile(extractBefore(projectPath,"PROJECTS"),'Toolboxes/plottingFunctions');
%     addpath (genpath (Paths.toolboxes)); 
%     addpath (genpath (Paths.plottingFunctions));
    
    Paths.project                   = projectPath;
    
    Paths.dataDir                   = {fullfile(Paths.project,'Datasets/Raw/Musicians');...
                                       fullfile(Paths.project,'Datasets/Raw/NonMusicians')};
    Paths.LW                        = fullfile(Paths.project,'Datasets/LW');
    Paths.preprocessing             = fullfile(Paths.project,'2_preprocessing');
    Paths.analysis                  = fullfile(Paths.project,'3_analysis');
    Paths.figures                   = fullfile(Paths.project,'4_figures');
                  
    Paths.Preprocessing.elecLabels  = fullfile(Paths.preprocessing,'electrodelabels.csv');
    
    load(fullfile(Paths.project, 'Log.mat'));
%% Experiment's setting

    % condNames & groups
    Cfg.condNames   = {'ABAB1', 'AAAA', 'ABAB2', 'CDEF'};
    Cfg.groups      = {'Musicians','NonMusicians'};
    Cfg.trialDur    = 67.2;
    
%% EEG preprocessing

    % Butterworth filter
    Cfg.Preprocessing.HPfilter.filterType   = 'highpass';
    Cfg.Preprocessing.HPfilter.cutOff       = 0.1;
    Cfg.Preprocessing.HPfilter.Order        = 4;
    
    % Electrode labels
    option                                  = delimitedTextImportOptions ("NumVariables", 1);
    option.DataLines                        = [1, Inf]; % Get all the rows containing the electrode labels
    Cfg.Preprocessing.elecLabels.labels     = readtable (Paths.Preprocessing.elecLabels, option);       
    
    % Rereference
    Cfg.Preprocessing.Reref.Elec            = 'Mastoids'; % {'Mastoids','Mast Right', 'Mast Left', 'All'}
    Cfg.Preprocessing.multiRerefenceMethodNames ...
                                            = {'rerefMast';'rerefAllInclMast';'rerefAllExclMast'};
    
    % ICA
    Cfg.Preprocessing.ICA.nIC              = 60;
    
    % Average
    Cfg.Preprocessing.NotAverage.lwSave    = false;

%% EEG Analysis

    Cfg.zscore.ampBasedZscoreBoundaries    = [1, 30];  % in Hz
    Cfg.zscore.iElec                       = 68;       % only compute the fronto-central pool of electrodes
    
    % FFT
    Cfg.Analysis.FFT.lwSave                = 0;
    
    % BlSNR
    Cfg.Analysis.Bl_snr.lowBin             = 2;
    Cfg.Analysis.Bl_snr.highBin            = 5;
    Cfg.Analysis.Bl_snr.lwSave             = 1;
    
    % Findpeak Vs zSNR
    Cfg.Analysis.PeakDetectionMethod       = 'findpeaks'; %'zSNR';
    Cfg.Analysis.PeakDetectionNSideBins    = [1 5]; % si on prend le 2 à 5, souvent il sélectionne les pentes des pics puisque il n'y a pas les pics adjacents ne sont pas pris en compte. La résolution spectrale est suffisamment bonne et il ne semble pas y avoir de spectral leakage.
    Cfg.Analysis.lowerLimit                = 1.23; 
    Cfg.Analysis.upperLimit                = 30;
    Cfg.Analysis.zSNRSign                  = 0.00001; % 0.0001
    
    Cfg.Analysis.Stats.percSign            = 0.05;
    Cfg.Analysis.Stats.sign                = norminv(1 - Cfg.Analysis.Stats.percSign); % p-value -> zscore

    

    %% Tapping preprocessing

    Cfg.Preprocessing.Tapping.Bl_snr.lowBin     = 2;
    Cfg.Preprocessing.Tapping.Bl_snr.highBin    = 5;

    Paths.('Musicians').TappingTaps             = fullfile(Paths.project, 'datasets/LW/Musicians/Preprocessing/Tapping/CheckTaps/TappingTaps.mat');
    Paths.('NonMusicians').TappingTaps          = fullfile(Paths.project, 'datasets/LW/NonMusicians/Preprocessing/Tapping/CheckTaps/TappingTaps.mat');
 
    Paths.('Musicians').PreprocessedData        = fullfile(Paths.project, 'datasets/LW/Musicians/Preprocessing/Tapping/PreprocessedData.mat');
    Paths.('NonMusicians').PreprocessedData     = fullfile(Paths.project, 'datasets/LW/NonMusicians/Preprocessing/Tapping/PreprocessedData.mat');


    %% Taping Analysis

    Paths.('Musicians').ITI                      = fullfile(Paths.project, 'datasets/LW/Musicians/Analysis/Tapping/ITI/ITI.mat');
    Paths.('NonMusicians').ITI                   = fullfile(Paths.project, 'datasets/LW/NonMusicians/Analysis/Tapping/ITI/ITI.mat');

    Paths.('Musicians').TappingAnalysis          = fullfile(Paths.project, 'datasets/LW/Musicians/Analysis/Tapping/TappingAnalysis.mat');
    Paths.('NonMusicians').TappingAnalysis       = fullfile(Paths.project, 'datasets/LW/NonMusicians/Analysis/Tapping/TappingAnalysis.mat');

    Paths.('Musicians').Raw                      = fullfile(Paths.project, 'datasets/Raw/Musicians');
    Paths.('NonMusicians').Raw                   = fullfile(Paths.project, 'datasets/Raw/NonMusicians');

    Cfg.gridIOI = 0.2; % stim.gridIOI;
    warning('GridIOI of 0.2 is currently hardcoded in getCfg')
    
%% Plotting

    Cfg.Figure.fontsize     = 15;
    Cfg.Figure.linewidth    = 1.5;
    Cfg.Figure.resolution   = '-r600';
    Cfg.Figure.plotNames    = {'Rhythm ABAB1','Rhythm AAAA','Rhythm ABAB2','Rhythm CDEF'};
    
    % Colors
    Cfg.Figure.lightGrey       = [0.6 0.6 0.6];
    Cfg.Figure.metRelColor     = [238 44  44 ]/255;
    Cfg.Figure.metUnrelColor   = [24  116 205]/255;
    Cfg.Figure.orange          = [251 133 0  ]/255;
    Cfg.Figure.lightOrange     = [255 198 133]/255;
    Cfg.Figure.purple          = [83  55  71 ]/255;
    Cfg.Figure.lightPurple     = [198 169 186]/255;
    Cfg.Figure.green           = [105 153 93 ]/255;
    Cfg.Figure.lightGreen      = [185 208 179]/255;
    
    
    % axe position
    Cfg.Figure.fftPos_individual   = {[0.32 0.64 0.18 0.24]; [0.32 0.22 0.18 0.24]; [0.55 0.64 0.18 0.24];  [0.55 0.22 0.18 0.24]};
    Cfg.Figure.zPos_individual     = {[0.08 0.64 0.18 0.26]; [0.08 0.22 0.18 0.26]; [0.8  0.64 0.18 0.26];  [0.8  0.22 0.18 0.26]};
    Cfg.Figure.condTitle           = {[0.08 0.88 0.42 0.05]; [0.08 0.48 0.42 0.05]; [0.55 0.88 0.42 0.05];  [0.55 0.48 0.42 0.05]};
    
    Cfg.Figure.zScoreMetRelLabels  = {[0.0783,0.5446,0.0894,0.087]; [0.0783,0.0923,0.0894,0.087]; [0.7975,0.5446,0.0894,0.087]; [0.7975,0.0923,0.0894,0.087]};
    Cfg.Figure.zScoreMetUnrelLabels= {[0.1783,0.5446,0.0894,0.087]; [0.1783,0.0923,0.0894,0.087]; [0.8991,0.5446,0.0894,0.087]; [0.8991,0.0923,0.0894,0.087]};

%% Stim Analysis

    Cfg.Stim.gridIOI = 0.2;
    
    Cfg.Models.name = {'Auditory Nerve','Inferior Colliculus','Envelope','None'};

end

