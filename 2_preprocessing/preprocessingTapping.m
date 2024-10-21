%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% NonRepComplex %%%%%%%%%%%%%%%
%%%%%% Data Preprocessing - Tapping %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths + Preprocessing parameters

%%%%%%%%%%%%%%%%%%
% Initial setups %
%%%%%%%%%%%%%%%%%%


clear all; clc

% Set path and load the Cfg and Paths files
[projectPath,~,~]   = fileparts(mfilename('fullpath'));
projectPath         = extractBefore(projectPath,"Preprocessing");
addpath (genpath (projectPath));

[Cfg, Paths, Log] = getCfg(projectPath);

LW_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condNames = Cfg.condNames;
groups    = Cfg.groups;

[Groups2test,~] = listdlg('PromptString',{'Select the group of participant you want to preprocess'},...
    'ListString',groups);

%% Tapping preprocessing


for iGroups = Groups2test
    group          = groups{iGroups};

    if exist(Paths.(group).PreprocessedData)
        Taps       = load(Paths.(group).PreprocessedData).Taps;
        TapAnomaly = load(Paths.(group).PreprocessedData).TapAnomaly;
    elseif exist(Paths.(group).TappingTaps)
        Taps       = load(Paths.(group).TappingTaps).TappingTaps;
        TapAnomaly = load(Paths.(group).TappingTaps).TapAnomaly;
    else
        Taps       = struct();
        TapAnomaly = struct();
    end

    % get subject directory
    %cd(Paths.dataDir{iGroups})
    subDir      = dir(fullfile(Paths.dataDir{iGroups},'ep_edit*.lw6'));
    subjects    = {subDir(1:end).name};

    [subjects2test,~] = listdlg('PromptString',{'Select the participants you want to preprocess'},...
        'ListString',subjects);

    for iSubjects = subjects2test

        % load LW file
        subName     = subjects{iSubjects};
        option      = struct('filename',subName);
        [lwdata.header, lwdata.data]      = CLW_load(fullfile(Paths.(group).Raw,subName));

        subName     = strrep(strrep(subName,'ep_edit ',''),'.lw6','');

        fprintf("\n-----------------------------\n Start of the preprocessing for participant "+subName+" \n\n" )

        % Removal of unused electrodes %

        [lwdata, Log] = tapRemovalUnusedElectrodes(lwdata, group, subName, Cfg, Log, Paths);

        % Segmentation
        [lwdataset, Log] = tapEpochSegmentation(lwdata, group, subName, Cfg, Log, Paths);

        % FFT
        [lwdataset, lwdatasetFFT, Log] = tapFFT(lwdataset, group, subName,  Cfg, Log, Paths);

        % BlSNR
        [lwdatasetFFT, lwdatasetBlSNR, Log] = tapBlSNR_FFT(lwdatasetFFT, group, subName, Cfg, Log, Paths);

        % CheckTap
        [lwdataset, Taps, TapAnomaly, Log] = tapCheckTap(lwdataset, group, subName, Taps, TapAnomaly, Cfg, Log, Paths);

        % saving preprocess data
        if exist(Paths.(group).PreprocessedData)
            disp("Saving Data");
            save(Paths.(group).PreprocessedData, 'Taps', 'TapAnomaly', 'Cfg', 'Paths', 'Log', '-append')
        else
            disp("Saving Data");
            save(Paths.(group).PreprocessedData, 'Taps', 'TapAnomaly', 'Cfg', 'Paths', 'Log')

        end
        disp("Participant:"+subName+" finished");
    end
end


