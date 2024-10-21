%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% NonRepComplex %%%%%%%%%%%%%%%
%%%%%%%% Data Preprocessing - EEG %%%%%%%%
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
  
    global Cfg Log Paths

    [Cfg, Paths, Log] = getCfg(projectPath); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    condNames = Cfg.condNames;
    groups    = Cfg.groups;
    
    [Groups2test,~] = listdlg('PromptString',{'Select the group of participant you want to preprocess'},...
                              'ListString',groups);
                          
%% EEG preprocessing  


for iGroups = Groups2test    
    
    % get subject directory
    cd(Paths.dataDir{iGroups})
    subDir      = dir(fullfile(Paths.dataDir{iGroups},'ep_edit*.lw6'));
    subjects    = {subDir(1:end).name};
    
    [subjects2test,~] = listdlg('PromptString',{'Select the participants you want to preprocess'},...
                              'ListString',subjects);     
        
    for iSubjects = subjects2test  
        
        % load LW file
        subName     = subjects{iSubjects};
        option      = struct('filename',subName);
        lwdata      = FLW_load.get_lwdata (option);
         
        subName     = strrep(strrep(subName,'ep_edit ',''),'.lw6','');
        group       = groups{iGroups};
               
        fprintf('\n-----------------------------\n Start of the preprocessing \n     %s selected \n\n',subName)

        % High-Pass Filter
        lwdata = hpFilter(lwdata, subName, group);

        % New electrode labels and removal of unused electrodes
        lwdata = electrodes(lwdata, subName, group);

        % Interpolation
        lwdata = interpolation(lwdata, subName, group);

        % Epoch Segmentation
        lwdata = epochSegmentation(lwdata, subName, group);

        % Rereferencing
        rereference (subName, group);

        % ICA
        lwdataset = icaComputation (subName, group);

        icsRemoval(subName, group)

        % Averaging and saving 
        averageAndSave(subName, group)        
        
    end
end
