% Add electrodes' location for topoplots %

% This script is only to use in the musician group from participant 007 to 
% 031 and in the non-musicians group until participant 22


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
    
    
% parameters 
    groups = Cfg.groups;   
                          
                          
   

%% Electrode locations 

for iGroups = 2% 1:length(groups)
    
    % get subject directory
    cd(Paths.dataDir{iGroups})
    subDir      = dir(fullfile(Paths.dataDir{iGroups},'ep_edit*.lw6'));
    subjects    = extractBetween({subDir(1:end).name},'ep_edit ','.lw6');
       
    [subjects2test,~] = listdlg('PromptString',{'Select the participants you want to preprocess'},...
                              'ListString',subjects);      
    for iSubjects = subjects2test
        
        % go through each folder...
        cd(Paths.LW)
        error('go through each folder...')
        
        
        
        option      = struct('filename',fullfile(ICAFolder,ICAFileName));
        lwdata      = FLW_load.get_lwdata(option); 

        fprintf ('  -> Add electrode location \n')

        addpath (genpath (Paths.lw6)); 

        locs = readlocs(Paths.lw6ElecLoc);

        for iLocs = 1:length(locs)
            locs(iLocs).topo_enabled = 1;
        end

        [lwdata.header] = RLW_edit_electrodes(lwdata.header,locs);

        rmpath(Paths.lw6);

        fprintf ('  ==> electrode location added \n')
    
    end
end
