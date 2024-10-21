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
  
    global Cfg Log Paths
    
    [Cfg, Paths, Log] = getCfg(projectPath);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    condNames = Cfg.condNames;
    groups    = Cfg.groups;
    
    [Groups2test,~] = listdlg('PromptString',{'Select the group of participant you want to preprocess'},...
                              'ListString',groups);


%% load & change tapping triggers

for iGroups = Groups2test
    
    % get subject directory
    cd(Paths.dataDir{iGroups})
    subDir = dir(fullfile(Paths.dataDir{iGroups},'sub*.lw6'));
    
    subjects = {subDir(1:end).name};
    
    [subjects2test,~] = listdlg('PromptString',{'Select the participants you want to preprocess'},...
                              'ListString',subjects);    
    
    for iSubjects = subjects2test
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % load participant in LW %
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,subName,~] = fileparts(subDir(iSubjects).name);
        
        % experiment design and condition order
        Log.(groups{iGroups}).(subName).Design      = struct2dataset (bids.util.tsvread(fullfile(Paths.dataDir{iGroups},[subName,'.tsv'])));
        Log.(groups{iGroups}).(subName).condOrder   = unique(Log.(groups{iGroups}).(subName).Design.condName,'stable');
        
        fprintf('------------------------------------------------- \n  -> %s to select \n',subName) 
        fprintf ('  -> Open the continuous data viewer, and click on the "send events \n     to workspace" button (right below the Events table) \n')    
        letswave7  
        
        %%%%%%%%%%%%%%%%%%%
        % Change triggers %
        %%%%%%%%%%%%%%%%%%%
        
        pause(10)   
        while ~exist('events','var')

            % etime(datevec(now),datevec(newFile.datenum)) < 30
            
            % Terminate if no "events" variable is found in the workspace after 60seconds
            iniTime = clock;
            if etime (clock, iniTime) > 60
                error ('Runtime error : Events were not correctly loaded in the workspace...')
            end
            disp('Come on, Chop, Chop !!!')
            pause(5)
        end; clear iniTime     
        
        % extract indexes
        listenIdx   = find(strcmp({events.code}, '1')==1);
        tapIdx      = find(strcmp({events.code}, '2')==1); 
        allIdx      = sort([listenIdx, tapIdx]);
        
        j=1;
        for iIdx = 1:length(allIdx)                   
            
            if Log.(groups{iGroups}).(subName).Design(j,:).terminatedTrial          % terminated trial  
                disp('terminated trial')
                

            elseif iIdx < length(Log.(groups{iGroups}).(subName).Design) && ...
                Log.(groups{iGroups}).(subName).Design(j+1,:).repeatedPrevTrial     % repeated trial
                disp('repeated trial')

                
            elseif strcmp(Log.(groups{iGroups}).(subName).Design(j,:).trial_type,'tap') % tapping trial
                for iCond = 1:length(condNames)
                    if strcmp(Log.(groups{iGroups}).(subName).Design(j,:).condName,condNames{iCond})
                        events(allIdx(iIdx)).code = [condNames{iCond},'_tap'];
                    end
                end
                
            elseif strcmp(Log.(groups{iGroups}).(subName).Design(j,:).trial_type,'listen') % listening trial
                for iCond = 1:length(condNames)
                    if strcmp(Log.(groups{iGroups}).(subName).Design(j,:).condName,condNames{iCond})
                        events(allIdx(iIdx)).code = [condNames{iCond},'_listen'];
                    end
                end
            end
            
            j=j+1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % load updated events %
        %%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf ('  -> Load the new events in Letswave (button next to the previous one) \n     and save it (without overwriting) \n')    

        loaded_events = 'false';
        
        while strcmp(loaded_events,'false') 
            newFile = dir(['ep_edit ',subName,'.lw6']);
            if ~isempty(newFile) && etime(datevec(now),datevec(newFile.datenum)) < 30 % new ep_edit less than 30seconds before
                loaded_events = 'true';
            else
                pause(5)
            end
        end
        
        pause(3) % time for LW to save the ep_edit
        close all force
        clear events
        

    end
    
    [~, neworder]         = sort(lower(fieldnames(Log.(groups{iGroups}))));
    Log.(groups{iGroups}) = orderfields(Log.(groups{iGroups}), neworder);    
    
end

cd(Paths.project)
save('Log.mat','Log')


