function lwdataset = icaComputation (subName, group)

    
    global Cfg Log Paths
    
    condNames   = Cfg.condNames;
    nIC         = Cfg.Preprocessing.ICA.nIC;
    
    fprintf ('  -> ICA \n')
    
    
%%%%%%%%%%%%%%%%%%%%%    
% Create LW dataset % 
%%%%%%%%%%%%%%%%%%%%%

    cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Rereference')); 
    % lwDir   = dir (['ep_',condNames{iCond},'*',subName,'*.lw6']);
    % lwDirTemp  = dir (['reref*',subName,'.lw6']);
    
    for iCond  = 1:length(condNames)
        
        lwDir                  = dir (['reref*',condNames{iCond},'*',subName,'*.lw6']);
        if length(lwDir)>1

            % Find the latest file in the LW directory
            for i = 1: length(lwDir)
                time(i,1) = lwDir(i).datenum;            
            end
            
            [~,timeIdx] = max(time);
            lwDir       = lwDir(timeIdx,:);    
            fprintf ('  -> Last rereferenced version selected \n') 
        end
        
        option.filename{iCond} = char({fullfile(lwDir.folder,lwDir.name)});
    end
    
    lwdataset   = FLW_load.get_lwdataset(option);

%%%%%%%%%%%%%%%%%%%    
% ICA Computation %  
%%%%%%%%%%%%%%%%%%%

    movefile(Paths.lw6,Paths.lw6OutsideMatlabPath)

    cd(fullfile(Paths.LW,group,'Preprocessing/EEG/ICA')); 
    
    option      = struct('ICA_mode',2, ...
                         'algorithm',1, ...
                         'num_ICs',nIC, ...
                         'suffix','ica_merged','is_save',1);
    lwdataset   = FLW_compute_ICA_merged.get_lwdataset(lwdataset,option);
    
    movefile(Paths.lw6OutsideMatlabPath,Paths.lw6)
    
    
%%%%%%%%%%%%%%%%%%%
% Update Log file %  
%%%%%%%%%%%%%%%%%%%

    Log.(group).(subName).Preprocessing.ICA.nIC   =  Cfg.Preprocessing.ICA.nIC;

    fprintf ('  ==> %i ICs calculated \n',nIC)

end

