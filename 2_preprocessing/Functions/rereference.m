function rereference (subName, group)

    global Cfg Log Paths
    
    condNames = Cfg.condNames;
    
    fprintf ('  -> Rereference \n')
    
    
%%%%%%%%%%%%%%%%%%%%%    
% Define electrodes % 
%%%%%%%%%%%%%%%%%%%%% 
    
    if strcmp(Cfg.Preprocessing.Reref.Elec,'Mastoids')
        refElec   = {'Mast1','Mast2'};
        
    elseif strcmp(Cfg.Preprocessing.Reref.Elec,'Mast Right')  
        refElec   = {'Mast2'};
        
    elseif strcmp(Cfg.Preprocessing.Reref.Elec,'Mast Left')   
        refElec   = {'Mast1'};
        
    else
        refElec   = table2cell(Cfg.Preprocessing.elecLabels.labels);
    end
    
    applyElec     = table2cell(Cfg.Preprocessing.elecLabels.labels);
    
    % If we don't want to apply it on the mastoids:
    % applyElec     = setdiff(applyElec,{'Mast1','Mast2'}) 
    
    
%%%%%%%%%%%%%%%    
% Rereference %
%%%%%%%%%%%%%%% 

    for iCond = 1:size(condNames,2)
        
        % Load each segmented lwdata
        cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Segmentation'));
        
        lwDir   = dir (['ep_',condNames{iCond},'*',subName,'*.lw6']);
        
        if length(lwDir)>1   
            % Find the latest file in the LW directory
            for i = 1: length(lwDir)
                time(i,1) = lwDir(i).datenum;            
            end
            
            [~,timeIdx] = max(time);
            lwDir       = lwDir(timeIdx,:);    
            fprintf ('  -> Last segmented version selected \n') 
        end
                       
        option  = struct ('filename',fullfile(lwDir.folder,lwDir.name));
        lwdata  = FLW_load.get_lwdata (option);
        
        cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Rereference'));      

        %lw7
%         option  = struct('reference_list',{refElec}, ...
%                          'apply_list',{applyElec}, ...
%                          'suffix','rerefLW7','is_save',1);
%         lwdata  = FLW_rereference.get_lwdata(lwdata,option);         
        
        % lw6
        [lwdataReref.header, lwdataReref.data] = RLW_rereference(lwdata.header, lwdata.data, 'apply_list',applyElec, 'reference_list', refElec);
        lwdataReref.header.name = ['reref ',lwdata.header.name];
        CLW_save (fullfile(Paths.LW,group,'Preprocessing/EEG/Rereference'), ...
                  lwdataReref.header, lwdataReref.data);
        
     end 
    
%%%%%%%%%%%%%%%%%%%
% Update Log file %  
%%%%%%%%%%%%%%%%%%%

    Log.(group).(subName).Preprocessing.Reref.refElec   =  Cfg.Preprocessing.Reref.Elec;
    Log.(group).(subName).Preprocessing.Reref.applyElec =  applyElec;
    
    fprintf ('  ==> Signal rereferenced \n')
    
end

