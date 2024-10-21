% Multirereferencing %



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
    
    if isfolder(Paths.lw6)
        movefile(Paths.lw6,Paths.lw6OutsideMatlabPath)
    end
   
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    condNames   = Cfg.condNames;
    groups      = Cfg.groups;
    
    % rereferencing methods
    methodNames{1}   = Cfg.Preprocessing.multiRerefenceMethodNames{1}; 
    refElecMeth{1}   = {'Mast1','Mast2'};
    applyElecMeth{1} = table2cell(Cfg.Preprocessing.elecLabels.labels);    

    methodNames{2}   = Cfg.Preprocessing.multiRerefenceMethodNames{2}; 
    refElecMeth{2}   = table2cell(Cfg.Preprocessing.elecLabels.labels);
    applyElecMeth{2} = table2cell(Cfg.Preprocessing.elecLabels.labels);  

    methodNames{3}   = Cfg.Preprocessing.multiRerefenceMethodNames{3}; 
    refElecMeth{3}   = table2cell(Cfg.Preprocessing.elecLabels.labels(1:end-2,:));
    applyElecMeth{3} = table2cell(Cfg.Preprocessing.elecLabels.labels(1:end-2,:)); 
    
    [Groups2test,~]  = listdlg('PromptString',{'Select the group of participant you want to preprocess'},...
                               'ListString',groups);

%% MultiRereference

    for iGroups = Groups2test  
        
        % get subject directory
        cd(Paths.dataDir{iGroups})
        subDir      = dir(fullfile(Paths.dataDir{iGroups},'ep_edit*.lw6'));
        subjects    = extractBetween({subDir(1:end).name},'ep_edit ','.lw6');


        [subjects2test,~] = listdlg('PromptString',{'Select the participants you want to preprocess'},...
                                  'ListString',subjects);   
        
        preprocessedSubjects{iGroups} = subjects2test;

        for iSubjects = subjects2test   
            
            group       = groups{iGroups};
            subName     = subjects{iSubjects};
            
            for iCond = 1:size(condNames,2)

                % Load each segmented lwdata
                cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Segmentation'));

                lwDir   = dir (['ep_',condNames{iCond},'*',subName,'.lw6']);

                if length(lwDir)>1   
                    % Find the latest file in the LW directory
                    for i = 1: length(lwDir)
                        time(i,1) = lwDir(i).datenum;            
                    end

                    [~,timeIdx] = max(time);
                    lwDir       = lwDir(timeIdx,:);    
                    fprintf ('  -> Last segmented version selected \n') 
                end            
            
                for iMethod = 1:length(methodNames)
                    option  = struct ('filename',fullfile(lwDir.folder,lwDir.name));
                    lwdata  = FLW_load.get_lwdata (option);

                    cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Rereference/MultiRereference'));  
                    
                    % rereferencing
                    suffix      = methodNames{iMethod};
                    refElec     = refElecMeth{iMethod};
                    applyElec   = applyElecMeth{iMethod};

                    option      = struct('reference_list',{refElec}, ...
                                         'apply_list',{applyElec}, ...
                                         'suffix',suffix,'is_save',1);
                    lwdata      = FLW_rereference.get_lwdata(lwdata,option); 
               
                end
            end
            
            %%%%%%%
            % ICA %
            %%%%%%%
            
            % create lwdataset
                        
            for iMethod = 1:length(methodNames)
                for iCond = 1: length(condNames)

                    cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Rereference/MultiRereference'));
                    
                    lwDir = dir ([methodNames{iMethod},'*',condNames{iCond},'*',subName,'.lw6']);
                    
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

                    cd(fullfile(Paths.LW,group,'Preprocessing/EEG/ICA/MultiRereference')); 
                    

                    option      = struct('ICA_mode',2, ...
                                         'algorithm',1, ...
                                         'num_ICs',Cfg.Preprocessing.ICA.nIC, ...
                                         'suffix','ica_merged','is_save',1);
                    lwdataset   = FLW_compute_ICA_merged.get_lwdataset(lwdataset,option);
            end
       end
    end
        
%% Remove unnecessary ICs manually

error('Stop the script to manually remove the ICs')

%% Average and save
for iMethod = 1:length(methodNames)
    for iGroups = Groups2test 
        for iSubjects = preprocessedSubjects{iGroups}
            for iCond = 1: length(condNames)

                group       = groups{iGroups};
                subName     = subjects{iSubjects};

                %%%%%%%%%%%%%%%%%%%%%%%
                % Load sp_filter data %
                %%%%%%%%%%%%%%%%%%%%%%%

                cd(fullfile(Paths.LW,group,'Preprocessing/EEG/ICA/MultiRereference'))
                subDir = dir(['sp_filter*',methodNames{iMethod},'*',(condNames{iCond}),'*',(subName),'.lw6']);  

                if length(subDir)>1
                    [ICAFileName, ICAFolder] = uigetfile(['sp_filter*',(condNames{iCond}),'*',(subName),'*lw6']);

                elseif isempty(subDir)
                    error('No sp_filter file found. Please look for it manually or select the rereference file if no ICA was computed')

                else     
                    ICAFileName    = subDir.name;
                    ICAFolder      = subDir.folder;
                end            

                option      = struct('filename',fullfile(ICAFolder,ICAFileName));
                lwdata      = FLW_load.get_lwdata(option); 


                %%%%%%%%%%%%%%%%%%%%%%
                % Average all trials %
                %%%%%%%%%%%%%%%%%%%%%%

                cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Average/MultiRereference'))

                option = struct('operation','average','suffix','avg','is_save',1);
                lwdata = FLW_average_epochs.get_lwdata(lwdata,option);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Store the averaged preprocessing % 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Musicians
                if strcmp(group,Cfg.groups{1})

                    if strcmp(methodNames{iMethod},'rerefMast')

                        PreprocessedEEG_Musicians_rerefMast.(condNames{iCond}).(subName) = lwdata;
    
                        % reorder the subjects
                        [~, neworder]   = sort(lower(fieldnames(PreprocessedEEG_Musicians_rerefMast.(condNames{iCond}))));
                        PreprocessedEEG_Musicians_rerefMast.(condNames{iCond}) = orderfields(PreprocessedEEG_Musicians_rerefMast.(condNames{iCond}), neworder);
                    
                    elseif strcmp(methodNames{iMethod},'rerefAllInclMast')
      
                        PreprocessedEEG_Musicians_rerefAllInclMast.(condNames{iCond}).(subName) = lwdata;
    
                        % reorder the subjects
                        [~, neworder]   = sort(lower(fieldnames(PreprocessedEEG_Musicians_rerefAllInclMast.(condNames{iCond}))));
                        PreprocessedEEG_Musicians_rerefAllInclMast.(condNames{iCond}) = orderfields(PreprocessedEEG_Musicians_rerefAllInclMast.(condNames{iCond}), neworder);

                    elseif strcmp(methodNames{iMethod},'rerefAllExclMast')
      
                        PreprocessedEEG_Musicians_rerefAllExclMast.(condNames{iCond}).(subName) = lwdata;
    
                        % reorder the subjects
                        [~, neworder]   = sort(lower(fieldnames(PreprocessedEEG_Musicians_rerefAllExclMast.(condNames{iCond}))));
                        PreprocessedEEG_Musicians_rerefAllExclMast.(condNames{iCond}) = orderfields(PreprocessedEEG_Musicians_rerefAllExclMast.(condNames{iCond}), neworder);
                    end

                % NonMusicians    
                elseif strcmp(group,Cfg.groups{2})

                    if strcmp(methodNames{iMethod},'rerefMast')

                        PreprocessedEEG_NonMusicians_rerefMast.(condNames{iCond}).(subName) = lwdata;
    
                        % reorder the subjects
                        [~, neworder]   = sort(lower(fieldnames(PreprocessedEEG_NonMusicians_rerefMast.(condNames{iCond}))));
                        PreprocessedEEG_NonMusicians_rerefMast.(condNames{iCond}) = orderfields(PreprocessedEEG_NonMusicians_rerefMast.(condNames{iCond}), neworder);
                    
                    elseif strcmp(methodNames{iMethod},'rerefAllInclMast')
      
                        PreprocessedEEG_NonMusicians_rerefAllInclMast.(condNames{iCond}).(subName) = lwdata;
    
                        % reorder the subjects
                        [~, neworder]   = sort(lower(fieldnames(PreprocessedEEG_NonMusicians_rerefAllInclMast.(condNames{iCond}))));
                        PreprocessedEEG_NonMusicians_rerefAllInclMast.(condNames{iCond}) = orderfields(PreprocessedEEG_NonMusicians_rerefAllInclMast.(condNames{iCond}), neworder);

                    elseif strcmp(methodNames{iMethod},'rerefAllExclMast')
      
                        PreprocessedEEG_NonMusicians_rerefAllExclMast.(condNames{iCond}).(subName) = lwdata;
    
                        % reorder the subjects
                        [~, neworder]   = sort(lower(fieldnames(PreprocessedEEG_NonMusicians_rerefAllExclMast.(condNames{iCond}))));
                        PreprocessedEEG_NonMusicians_rerefAllExclMast.(condNames{iCond}) = orderfields(PreprocessedEEG_NonMusicians_rerefAllExclMast.(condNames{iCond}), neworder);
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%    
        % Save all preprocessing % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        cd(Paths.preprocessing)
    
        if strcmp(group,Cfg.groups{1}) 
            
            

            if strcmp(methodNames{iMethod},'rerefMast')
                save('PreprocessedEEG_Musicians_rerefMast.mat','PreprocessedEEG_Musicians_rerefMast')
            elseif strcmp(methodNames{iMethod},'rerefAllInclMast')
                save('PreprocessedEEG_Musicians_rerefAllInclMast.mat','PreprocessedEEG_Musicians_rerefAllInclMast')
            elseif strcmp(methodNames{iMethod},'rerefAllExclMast')
                save('PreprocessedEEG_Musicians_rerefAllExclMast.mat','PreprocessedEEG_Musicians_rerefAllExclMast')             
            end

        elseif strcmp(group,Cfg.groups{2})
            if strcmp(methodNames{iMethod},'rerefMast')
                save('PreprocessedEEG_NonMusicians_rerefMast.mat','PreprocessedEEG_NonMusicians_rerefMast')
            elseif strcmp(methodNames{iMethod},'rerefAllInclMast')
                save('PreprocessedEEG_NonMusicians_rerefAllInclMast.mat','PreprocessedEEG_NonMusicians_rerefAllInclMast')
            elseif strcmp(methodNames{iMethod},'rerefAllExclMast')
                save('PreprocessedEEG_NonMusicians_rerefAllExclMast.mat','PreprocessedEEG_NonMusicians_rerefAllExclMast')             
            end
        end
    
        fprintf ('  ==> Preprocessed EEG structures saved \n')
    end
end


movefile(Paths.lw6OutsideMatlabPath,Paths.lw6)
