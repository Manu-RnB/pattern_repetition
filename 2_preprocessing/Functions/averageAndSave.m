function averageAndSave(subName, group)


    global Cfg Log Paths
    
    fprintf ('  -> Averaging and saving \n')
    
    condNames   = Cfg.condNames;
    
 
    % load previous data structure
    if strcmp(group,Cfg.groups{1})
        load(fullfile(Paths.preprocessing,'PreprocessedEEG_Musicians.mat'))
    elseif strcmp(group,Cfg.groups{2})
        load(fullfile(Paths.preprocessing,'PreprocessedEEG_NonMusicians.mat'))
    end
    


    for iCond = 1: size(condNames,2)
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % Load sp_filter data %
        %%%%%%%%%%%%%%%%%%%%%%%
        cd(fullfile(Paths.LW,group,'Preprocessing/EEG/ICA'))
        subDir = dir(['sp_filter*',(condNames{iCond}),'*',(subName),'*lw6']);
            
        if length(subDir)>1
            [ICAFileName, ICAFolder] = uigetfile(['sp_filter*',(condNames{iCond}),'*',(subName),'*lw6']);
            
        elseif isempty(subDir)
            warning('No sp_filter file found. Please look for it manually or select the rereference file if no ICA was computed')
            
            cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Rereference'))
            
            rerefDir = dir(['reref*',(condNames{iCond}),'*',(subName),'*lw6']);
            
            if length(rerefDir) == 1 && strcmp(questdlg(['Would you like to load this file: ',rerefDir.name]),'Yes')
                ICAFileName     = rerefDir.name;
                ICAFolder       = rerefDir.folder;
            else
                [ICAFileName, ICAFolder] = uigetfile(['reref*',(condNames{iCond}),'*',(subName),'*lw6']);
            end
        else     
            ICAFileName    = subDir.name;
            ICAFolder      = subDir.folder;
        end
       
        option      = struct('filename',fullfile(ICAFolder,ICAFileName));
        lwdata      = FLW_load.get_lwdata(option);  

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store the unaveraged preprocessing % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         if Cfg.Preprocessing.NotAverage.lwSave
%             EEG_preprocessedNotAveraged.(condNames{iCond}).(subName) = lwdata; 
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Average all trials %
        %%%%%%%%%%%%%%%%%%%%%%
        
        cd(fullfile(Paths.LW,group,'preprocessing/EEG/Average'))
        
        option = struct('operation','average','suffix','avg','is_save',1);
        lwdata = FLW_average_epochs.get_lwdata(lwdata,option);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store the averaged preprocessing % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Musicians
        if strcmp(group,Cfg.groups{1})
            
            PreprocessedEEG_Musicians.(condNames{iCond}).(subName) = lwdata;
            
            % reorder the subjects
            [~, neworder]   = sort(lower(fieldnames(PreprocessedEEG_Musicians.(condNames{iCond}))));
            PreprocessedEEG_Musicians.(condNames{iCond}) ...
                            = orderfields(PreprocessedEEG_Musicians.(condNames{iCond}), neworder);
             
        % NonMusicians    
        elseif strcmp(group,Cfg.groups{2})
            
            PreprocessedEEG_NonMusicians.(condNames{iCond}).(subName) = lwdata;
            
            % reorder the subjects
            [~, neworder]   = sort(lower(fieldnames(PreprocessedEEG_NonMusicians.(condNames{iCond}))));
            PreprocessedEEG_NonMusicians.(condNames{iCond}) ...
                            = orderfields(PreprocessedEEG_NonMusicians.(condNames{iCond}), neworder);
        end
    end

    
%%%%%%%%%%%%%%%    
% ReadMe File %
%%%%%%%%%%%%%%%

    % Add a comment to the preprocessing
    commentQuest = questdlg ('Would you like to add a specific comment about this analysis?', ...
                                        'Add a comment','Yes','No','No');
    switch commentQuest
        case 'No'
            commentTxt = {' '};
        case 'Yes'
            commentTxt = inputdlg ('Please enter your comment here','ReadMe comment');
    end
    

    % Write ReadMe file
    try
        cd(fullfile(Paths.preprocessing,'ReadMe',group))
    catch 
        mkdir(fullfile(Paths.preprocessing,'ReadMe',group))
        cd(fullfile(Paths.preprocessing,'ReadMe',group))
    end
    ReadMe = fopen(['ReadMe_',subName,'.txt'],'w');
    
    fprintf(ReadMe,['\nSummary of the Preprocessing \n---------------------------\n\n', ...
                    'Subject ID : ', subName, '\n\n', ...
                    'Group : ', group, '\n\n', ...
                    'Date of analysis : ', datestr(now),'\n\n', ...
                    'Comment : ', commentTxt{1},'\n\n\n', ...
                    'Butterworth filter  \n------------------ \n \t filter type : ',Cfg.Preprocessing.HPfilter.filterType  , '\n \t low cutoff : ',num2str(Cfg.Preprocessing.HPfilter.cutOff) , '\n\t filter order : ',num2str(Cfg.Preprocessing.HPfilter.Order),'\n\n', ...
                    'Removal of unused channels \n--------------------------\n\n', ...
                    'Rereference to : ', Cfg.Preprocessing.Reref.Elec, '\n\n', ...
                    'Channel interpolation   \n----------------------\n \t interpolated channel(s) : ', Log.(group).(subName).Preprocessing.Interpolation.InterpolatedElectrode, '\n \t number of closest electrodes selected : ',Log.(group).(subName).Preprocessing.Interpolation.ClosestElectrodes,'\n\n', ...
                    'Independent component analysis \n------------------------------ \n \t number of ICs : ', num2str(Log.(group).(subName).Preprocessing.ICA.nIC),...
                    ]);
    
    fclose(ReadMe);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Save all preprocessing % 
%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd(Paths.preprocessing)

    if strcmp(group,Cfg.groups{1})
        save('PreprocessedEEG_Musicians.mat','PreprocessedEEG_Musicians','-append')
    elseif strcmp(group,Cfg.groups{2})
        save('PreprocessedEEG_NonMusicians.mat','PreprocessedEEG_NonMusicians','-append')
    end
    
    fprintf ('  ==> Preprocessed EEG structures saved \n')
    
    cd(Paths.project)
    save('Log.mat','Log','-append')
    

end

