function [subjects,subjectsIdx,output] = EEGSpectralAnalysis(analysisType, group, rereference)

% Spectral analysis of the EEG

%%%%%%%%%%%%%%%%%%
% Initial setups %
%%%%%%%%%%%%%%%%%%

    global Cfg Paths Log  

    condNames   = Cfg.condNames;
    groups      = Cfg.groups;

    if isfolder(Paths.lw6)
        movefile(Paths.lw6,Paths.lw6OutsideMatlabPath)
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spectral Analysis based on the Averaged Data %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(analysisType,'AveragedTrials') % AveragedTrialsSpectralAnalysis == 1 
       
        cd(Paths.preprocessing)                       
        matDir  = dir(['PreprocessedEEG_',group,'.mat']);         
        data    = load(fullfile(matDir.folder,matDir.name));
        field   = fieldnames(data);
        data    = data.(field{1});

        subjects    = fieldnames(data.(condNames{1}));
        subjectsIdx = 1:length(subjects);
                       
        cd(fullfile(Paths.LW,(group),'Analysis/EEG/Spectral Analysis/Averaged'))
        
        %%%%%%%    
        % FFT %
        %%%%%%%
    
        i=1;
    
        for iCond = 1:length(condNames)
            for iSubjects = subjectsIdx
                
                % Extraction of the data
                lwdata = data.(condNames{iCond}).(subjects{iSubjects});
    
                % Channel pools
                chanLabels   = {table2cell(Log.(group).(subjects{iSubjects}).Preprocessing.ElecLabels)};
                frontalPool  = {{'F1','Fz','F2','FC1','FCz','FC2','C1','Cz','C2'}};
                
                option       = struct('name','avgAll','channels', chanLabels,'suffix','new_chan','is_save',0);
                lwdata       = FLW_new_channel_averaged.get_lwdata(lwdata,option);
                
                option       = struct('name','avgFront','channels', frontalPool,'suffix','new_chan','is_save',0);
                lwdata       = FLW_new_channel_averaged.get_lwdata(lwdata,option);
                
                % update the elecLabel in log file
                Log.(group).(subjects{iSubjects}).Preprocessing.ElecLabels = cell2table({lwdata.header.chanlocs.labels}');
                 
                % FFT
                option  = struct('output','amplitude',...
                                 'half_spectrum',1,...
                                 'suffix','fft','is_save',Cfg.Analysis.FFT.lwSave);
                lwdata  = FLW_FFT.get_lwdata(lwdata,option);
                
                
                % Baseline correction
                option = struct('xstart',Cfg.Analysis.Bl_snr.lowBin,...
                                'xend',Cfg.Analysis.Bl_snr.highBin,...
                                'num_extreme',0,...
                                'operation','subtract',...
                                'suffix','bl_snr','is_save',Cfg.Analysis.Bl_snr.lwSave);
                lwdata = FLW_baseline_SNR.get_lwdata(lwdata,option);
                
                % save the fft (and squeeze empty dimensions)
                output.(condNames{iCond}).(subjects{iSubjects}).(analysisType).header    = lwdata.header;
                output.(condNames{iCond}).(subjects{iSubjects}).(analysisType).data      = squeeze(lwdata.data);       
                
                disp(['Condition = ',condNames{iCond}, ';  Subject = ', subjects{iSubjects}])
                
                % store the xstep to add it later in the AllParticipants
                xstep(i) = lwdata.header.xstep;
                
                % freqVec
                Log.(group).(subjects{iSubjects}).SpectralAnalysis.freqVec = 0 : xstep(i) : lwdata.header.datasize(6)*xstep(i)-xstep(i);
                
                i=i+1;
    
                clear lwdata
            end
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Spectral Analysis based on the Non-Averaged Data %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif strcmp(analysisType,'NotAveragedTrials') % Not Averaged Data

        % get subject directory
        if      strcmp(group, groups{1}); iGroups = 1;
        elseif  strcmp(group, groups{2}); iGroups = 2;
        end

        cd(Paths.dataDir{iGroups})
        subDir      = dir(fullfile(Paths.dataDir{iGroups},'ep_edit*.lw6'));
        subjects    = extractBetween({subDir(1:end).name},'ep_edit ','.lw6');
    
        [subjectsIdx,~] = listdlg('PromptString',{'Select the participants you want to preprocess'},...
                                    'ListString',subjects);
        i=1;

        for iSubjects = subjectsIdx
            for iCond = 1:length(condNames)

                subName     = subjects{iSubjects};
    
                %%%%%%%%%%%%%%%%%%%%%%%
                % Load sp_filter data %
                %%%%%%%%%%%%%%%%%%%%%%%
        
                cd(fullfile(Paths.LW,(group),'Preprocessing/EEG/ICA/MultiRereference'))
                subDir = dir(['sp_filter*',rereference,'*',(condNames{iCond}),'*',(subName),'.lw6']);  
        
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

                %%%%%%%%%%%%%%%%%
                % Channel pools % 
                %%%%%%%%%%%%%%%%%
                cd(fullfile(Paths.LW,group,'Analysis/EEG/Spectral Analysis/NotAveraged'))

                chanLabels   = {table2cell(Log.(group).(subjects{iSubjects}).Preprocessing.ElecLabels)};
                frontalPool  = {{'F1','Fz','F2','FC1','FCz','FC2','C1','Cz','C2'}};
                
                option       = struct('name','avgAll','channels', chanLabels,'suffix','new_chan','is_save',0);
                lwdata       = FLW_new_channel_averaged.get_lwdata(lwdata,option);
                
                option       = struct('name','avgFront','channels', frontalPool,'suffix','new_chan','is_save',0);
                lwdata       = FLW_new_channel_averaged.get_lwdata(lwdata,option);
                
                % update the elecLabel in log file
                Log.(group).(subjects{iSubjects}).Preprocessing.ElecLabels = cell2table({lwdata.header.chanlocs.labels}');
                
                %%%%%%%
                % FFT %
                %%%%%%%
                option  = struct('output','amplitude',...
                                 'half_spectrum',1,...
                                 'suffix','fft','is_save',Cfg.Analysis.FFT.lwSave);
                lwdata  = FLW_FFT.get_lwdata(lwdata,option);
                
                %%%%%%%%%%%%%%%%%%%%%%%
                % Baseline correction % 
                %%%%%%%%%%%%%%%%%%%%%%%

                option = struct('xstart',Cfg.Analysis.Bl_snr.lowBin,...
                                'xend',Cfg.Analysis.Bl_snr.highBin,...
                                'num_extreme',0,...
                                'operation','subtract',...
                                'suffix','bl_snr','is_save',1);
                lwdata = FLW_baseline_SNR.get_lwdata(lwdata,option); 

                %%%%%%%%%%%%%
                % Averaging %
                %%%%%%%%%%%%%

                option = struct('operation','average','suffix','avg','is_save',1);
                lwdata = FLW_average_epochs.get_lwdata(lwdata,option);


                % save the fft (and squeeze empty dimensions)
                output.(condNames{iCond}).(subjects{iSubjects}).(analysisType).header    = lwdata.header;
                output.(condNames{iCond}).(subjects{iSubjects}).(analysisType).data      = squeeze(lwdata.data);       
                
                disp(['Condition = ',condNames{iCond}, ';  Subject = ', subjects{iSubjects}])
                
                % store the xstep to add it later in the AllParticipants
                xstep(i) = lwdata.header.xstep;
                
                % freqVec
                Log.(group).(subjects{iSubjects}).SpectralAnalysis.freqVec = 0 : xstep(i) : lwdata.header.datasize(6)*xstep(i)-xstep(i);
                
                i=i+1;
    
                clear lwdata

            end
        end
    end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Mean of all participants %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computation 1: irrespective of the condition order group
% computation 2: depending on the condition order group (ABAB1 first or CDEF first)

    for iCond = 1:length(condNames)
        
        % Extraction of the data for each participant
        for iSubjects = subjectsIdx
            lwdata (:,:,iSubjects) = output.(condNames{iCond}).(subjects{iSubjects}).(analysisType).data;
        end

        % Mean of all participants (irrespective of their condition group)
            output.(condNames{iCond}).('AllParticipants').(analysisType).data              = mean(lwdata,3);
            output.(condNames{iCond}).('AllParticipants').(analysisType).header.datasize   = [1,size(lwdata,1),1,1,1,size(lwdata,2)];

            if length(unique(xstep)) == 1
                output.(condNames{iCond}).('AllParticipants').(analysisType).header.xstep      = unique(xstep);
            else 
                error('Different xstep in participants')
            end
            
        % Mean of all participants based on their condition group
            % preallocation
            nElec           = length(chanLabels{1});
            lwdataGroup1    = NaN(nElec,size(lwdata,2),length(subjectsIdx));
            lwdataGroup2    = NaN(nElec,size(lwdata,2),length(subjectsIdx));
            nSubGroup1      = 0;
            nSubGroup2      = 0;


            % Extraction of the data for each participant depending on their group
            for iSubjects = subjectsIdx
                if strcmp(Log.(group).(subjects{iSubjects}).condOrder{1},condNames{1}) % condition order group 1
                    lwdataGroup1 (:,:,iSubjects) = output.(condNames{iCond}).(subjects{iSubjects}).(analysisType).data;
                    nSubGroup1 = nSubGroup1+1;
                else
                    lwdataGroup2 (:,:,iSubjects) = output.(condNames{iCond}).(subjects{iSubjects}).(analysisType).data;
                    nSubGroup2 = nSubGroup2+1;
                end
            end

            % Mean of participants from group1
            output.(condNames{iCond}).('AllParticipantsGroup1').(analysisType).data              = nanmean(lwdataGroup1,3);
            output.(condNames{iCond}).('AllParticipantsGroup1').(analysisType).nSubjects         = nSubGroup1;

            % Mean of participants from group2
            output.(condNames{iCond}).('AllParticipantsGroup2').(analysisType).data              = nanmean(lwdataGroup2,3);
            output.(condNames{iCond}).('AllParticipantsGroup2').(analysisType).nSubjects         = nSubGroup2;
    end

    % Update the subjects variables for future use in the zscore function
    subjectsIdx = [subjectsIdx length(subjectsIdx)+1]; % Add the AllParticipants
    subjects    = fieldnames(output.(condNames{iCond}));
    
    Log.(group).AllParticipants.Preprocessing.ElecLabels            = Log.(group).(subjects{iSubjects}).Preprocessing.ElecLabels;
    Log.(group).('AllParticipantsGroup1').Preprocessing.ElecLabels  = Log.(group).(subjects{iSubjects}).Preprocessing.ElecLabels;
    Log.(group).('AllParticipantsGroup2').Preprocessing.ElecLabels  = Log.(group).(subjects{iSubjects}).Preprocessing.ElecLabels;


fprintf ('  ==> Spectral analysis computed \n')

cd(Paths.analysis)

movefile(Paths.lw6OutsideMatlabPath,Paths.lw6)


end

