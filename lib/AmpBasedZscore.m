function EEG_zscores = AmpBasedZscore(sFFT,freqVec,subjects,iSubjects,iCond,iChan,EEG_zscores,analysisType,varargin)



%% Initial Setup %%
    
    global cfg Cfg EEG_fft CochMod Paths

    if nargin == 9
        
        Options     = varargin{1};
        optionNames = fieldnames(Options);
        
        if ~isempty(find(strcmp(optionNames,'AmpBasedZscoreBoundaries')))
            applyAmpBasedZscoreBoundaries   = true;
            lowerLimit                      = Options.AmpBasedZscoreBoundaries(1);
            upperLimit                      = Options.AmpBasedZscoreBoundaries(2);
            rangeMeanStdForFindpeaks        = Options.AmpBasedZscoreBoundaries;
        else
            applyAmpBasedZscoreBoundaries = false;
        end
        
        if ~isempty(find(strcmp(optionNames,'rereference')))
            rereference   = Options.rereference;
        end        
        if ~isempty(find(strcmp(optionNames,'group')))
            group   = Options.group;
        end 
        
    end

    condNames   = Cfg.condNames;
    frex        = 1/(Cfg.Stim.gridIOI*12) * [1:12*20];
    zStruct     = EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan};
    
    

%% Preperation of the AmpBasedZscore at the individual level %%

    %%%%%%%%%%%%%
    % Findpeaks %
    %%%%%%%%%%%%%

    % Set threshold for peak detection
    if applyAmpBasedZscoreBoundaries
        minRange    = dsearchn(freqVec',lowerLimit);
        maxRange    = dsearchn(freqVec',upperLimit);
        meanFFT     = mean(sFFT(minRange:maxRange));
        stdFFT      = std (sFFT(minRange:maxRange)); 
    else
        meanFFT     = mean(sFFT);
        stdFFT      = std (sFFT); 
    end

    % le fait de trier en descendant est pratique pour la méthode 12 pour ne garder que les 4 intensités les plus fortes. 
    [pks, locs] = findpeaks (sFFT,freqVec,'SortStr','descend','MinPeakProminence',meanFFT+2*stdFFT,'MinPeakHeight',0+2*stdFFT);


    % set boundaries in detected locs 
    if applyAmpBasedZscoreBoundaries
        clear Options
        Options.pks         = pks;
        Options.lowerLimit  = lowerLimit;
        Options.upperLimit  = upperLimit;
        [locs,Details,pks]  = locsBoundaries(locs,Options) ; 
    end
   
    % met(un)rel classifier
    [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  All prominent peaks (with and without 5Hz hamronics) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    AllProminentPeakmethodName              = [{'Amp_EEG_WithinGroup_PerCond_IndividualPart - All Peaks'};...
                                               {'Amp_EEG_WithinGroup_PerCond_IndividualPart - All Peaks (without 5Hz harmonics)'}];
    
    
%     AllProminentPeakmethodName              = [{'AmpIndPartPerCond - All Peaks'};...
%                                                {'AmpIndPartPerCond - All Peaks (without 5Hz harmonics)'}];
                                        
    for nMethods = 1:2
        
        % retrieve the index of this zscore method
        methodName                           = AllProminentPeakmethodName{nMethods};
        zMethodIdx                           = find(cell2mat(cellfun(@(x) strcmp(methodName, x), {zStruct.methodName}, 'UniformOutput', 0)));
        if isempty(zMethodIdx);              zMethodIdx = size(zStruct,2) + nMethods; end

        % remove 5Hz harmonics in the second iteration of nMethods
        if nMethods == 2
            clear Options
            Options.pks         = pks;
            Options.remove5Hz   = true;
            Options.lowerLimit  = lowerLimit;  % remove f1 if not done already
            [locs,Details,pks]  = locsBoundaries(locs,Options);
           
            
            % classify metRel and metUnrel frequencies
            [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);
            
        end
        
        if exist('metRelIdx','var') && exist('metUnrelIdx','var')  
            % computation of the zscore
            amps                                 = pks;
            zscores                              = (amps - nanmean(amps)) ./ nanstd(amps);
            zscoresMetRel                        = mean(zscores(metRelIdx));
            zscoresMetUnrel                      = mean(zscores(metUnrelIdx)); 

            % save the zscore
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).methodName     = {methodName};
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).amps           = amps;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).freq           = locs;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).frexIdx        = frexIdx; % dsearchn(freqVec',locs');
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).metRelIdx      = metRelIdx;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).metUnrelIdx    = metUnrelIdx;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).zMetRel        = zscoresMetRel;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).zMetUnrel      = zscoresMetUnrel;
        end
        
        clear amps zscores zscoresMetRel zscoresMetUnrel metRelIdxWithout5Hz idx2remove metRelIdx metUnrelIdx
        
    end
    
    

%% Preperation of the individual AmpBasedZscore based on locations from the averaged participants %%


    %%%%%%%%%%%%%%%%%%%%%%%%
    %  All prominent peaks %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    AllProminentPeakmethodName = [{'Amp_EEG_WithinGroup_AcrossCond_IndividualPart - All Peaks - Peaks sign in at least one condition'};...
                                  {'Amp_EEG_WithinGroup_AcrossCond_IndividualPart - All Peaks (without 5Hz harmonics) - Peaks sign in at least one condition'};...
                                  {'Amp_EEG_WithinGroup_AcrossCond_IndividualPart - All Peaks (without 5Hz harmonics) - Peaks sign in every conditions'}];
    

    for nMethods = 1:length(AllProminentPeakmethodName)
        
        % find the methodIdx to recover the significant peaks in the spectrum of the averaged participants and across conditions
        ampMethodUsed           = [{'Amp_EEG_WithinGroup_AcrossCond_AllPart - All Peaks - Peaks sign in at least one condition'};...
                                   {'Amp_EEG_WithinGroup_AcrossCond_AllPart - All Peaks (without 0.4 5Hz harmonics) - Peaks sign in at least one condition'};...
                                   {'Amp_EEG_WithinGroup_AcrossCond_AllPart - All Peaks (without 0.4 5Hz harmonics) - Peaks sign in every conditions'}];

        zStructAllPart          = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iChan};
        methodIdxAllPart        = find(cell2mat(cellfun(@(x) strcmp(ampMethodUsed(nMethods), x), {zStructAllPart.methodName}, 'UniformOutput', 0)));
        
        % retrieve the peaks from the averaged participants spectrum and across all conditions
        frexIdx                 = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iChan}(methodIdxAllPart).frexIdx;
        freq                    = freqVec(frexIdx);
        metRelIdx               = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iChan}(methodIdxAllPart).metRelIdx;
        metUnrelIdx             = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iChan}(methodIdxAllPart).metUnrelIdx;

        % retrieve the index of the current zscore method
        methodName                  = AllProminentPeakmethodName{nMethods};
        zStruct                     = EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan};
        zMethodIdx                  = find(cell2mat(cellfun(@(x) strcmp({methodName}, x), {zStruct.methodName}, 'UniformOutput', 0)));
        if isempty(zMethodIdx);     zMethodIdx = size(zStruct,2) + 1 ; end 
        
        % compute the zscore
        amps            = sFFT(frexIdx);
        zscores         = (amps - nanmean(amps)) ./ nanstd(amps);
        zscoresMetRel   = mean(zscores(metRelIdx));
        zscoresMetUnrel = mean(zscores(metUnrelIdx));
                           
        % store zscores
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).methodName   = {methodName};
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).amps         = amps;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).freq         = freq;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).frexIdx      = frexIdx;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).metRelIdx    = metRelIdx;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).metUnrelIdx  = metUnrelIdx;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).zMetRel      = zscoresMetRel;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).zMetUnrel    = zscoresMetUnrel;

        clear frexIdx metRelIdx metUnrelIdx
    end 
    
    
    
%% Ampbased zscore avec les pics du stim et de l'EEG within group %%

    if exist(fullfile(Paths.project,'Stimuli',['CochMod_',group,'_',rereference,'_',analysisType,'.mat'])) % if it exists, it should have been loaded before
        
        mixzMethodNames = [{'MixEEGStim_WithinGroup_AcrossCond_IndividualPart - All Peaks (without 5Hz harmonics) - Peaks sign in at least one condition'};...
                           {'MixEEGStim_WithinGroup_AcrossCond_IndividualPart - All Peaks (without 5Hz harmonics) - Peaks sign in every conditions'}];

        for iMethods = 1:2
        
            % prepare CochMod structure
            if or(iCond  == 1,iCond == 3); iStim = 2;
            elseif iCond == 2; iStim = 1;
            elseif iCond == 4; iStim = 3;
            end
            
            % find the zStimIdx for the Stim
            stimNames    = fieldnames(CochMod);
            zStimName    = {'Amp_Stim - Peaks sign in at least one condition';...
                            'Amp_Stim - Peaks sign in every conditions'};
            zStimStruct  = CochMod.(stimNames{iStim}).AN.zScoreOut;
            zStimIdx     = find(cell2mat(cellfun(@(x) strcmp(zStimName{iMethods}, x), {zStimStruct.methodName}, 'UniformOutput', 0)));

            % find the zEEGIdx for the EEG
            zEEGName    = {'Amp_EEG_WithinGroup_AcrossCond_IndividualPart - All Peaks (without 5Hz harmonics) - Peaks sign in at least one condition';...
                           'Amp_EEG_WithinGroup_AcrossCond_IndividualPart - All Peaks (without 5Hz harmonics) - Peaks sign in every conditions'};

            zEEGStruct  = EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan};
            zEEGIdx     = find(cell2mat(cellfun(@(x) strcmp(zEEGName{iMethods}, x), {zEEGStruct.methodName}, 'UniformOutput', 0)));
            
            % find the zMethod for this zscoring
            zMethodIdx  = find(cell2mat(cellfun(@(x) strcmp(mixzMethodNames{iMethods}, x), {zEEGStruct.methodName}, 'UniformOutput', 0)));
            if isempty(zMethodIdx);     zMethodIdx = size(zEEGStruct,2) + 1; end
            
            % concatenate all frequencies of interest
            stimFreq = CochMod.(stimNames{iStim}).AN.zScoreOut(zStimIdx).freq;
            if size(stimFreq,1) ~=1; stimFreq = stimFreq'; end
            freq        = [EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zEEGIdx).freq, ...
                           stimFreq];
             
            frexIdx     = unique(dsearchn(freqVec',freq'));
            locs        = freqVec(frexIdx);
            [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);

            
            % compute the zscore
            amps            = sFFT(frexIdx);
            zscores         = (amps - nanmean(amps)) ./ nanstd(amps);
            zscoresMetRel   = mean(zscores(metRelIdx));
            zscoresMetUnrel = mean(zscores(metUnrelIdx));

            % store zscores
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).methodName   = {mixzMethodNames{iMethods}};
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).amps         = amps;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).freq         = locs;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).frexIdx      = frexIdx;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).metRelIdx    = metRelIdx;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).metUnrelIdx  = metUnrelIdx;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).zMetRel      = zscoresMetRel;
            EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).zMetUnrel    = zscoresMetUnrel;

            clear frexIdx metRelIdx metUnrelIdx
        end
    else
        warning('Merged peaks from stim and EEG not computed as the CochMod.mat couldn''t be found')
    end


end

