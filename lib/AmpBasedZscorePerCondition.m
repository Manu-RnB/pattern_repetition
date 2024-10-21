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

    
    
%% Preperation of the individual AmpBasedZscore based on locations from the averaged participants per condition %%
    
newMethodName       = {'Amp_EEG_WithinGroup_PerCond_IndividualPart - All Peaks';...
                       'Amp_EEG_WithinGroup_PerCond_IndividualPart - All Peaks (without 0.4 5Hz harmonics)'};

previousMethodName  = {'Amp_EEG_WithinGroup_PerCond_AllPart - All Peaks';...
                       'Amp_EEG_WithinGroup_PerCond_AllPart - All Peaks (without 0.4 5Hz harmonics)'};

 for nMethods = 1:length(newMethodName)  
     
        zStructAllPart          = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iChan};
        methodIdxAllPart        = find(cell2mat(cellfun(@(x) strcmp(previousMethodName(nMethods), x), {zStructAllPart.methodName}, 'UniformOutput', 0)));
        
        % retrieve the peaks from the averaged participants spectrum and across all conditions
        frexIdx                 = zStructAllPart(methodIdxAllPart).frexIdx;
        freq                    = freqVec(frexIdx);
        metRelIdx               = zStructAllPart(methodIdxAllPart).metRelIdx;
        metUnrelIdx             = zStructAllPart(methodIdxAllPart).metUnrelIdx;
     
        % retrieve the index of the current zscore method
        zStruct                     = EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan};
        zMethodIdx                  = find(cell2mat(cellfun(@(x) strcmp({newMethodName{nMethods}}, x), {zStruct.methodName}, 'UniformOutput', 0)));
        if isempty(zMethodIdx);     zMethodIdx = size(zStruct,2) + 1 ; end 

        % compute the zscore
        amps            = sFFT(frexIdx);
        zscores         = (amps - nanmean(amps)) ./ nanstd(amps);
        zscoresMetRel   = mean(zscores(metRelIdx));
        zscoresMetUnrel = mean(zscores(metUnrelIdx));
                           
        % store zscores
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).methodName   = newMethodName{nMethods};
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).amps         = amps;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).freq         = freq;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).frexIdx      = frexIdx;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).metRelIdx    = metRelIdx;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).metUnrelIdx  = metUnrelIdx;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).zMetRel      = zscoresMetRel;
        EEG_zscores.(condNames{iCond}).(subjects{iSubjects}).(analysisType){:,iChan}(zMethodIdx).zMetUnrel    = zscoresMetUnrel;
        
        clear frexIdx metRelIdx metUnrelIdx
 end


end

