function [EEG_zscores] = peakDetectionAndOverallZscore(EEG_zscores, EEG_fft, group, analysisType, varargin)

% computes the AmpBasedZscore at the group level : 
%       - for each condition seperately
%       - across conditions 


%% initial Setup 

    fprintf('\n\n-----------------------------\n Zscoring on data from all participants \n')
    
    % set global and main variables 
    global Cfg Paths Log 
    
    condNames       = Cfg.condNames;
    subjects        = fieldnames(Log.(group));
    gridIOI         = Cfg.Stim.gridIOI;
    frex            = 1/(gridIOI*12) * [1:12*20];
    nElec           = length(table2cell(Log.(group).AllParticipants.Preprocessing.ElecLabels));
    
    
    % optional arguments
    if nargin > 3
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
        
        if ~isempty(find(strcmp(optionNames,'iElec')))
            selectedElec = Options.iElec;
        else
            selectedElec = [1:68];
        end

        if ~isempty(find(strcmp(optionNames,'PeakDetection')))
            PeakDetectionMethod     = Options.PeakDetection.Method;
            PeakDetectionNSideBins  = Options.PeakDetection.nSideBins;
        else
            PeakDetectionMethod = 'findpeaks';
        end        
        
    else
        applyAmpBasedZscoreBoundaries = false;
    end    
    
    
    % define freqVec (+ check that all participants have an identical one)
%     for iSub = 1:length(subjects)
%         freqVec(iSub,:) = Log.(group).(subjects{iSub}).SpectralAnalysis.freqVec;
%     end
    
%     if isempty(find(diff(freqVec,1))) 
%         clear freqVec
        freqVec = Log.(group).(subjects{1}).SpectralAnalysis.freqVec;
        freqRes = freqVec(2) - freqVec(1);
%     end
warning('freqVec dimensions not double checked')
    
clear Options

%% AmpBasedZscore at the group level for each condition seperately %%

    for iCond = 1:length(condNames)
        for iElec = selectedElec
            
            %%%%%%%%%%%%%%%%%%%%%
            % Findpeaks or zSNR %
            %%%%%%%%%%%%%%%%%%%%%
            
            % extract fft + set threshold for peak detection
            sFFT        = EEG_fft.(condNames{iCond}).AllParticipants.(analysisType).data(iElec,:);
            
            if applyAmpBasedZscoreBoundaries % determine the mean and std according to a certain range
                minRange    = dsearchn(freqVec',lowerLimit);
                maxRange    = dsearchn(freqVec',upperLimit);
                meanFFT     = mean(sFFT(minRange:maxRange));
                stdFFT      = std (sFFT(minRange:maxRange)); 
            else
                minRange    = 1;
                maxRange    = length(freqVec);
                meanFFT     = mean(sFFT);
                stdFFT      = std (sFFT); 
            end
            

            clear pks locs
            
            if strcmp(PeakDetectionMethod,'zSNR') % advantage of looking at local maxima
                
                iLocs = 1;
                iPks  = 1;
                percSign  = Cfg.Analysis.zSNRSign; % 0.001 = value from Dzhelyova et al. 2020 but value with similar results to the findpeak function
                sign      = norminv(1 - percSign);
                
                % get adapted frequency vector indices (with restricted range and no possible negative frequencies)
                adaptedIdx      = minRange:maxRange;
                
                if minRange < PeakDetectionNSideBins(2)
                    adaptedIdx = adaptedIdx(1   + PeakDetectionNSideBins(2) : end - PeakDetectionNSideBins(2));
                end
              
                % select segment of interest, extract amplitudes, and compute zscore
                for iFreq = 1:length(adaptedIdx)
                    
                    % select segment 
                    AmpCentralBin   = sFFT(adaptedIdx(iFreq));
                    segmentSideBins = [adaptedIdx(iFreq) -  PeakDetectionNSideBins(2) : adaptedIdx(iFreq) -  PeakDetectionNSideBins(1),...
                                       adaptedIdx(iFreq) +  PeakDetectionNSideBins(1) : adaptedIdx(iFreq) +  PeakDetectionNSideBins(2)];              
                    AmpsSideBins    =  sFFT(segmentSideBins);                               
   
                    peakEmergence   =  (AmpCentralBin - mean(AmpsSideBins))/std(AmpsSideBins); % Rossion et al., 2020 et détails dans Dzhelyova et al. 2020. 
                    %peakEmergenceBis=  zscore([AmpsSideBins,AmpCentralBin]);
                    
                    % isequal(peakEmergence,peakEmergenceBis(end))
                                   
                    if peakEmergence > sign % Cfg.Analysis.Stats.sign
                        locs(iLocs) = freqVec(adaptedIdx(iFreq)); 
                        pks(iPks)   = sFFT(adaptedIdx(iFreq));
                        
                        iLocs       = iLocs+1;
                        iPks        = iPks+1;
                                        
%                         freqqq = freqVec(segmentSideBins(1):segmentSideBins(end));
%                         amps   = sFFT(segmentSideBins(1):segmentSideBins(end));
% 
%                         plot([freqqq],zscore(amps))
% 
%                         line([freqqq(1) freqqq(end)],[Cfg.Analysis.Stats.sign Cfg.Analysis.Stats.sign])
%                         title(num2str(freqVec(adaptedIdx(iFreq))))
%                         %pause(0.5)
%                         close all
                        
                    end
                end
                
                if iLocs == 1
                    error('No peak found')
                end

                
            if iElec == 68 
                subplot(4,1,iCond)
                [pks_findpeaks,locs_findpeaks]      = findpeaks (sFFT,freqVec,'MinPeakProminence',meanFFT+2*stdFFT,...
                                                                 'MinPeakHeight',0+2*stdFFT);

                stem(freqVec,sFFT,'k','marker','none'); hold on
                scatter(locs,pks,'g')
                scatter(locs_findpeaks,pks_findpeaks+0.01,'r')
                xlim([0 30])
                ylim([0 0.3])
                if iCond == 4
                    pause(1)
                end
            end
                
            else    

                [pks,locs]      = findpeaks (sFFT,freqVec,'MinPeakProminence',meanFFT+2*stdFFT,...
                                                          'MinPeakHeight',0+2*stdFFT);
            end
            

                                                  
            % set boundaries in detected locs 
            if applyAmpBasedZscoreBoundaries
                clear Options
                Options.pks         = pks;
                Options.lowerLimit  = lowerLimit;
                Options.upperLimit  = upperLimit;
                [locs,Details,pks]  = locsBoundaries(locs,Options) ; 
            end

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  All prominent peaks (with and without 5Hz hamronics) %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            AllProminentPeakmethodName = [{'Amp_EEG_WithinGroup_PerCond_AllPart - All Peaks'};...
                                          {'Amp_EEG_WithinGroup_PerCond_AllPart - All Peaks (without 0.4 5Hz harmonics)'}];
                      
%             AllProminentPeakmethodName = [{'AmpAllPartPerCond - All Peaks'};...
%                                           {'AmpAllPartPerCond - All Peaks (without 0.4 5Hz harmonics)'}];
%           
            for nMethods = 1:2 
                % set boundaries for the findpeak detection
                if nMethods == 2
                    clear Options
                    Options.pks         = pks;
                    Options.remove5Hz   = true;
                    Options.lowerLimit  = lowerLimit;  % remove f1 if not done already
                    [locs,Details,pks]  = locsBoundaries(locs,Options);
                end

                % classify metRel and metUnrel frequencies
                [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);
             
                % zscore
                if exist('metRelIdx','var')     && ~isempty(metRelIdx) && ...
                   exist('metUnrelIdx','var')   && ~isempty(metUnrelIdx)
                    
                    % find the correct methodName in the EEG_zscore structure
                    methodName  = AllProminentPeakmethodName {nMethods};
                    zStruct     = EEG_zscores.(condNames{1}).AllParticipants.(analysisType){1,iElec};
                    methodIdx   = find(cell2mat(cellfun(@(x) strcmp({methodName}, x), {zStruct.methodName}, 'UniformOutput', 0)));
                    
                    if isempty(methodIdx);  methodIdx = size(zStruct,2)+1; end
                    
                    % compute zscore
                    amps            = pks;
                    zscores         = (amps - nanmean(amps)) ./ nanstd(amps);
                    zscoresMetRel   = mean(zscores(metRelIdx));
                    zscoresMetUnrel = mean(zscores(metUnrelIdx));
                    
                    % store zscores
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).methodName   = {methodName};
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).amps         = amps;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).freq         = locs;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).frexIdx      = frexIdx;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).metRelIdx    = metRelIdx;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).metUnrelIdx  = metUnrelIdx;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).zMetRel      = zscoresMetRel;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).zMetUnrel    = zscoresMetUnrel;
                end
                
                clear metRelIdx metUnrelIdx amps idx2remove
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % plot findpeak results %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            
%             if iElec == 68
%                 
%                 % retrieve data
%                 metRelIdx   = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).metRelIdx;
%                 metUnrelIdx = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).metUnrelIdx;
% 
%                 % plot
%                 disp(iCond)
%                 subplot(4,1,iCond)
%                 plot(freqVec,sFFT,'k');hold on; 
%                 stem(locs(metRelIdx), pks(metRelIdx), 'Color','r','LineWidth',1.25,'marker','none')
%                 stem(locs(metUnrelIdx), pks(metUnrelIdx), 'Color','b','LineWidth',1.25,'marker','none')
%                 
%                 % layout
%                 set(gca,'lineWidth',1.25,'TickDir','out','xlim',[0 11],'XTick',[0:0.5:11],'ylim',[-0.5 1])
%                 box off
%                 title(condNames{iCond})
%                 
%                 % save
%                 if iCond == 4
%                     set(gcf,'PaperPositionMode','auto')
%                     print('PeakDetectionAllParticipants.jpg','-djpeg','-r300') 
%                     savefig('PeakDetectionAllParticipants.fig')
%                 end
%                 
%                 clear metRelIdx metUnrelIdx 
%             end
        end
    end
            
%% AmpBasedZscore at the group level across condition %%

    % prepare the comparison of the peak detection in all conditions (row
    % with freqVec and boleans with true if peak detected at the frequency)
    acrossCondLocs        = zeros(length(freqVec),length(condNames)+2); % first extra for the freqVec and second extra for the sum
    acrossCondLocs(:,1)   = freqVec;

    AllProminentPeakmethodName = [{'Amp_EEG_WithinGroup_AcrossCond_AllPart - All Peaks - Peaks sign in at least one condition'};...
                                  {'Amp_EEG_WithinGroup_AcrossCond_AllPart - All Peaks (without 0.4 5Hz harmonics) - Peaks sign in at least one condition'};...
                                  {'Amp_EEG_WithinGroup_AcrossCond_AllPart - All Peaks (without 0.4 5Hz harmonics) - Peaks sign in every conditions'}];
    
%     AllProminentPeakmethodName = [{'AmpAllPartAcrossCond - All Peaks - Peaks sign in at least one condition'};...
%                                   {'AmpAllPartAcrossCond - All Peaks (without 0.4 5Hz harmonics) - Peaks sign in at least one condition'};...
%                                   {'AmpAllPartAcrossCond - All Peaks (without 0.4 5Hz harmonics) - Peaks sign in every conditions'}];
                              
    for iElec = selectedElec      
        locsStart   = 1; 
        for iCond = 1:length(condNames)
            
            clear templocs
            
            % retrieve previous freq detected in all particiants for each condition seperately
            ampMethodUsed       = 'Amp_EEG_WithinGroup_PerCond_AllPart - All Peaks';
            zStruct             = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}; 
            previousMethodIdx   = find(cell2mat(cellfun(@(x) strcmp({ampMethodUsed}, x), {zStruct.methodName}, 'UniformOutput', 0)));
            
            if ~isempty (previousMethodIdx)
                templocs    = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(previousMethodIdx).freq;
                nlocs       = length(templocs);
                locsEnd     = locsStart + nlocs-1;

                locs(locsStart : locsEnd) = templocs;

                locsStart   = locsStart + nlocs;
                
                % peaks significant in every conditions
                for iLocs = 1:nlocs
                    [~,idx] = min(abs(freqVec-templocs(iLocs)));
                    acrossCondLocs(idx,iCond+1) = 1;
                end
            end
        end
        
        % if no locs found
        if isempty(locs);  error('missing locs in this cond and elec');  end
        
        % else only keep unique locs
        locs = unique(sort(locs));
            
        % prepare method with and without the 1st, 12th and harmonics
        for nMethods = 1:3
            
            % set boundaries for the findpeak detection
            if nMethods == 2
                clear Options
                Options.remove5Hz   = true;
                Options.lowerLimit  = lowerLimit;  % remove f1 if not done already
                [locs,Details]      = locsBoundaries(locs,Options);
            end
            
            % Third condition : only keep the peaks significantly present in all conditions
            if nMethods == 3
                clear idx
                acrossCondLocs(:,iCond+2) = sum(acrossCondLocs(:,[2:iCond+1]),2);
                idx     = find(acrossCondLocs(:,end)==4);
                locsIdx = dsearchn(locs',freqVec(idx)');
                locs    = locs(locsIdx);
            end
            
            % classify metRel and metUnrel frequencies
            [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);
                    
            % zscore
            if exist('metRelIdx','var') && exist('metUnrelIdx','var') 
                for iCond = 1:length(condNames)
                
                    % find the correct methodName in the EEG_zscore structure
                    methodName  = AllProminentPeakmethodName {nMethods};
                    zStruct     = EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){1,iElec};
                    methodIdx   = find(cell2mat(cellfun(@(x) strcmp({methodName}, x), {zStruct.methodName}, 'UniformOutput', 0)));

                    if isempty(methodIdx);  methodIdx = size(zStruct,2)+1; end

                    % retrieve EEG amplitude at the new locations
                    frexIdx         = dsearchn(freqVec',locs');
                    amps            = EEG_fft.(condNames{iCond}).AllParticipants.(analysisType).data(iElec,frexIdx);
                    
                    % compute zscore
                    zscores         = (amps - nanmean(amps)) ./ nanstd(amps);
                    zscoresMetRel   = mean(zscores(metRelIdx));
                    zscoresMetUnrel = mean(zscores(metUnrelIdx));
                
                    % store zscores
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).methodName   = {methodName};
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).amps         = amps;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).freq         = locs;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).frexIdx      = frexIdx;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).metRelIdx    = metRelIdx;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).metUnrelIdx  = metUnrelIdx;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).zMetRel      = zscoresMetRel;
                    EEG_zscores.(condNames{iCond}).AllParticipants.(analysisType){:,iElec}(methodIdx).zMetUnrel    = zscoresMetUnrel;
                    
                    Log.AllParticipants.SpectralAnalysis.freqVec = freqVec;
                end
            end
            clear metRelIdx metUnrelIdx amps idx2remove
        end
    end
end

