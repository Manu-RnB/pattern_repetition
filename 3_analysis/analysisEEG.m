%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% NonRepComplex %%%%%%%%%%%%%%
%%%%%%%%%%% Data Analysis - EEG %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script structure
% ----------------

% 1) FFT computation
% 2) Selection of prominent peak in the stimulus
% 3) Selection of prominent peak in the EEG (group-averaged fft per condition)
% 4) Merging of the frequencies of interest (FOI) of the stim and EEG
% 5) Mix Stim/EEG zscore with boundaries and a restriction of the # of FOI
%    a) [1.25 5]Hz
%    b) [1.25 30]Hz
% 6) Stim based zscore with [1.25 5]Hz boundaries and a restriction of the # of FOI
% 7) Zscore of the stimulus for each method
% 8) Figure (Comparison of all methods)
% 9) Figure (spectral domain visualization)
% 9) Prepare table for RStudio


%% Paths + Preprocessing parameters


% Initial setups 
% --------------
  
    global Cfg Log Paths CochMod
    
    [Cfg, Paths, Log] = getCfg(projectPath); 
    
    printFig = true;
    
% Analysis parameters
% -------------------
    stimNames   = {'AAAA','ABAB','CDEF'};
    condNames   = Cfg.condNames;
    groups      = Cfg.groups;
    groupNames  = {'Musicians', 'Non-musicians'};
    condOrderNames ...
                = {'maximum prior context', 'minimal prior context'};
    electrode   = 68; % Fronto-central pool
    methodNames = {'MixLargeRange','MixShortRange','StimBased', 'PreviousMix'};
    frex        = 1/(Cfg.Stim.gridIOI*12) * [1:12*20];
    
    
    
     
%% FFT

% Description
% -----------

% 1) Load data
% 2) Pools of electrodes (Avg, fronto-central)
% 3) FFT & Baseline correction
% 4) Grand averages

for  iGroup = 1:length(groups)
    
%     % parameters
%     group           = groups{iGroup};
%     analysisType    = 'AveragedTrials';
%     rereference     = 'rerefMast';
%     
%     [subjects,subjectsIdx,output] ...
%                     = EEGSpectralAnalysis(analysisType, group, rereference);
%     EEG.FFT.(group) = output;
%     
%     idx2keep        = find(contains(subjects,'sub'));
%     subjects        = subjects(idx2keep);
%      
%     EEG.FFT.(group).subjects = subjects; 
end

%% Selection of prominent peak in the stimulus

% Description
% -----------

% 1) Load the cochlear model
% 2) Mean of all auditory nerve fibers
% 3) FFT
% 4) Selection of prominent peaks
% 5) set boundaries (mainly to remove 5Hz and multiples)

figure('Color', [1 1 1], 'Position',[1 1 1200 900],'name','Selection of prominent peak in the stim')


for iStim = 1:length(stimNames)
     
    % load cochlear model
    cd(fullfile(projectPath,'1_stimulus'))
    CochModDir  = dir(['UREAR_BMF[64]*',stimNames{iStim},'*F3*']);
    data        = load([CochModDir.folder,filesep,CochModDir.name]);

    % Mean of all fibers 
    data.output.AN.meanSOut = mean(data.output.AN.an_sout,1);
        
    % FFT
    [freqVec, sFFT] = FFT(data.output.AN.meanSOut, data.output.AN.fs);
    
    CochMod.(stimNames{iStim}).F3.AN.freqVec  	= freqVec;
    CochMod.(stimNames{iStim}).F3.AN.sFFT       = sFFT;        
    
    % selection of prominent peaks (threshold = mean + 2*std)
    rangelimits = dsearchn(freqVec',[1.249 30]');
    range       = rangelimits(1):rangelimits(2);
    
    meanFFT     = mean(sFFT(range));
    stdFFT      = std(sFFT(range));
    
    [pks,locs]  = findpeaks (sFFT,freqVec,'MinPeakProminence',meanFFT+2*stdFFT);
                                      
    CochMod.(stimNames{iStim}).F3.AN.peaksDetection.locs = locs;
    CochMod.(stimNames{iStim}).F3.AN.peaksDetection.pks  = pks;
    
    % figure
    subplot(3,1,iStim)
    stem(freqVec,sFFT,'k','marker','none'); hold on
    scatter(locs,pks,'r')
    xlim([0 10])
    title((stimNames{iStim}))
    box off
   
    clear pks locs
end

%% Selection of prominent peak in the EEG

% Description
% -----------

% 1) Load the spectral analysis of the EEG
% 2) Selection of prominent peaks

iPlot = 1;
figure('Color', [1 1 1], 'Position',[1 1 1200 900],'name','Selection of prominent peak in the EEG')

if ~exist('EEG.FFT')
    load(fullfile(projectPath,'3_analysis/SpectralAnalysis_Musicians.mat'));
    load(fullfile(projectPath,'3_analysis/SpectralAnalysis_NonMusicians.mat'));
end

for iCond = 1:length(condNames)
    for iGroup = 1:length(groups)
    
        group       = groups{iGroup};
        condName    = condNames{iCond};
             
        % load magnitude spectrum
        if exist('EEG.FFT')
            sFFT        = EEG.FFT.(group).(condName).AllParticipants.AveragedTrials.data(electrode,:);
            subjects    = EEG.FFT.(group).subjects;
        else
            if strcmp(group,'Musicians')
                sFFT        = AnalyzedEEG_Musicians_blsnrFFT.(condName).AllParticipants.AveragedTrials.data(electrode,:);
                subjects    = fieldnames(AnalyzedEEG_Musicians_blsnrFFT.(condName));
            elseif strcmp(group,'NonMusicians')
                sFFT        = AnalyzedEEG_NonMusicians_blsnrFFT.(condName).AllParticipants.AveragedTrials.data(electrode,:);
                subjects    = fieldnames(AnalyzedEEG_NonMusicians_blsnrFFT.(condName));
            end
            
        end
        
        freqVec     = Log.(group).(subjects{1}).SpectralAnalysis.freqVec;
        
        % define range limits for the mean and std computation
        rangelimits = dsearchn(freqVec',[1.249 30]');
        range       = rangelimits(1):rangelimits(2);
        meanFFT     = mean(sFFT(range));
        stdFFT      = std(sFFT(range));
        
        % peak detection
        [pks, locs] = findpeaks (sFFT,freqVec,'MinPeakProminence',meanFFT+2*stdFFT,'MinPeakHeight',0+2*stdFFT); % 
        
        % storing                          
        EEG.peakDetection.(group).(condName).pks    = pks;
        EEG.peakDetection.(group).(condName).locs   = locs;
        
        subplot(4,2,iPlot)
        stem(freqVec,sFFT,'k','marker','none'); hold on
        scatter(locs,pks,'r')
        xlim([0 10])
        title(condName)
        box off
        
        
        iPlot = iPlot + 1;
    end
end

annotation('textbox',[0.06 0.95 0.49 0.05] ,'String','Musicians','EdgeColor','none',...
                 'HorizontalAlignment','center','fontsize',16,'BackgroundColor','w','FontWeight','bold') 
annotation('textbox',[0.5 0.95 0.49 0.05] ,'String','Non-musicians','EdgeColor','none',...
                 'HorizontalAlignment','center','fontsize',16,'BackgroundColor','w','FontWeight','bold') 


%% Merging of the frequencies of interest (FOI) of the stim and EEG

% Description
% -----------

% 1) retrieve and merge the prominent peaks from the stim and EEG
% 2) only keep FOI that are common in both groups

for iCond = 1:length(condNames)
    
    % retrieve the corresponding stimName
    condName = condNames{iCond};
    
    if iCond == 1 || iCond == 3;    stimName = stimNames{2};
    else;                           stimName = condName;
    end
        
    % check common FOI across groups
    tempFOI = zeros(300,2);
    
    for iGroup = 1:length(groups)
        group   = groups{iGroup};
        nLocs   = length(EEG.peakDetection.(group).(condName).locs);
        tempFOI(1:nLocs,iGroup) ...   
                = EEG.peakDetection.(group).(condName).locs;
    end
    locsEEG     = intersect(tempFOI(:,1),tempFOI(:,2)); 
    locsEEG(locsEEG ==0) = [];

    % to check which FOI was kept after merging the two grand averages
%     disp([condName,' & ',stimName]) 
%     [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locsEEG',frex);
%     disp(locsEEG(metRelIdx))
%     disp(length(metRelIdx))
%     disp(locsEEG(metUnrelIdx))
%     disp(length(metUnrelIdx))
    
    % retrieve the prominent peaks in the stim and EEG  (locs) + store
    locsStim                    = CochMod.(stimName).F3.AN.peaksDetection.locs;
    Zscores.AllFOI.(condName)   = unique(round([locsStim, locsEEG'],3));       
    testLocs = NaN(50,2);
    testLocs(1:length(locsStim),1) = locsStim;
    testLocs(1:length(locsEEG),2) = locsEEG';
    setdiff(round(locsEEG,2),round(locsStim,2)) % find added values from the EEG to the stim => in AAAA, only the 2.5Hz is added by the EEG, which explains why the main method, and the one based on the stimulus have exactly the same FOI.
    
    
    % test to check the number of FOI in the stim 
    testStimNlocs   = locsStim;
    idx2remove       = find(testStimNlocs>5.1);
    testStimNlocs(idx2remove) = [];

end


%% Mix Stim/EEG zscore with boundaries and a restriction of the # of FOI
%    a) [1.25 5]Hz
%    b) [1.25 30]Hz

% Description
% -----------

% 1) retrieve the frequencies of interest
% 2) set boundaries
% 3) retrieve the magnitude spectrum of the stim
% 4) classify metRel and metUnrel frequencies
% 5) restrict the number of FOI based on the highest magnitudes in the stim
% 6) retrieve the magnitude spectrum of the EEG
% 7) compute zscore

% parameters
    lowerLimit = [1.249    1.249];
    upperLimit = [30       5];
    
        
for iCond = 1:length(condNames)
    for zMethod = 1:2    
           
        % retrieve the frequencies of interest
        condName    = condNames{iCond}; 
        locs        = Zscores.AllFOI.(condName);

        % set boundaries 
        Options.lowerLimit  = lowerLimit(zMethod);
        Options.upperLimit  = upperLimit(zMethod);
        Options.remove5Hz   = true;
        [locs,Details]      = locsBoundaries(locs,Options); 
      
        % retrieve the magnitude spectrum of the cochlear model
        if iCond == 1 || iCond == 3;    stimName = stimNames{2};
        else;                           stimName = condName;
        end
        
        freqVec = CochMod.(stimName).F3.AN.freqVec;
        sFFT    = CochMod.(stimName).F3.AN.sFFT; % CochMod.(stimNames{iStim}).F3.AN.sFFT;
        
        % classify metRel and metUnrel frequencies
        [frexIdx,metRelIdx,metUnrelIdx] ...
                    = freqClassifier(freqVec,locs,frex);
        
        % manually select the lowest number of FOI across conditions
        fprintf('Group: %s, Condition: %s, Method: %i, nMetRel: %i, nMetUnrel: %i \n', group, condName, zMethod, length(metRelIdx), length(metUnrelIdx)) 
        
        % restrict the number of FOI (hardcoded)
        if      zMethod     == 1
                nMetRel     = 3; 
                nMetUnrel   = 24;
        elseif  zMethod     ==2
                nMetRel     = 3;
                nMetUnrel   = 12;          
        end
               
        % only keep the metUnrel FOI with the largest amplitudes in the stimulus (per condition)
        [~,maxIdxUnrel]  = maxk(sFFT(frexIdx(metUnrelIdx)),nMetUnrel);
        
        % only keep the metRel FOI with the largest amplitudes in the stimulus (per condition)
        % make sure it selects the three metRel frequencies below 5Hz
        metRelIdx2KeepForSure   = find(ismember(round(freqVec(frexIdx(metRelIdx)),2),[1.25, 2.5, 3.75]));
        remainingMetRelFreq     = 1:nMetRel;
        remainingMetRelFreq(metRelIdx2KeepForSure) = [];
        
        [~,maxIdxRel]  = maxk(sFFT(frexIdx(metRelIdx(remainingMetRelFreq))),nMetRel-3);

        % merge selected metRel and metUnrel FOI
        locs = sort([freqVec(frexIdx(metUnrelIdx(maxIdxUnrel))),...
                     freqVec(frexIdx(metRelIdx(metRelIdx2KeepForSure)))]);
    	
        if ~isempty(maxIdxRel)
            locs = sort([locs, freqVec(frexIdx(metRelIdx(maxIdxRel)))]);
        end
      

        for iGroup = 1:length(groups)
             group       = groups{iGroup};
              
            if exist('EEG.FFT')
                subjects    = EEG.FFT.(group).subjects;
            elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
                subjects    = fieldnames(AnalyzedEEG_Musicians_blsnrFFT.(condName));
            elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
                subjects    = fieldnames(AnalyzedEEG_NonMusicians_blsnrFFT.(condName));
            end
            
            idx2keep        = find(contains(subjects,'sub'));
            subjects        = subjects(idx2keep);
            
             for iSubject = 1:length(subjects)

                % load magnitude spectrum
                freqVec     = Log.(group).(subjects{iSubject}).SpectralAnalysis.freqVec;
                if exist('EEG.FFT')
                    sFFT        = EEG.FFT.(group).(condName).(subjects{iSubject}).AveragedTrials.data(electrode,:);
                elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
                    sFFT        = AnalyzedEEG_Musicians_blsnrFFT.(condName).(subjects{iSubject}).AveragedTrials.data(electrode,:);
                elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
                    sFFT        = AnalyzedEEG_NonMusicians_blsnrFFT.(condName).(subjects{iSubject}).AveragedTrials.data(electrode,:);
                end                
                
                
                % classify metRel and metUnrel frequencies    
                [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);

                % computation of the zscore 
                amps                = sFFT(frexIdx);
                zscores             = (amps - mean(amps)) ./ std(amps);
                zscoresMetRel       = mean(zscores(metRelIdx));
                zscoresMetUnrel     = mean(zscores(metUnrelIdx));

                Zscores.(methodNames{zMethod}).(group).(condName).amps(iSubject,:)          = amps;
                Zscores.(methodNames{zMethod}).(group).(condName).locs(iSubject,:)          = locs;
                Zscores.(methodNames{zMethod}).(group).(condName).frexIdx(iSubject,:)       = frexIdx;
                Zscores.(methodNames{zMethod}).(group).(condName).metRelIdx(iSubject,:)     = metRelIdx;
                Zscores.(methodNames{zMethod}).(group).(condName).metUnrelIdx(iSubject,:)   = metUnrelIdx;
                Zscores.(methodNames{zMethod}).(group).(condName).zMetRel(iSubject,:)       = zscoresMetRel;
                Zscores.(methodNames{zMethod}).(group).(condName).zMetUnrel(iSubject,:)     = zscoresMetUnrel;
             end
        end

        clear Options zscores freq locs freqIdx metRelIdx metUnrelIdx amps
    end
end



%% Zscore of the stimulus for the two first methods

% Description
% -----------

% 1) retrieve the frequencies of interest
% 2) retrieve the magnitude spectrum of the stim
% 3) classify metRel and metUnrel frequencies
% 4) compute zscore

for iCond = 1:length(condNames)
    for zMethod = 1:2
        
        % test
        locsMus        = Zscores.(methodNames{zMethod}).Musicians.(condName).locs(1,:);
        locsNonMus     = Zscores.(methodNames{zMethod}).NonMusicians.(condName).locs(1,:);
        %isequal(locsMus,locsNonMus)
        
        % retrieve FOI
        condName    = condNames{iCond};
        locs        = Zscores.(methodNames{zMethod}).(group).(condName).locs(1,:);
        
        % retrieve the magnitude spectrum of the cochlear model
        if iCond == 1 || iCond == 3;    stimName = stimNames{2};
        else;                           stimName = condName;
        end
        
        freqVec = CochMod.(stimName).F3.AN.freqVec;
        sFFT    = CochMod.(stimName).F3.AN.sFFT;
        
        % classify metRel and metUnrel frequencies
        [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);
        
        % computation of the zscore 
        amps                = sFFT(frexIdx);
        % if iCond == 1 || iCond == 3; disp(amps); end 
        zscores             = (amps - mean(amps)) ./ std(amps);
        zscoresMetRel       = mean(zscores(metRelIdx));
        zscoresMetUnrel     = mean(zscores(metUnrelIdx));        
        
        % store
        Zscores.(methodNames{zMethod}).cochMod.(stimName).locs          = locs;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).amps          = amps;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).frexIdx       = frexIdx;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).metRelIdx     = metRelIdx;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).metUnrelIdx   = metUnrelIdx;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).zMetRel       = zscoresMetRel;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).zMetUnrel     = zscoresMetUnrel;   
        
    end
end

%% Stim based zscore with [1.25 30]Hz boundaries and a restriction of the # of FOI

% Description
% -----------

% 1) retrieve the frequencies of interest in the stim for each condition
% 2) set boundaries
% 3) find and store the number of met(un)rel FOI for each condition
% 4) find the lowest number of FOI across condition = n
% 5) find the n highest peaks in the stim magnitude spectrum
% 6) compute the stim zscore
% 7) compute the EEG zscore

zMethod = 3;

for iCond = 1:length(condNames)

    % retrieve the corresponding stimName
    condName = condNames{iCond};
    
    if iCond == 1 || iCond == 3;    stimName = stimNames{2};
    else;                           stimName = condName;
    end
       
    % find prominent peaks in the stimulus
    nLocs   = length(CochMod.(stimName).F3.AN.peaksDetection.locs);
    locs    = CochMod.(stimName).F3.AN.peaksDetection.locs; 
    pks     = CochMod.(stimName).F3.AN.peaksDetection.pks;
    
    % set boundaries 
    Options.lowerLimit  = 1.249; % Hz
    Options.upperLimit  = 30;
    Options.remove5Hz   = true;
    Options.pks         = pks;
    [locs,Details,pks]  = locsBoundaries(locs,Options); 
    
    % classify metRel and metUnrel frequencies    
    freqVec             = CochMod.(stimName).F3.AN.freqVec;
    [frexIdx,metRelIdx,metUnrelIdx] ...
                        = freqClassifier(freqVec,locs,frex); 
        
    % find and store the number of met(un)rel FOI for each condition
    nMetRel(iCond)      = length(metRelIdx);
    nMetUnrel(iCond)    = length(metUnrelIdx); 
    

    subplot(1,4,iCond)
    stem(freqVec, CochMod.(stimName).F3.AN.sFFT,'k','marker','none'); hold on
    scatter(locs,pks,'r')
    xlim([0 15])
    clear locs Options
end

% find the lowest number of FOI in each category
nMetRel     = min(nMetRel);
nMetUnrel   = min(nMetUnrel);

for iCond = 1:length(condNames)
    
    % retrieve the corresponding stimName
    condName = condNames{iCond};
    
    if iCond == 1 || iCond == 3;    stimName = stimNames{2};
    else;                           stimName = condName;
    end
    
    % load cochlear model
    freqVec         = CochMod.(stimName).F3.AN.freqVec;
    sFFT            = CochMod.(stimName).F3.AN.sFFT; 
    locs            = CochMod.(stimName).F3.AN.peaksDetection.locs;
    
    % set boundaries & categories FOI
    Options.lowerLimit  = 1.249; % Hz
    Options.upperLimit  = 30;
    Options.remove5Hz   = true;
    [locs,Details]      = locsBoundaries(locs,Options); 
    
    [frexIdx,metRelIdx,metUnrelIdx] ...
                    = freqClassifier(freqVec,locs,frex); 
           
    % only keep the metUnrel FOI with the largest amplitudes in the stimulus (per condition)
    [~,maxIdxUnrel]  = maxk(sFFT(frexIdx(metUnrelIdx)),nMetUnrel);
        
    % only keep the metRel FOI with the largest amplitudes in the stimulus (per condition)
    % make sure it selects the three metRel frequencies below 5Hz
    if nMetRel > 3 
        metRelIdx2KeepForSure       = find(ismember(round(freqVec(frexIdx(metRelIdx)),2),[1.25, 2.5, 3.75]));
        additionnalMetRel2select    = 1:nMetRel;
        additionnalMetRel2select(metRelIdx2KeepForSure) = [];
        
        [~,maxIdxRel]  = maxk(sFFT(frexIdx(metRelIdx(remainingMetRelFreq))),nMetRel-3);
        
        % merge selected metRel and metUnrel FOI
        locs = sort([freqVec(frexIdx(metUnrelIdx(maxIdxUnrel))),...
                     freqVec(frexIdx(metRelIdx(metRelIdx2KeepForSure))),...
                     freqVec(frexIdx(metRelIdx(maxIdxRel)))]);   
    else
        frexIdx_metRelIdx   = dsearchn(freqVec',[1.25, 2.5, 3.75]');
        locs                = sort([freqVec(frexIdx(metUnrelIdx(maxIdxUnrel))),...
                                    freqVec(frexIdx_metRelIdx)]);
    end
    
    % compute the stim zscore
    % -----------------------

        % classify metRel and metUnrel frequencies    
        [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);

        % computation of the zscore 
        amps                = sFFT(frexIdx);
        zscores             = (amps - mean(amps)) ./ std(amps);
        zscoresMetRel       = mean(zscores(metRelIdx));
        zscoresMetUnrel     = mean(zscores(metUnrelIdx));

        Zscores.(methodNames{zMethod}).cochMod.(stimName).locs          = locs;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).amps          = amps;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).frexIdx       = frexIdx;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).metRelIdx     = metRelIdx;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).metUnrelIdx   = metUnrelIdx;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).zMetRel       = zscoresMetRel;
        Zscores.(methodNames{zMethod}).cochMod.(stimName).zMetUnrel     = zscoresMetUnrel;         
        

    % compute the EEG zscore
    % ----------------------
    
    for iGroup = 1:length(groups)
         group       = groups{iGroup};
        if exist('EEG.FFT')
            subjects    = EEG.FFT.(group).subjects;
        elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
            subjects    = fieldnames(AnalyzedEEG_Musicians_blsnrFFT.(condName));
        elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
            subjects    = fieldnames(AnalyzedEEG_NonMusicians_blsnrFFT.(condName));
        end

        idx2keep        = find(contains(subjects,'sub'));
        subjects        = subjects(idx2keep);
    
         for iSubject = 1:length(subjects)

                % load magnitude spectrum
                freqVec     = Log.(group).(subjects{iSubject}).SpectralAnalysis.freqVec;
                if exist('EEG.FFT')
                    sFFT        = EEG.FFT.(group).(condName).(subjects{iSubject}).AveragedTrials.data(electrode,:);
                elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
                    sFFT        = AnalyzedEEG_Musicians_blsnrFFT.(condName).(subjects{iSubject}).AveragedTrials.data(electrode,:);
                elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
                    sFFT        = AnalyzedEEG_NonMusicians_blsnrFFT.(condName).(subjects{iSubject}).AveragedTrials.data(electrode,:);
                end          
                
                % load the locs from the stim based method
                locs        = Zscores.(methodNames{zMethod}).cochMod.(stimName).locs;
                
                % classify metRel and metUnrel frequencies    
                [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex);

                % computation of the zscore 
                amps                = sFFT(frexIdx);
                zscores             = (amps - mean(amps)) ./ std(amps);
                zscoresMetRel       = mean(zscores(metRelIdx));
                zscoresMetUnrel     = mean(zscores(metUnrelIdx));

                Zscores.(methodNames{zMethod}).(group).(condName).amps(iSubject,:)          = amps;
                Zscores.(methodNames{zMethod}).(group).(condName).locs(iSubject,:)          = locs;
                Zscores.(methodNames{zMethod}).(group).(condName).frexIdx(iSubject,:)       = frexIdx;
                Zscores.(methodNames{zMethod}).(group).(condName).metRelIdx(iSubject,:)     = metRelIdx;
                Zscores.(methodNames{zMethod}).(group).(condName).metUnrelIdx(iSubject,:)   = metUnrelIdx;
                Zscores.(methodNames{zMethod}).(group).(condName).zMetRel(iSubject,:)       = zscoresMetRel;
                Zscores.(methodNames{zMethod}).(group).(condName).zMetUnrel(iSubject,:)     = zscoresMetUnrel;
         end
    end
        disp([methodNames{zMethod},' ',condName,' ',group])
        disp(locs(metRelIdx)')
        disp(length(metRelIdx))
        disp(locs(metUnrelIdx)')
        disp(length(metUnrelIdx))
end


%% Figures for method comparison (per group, condition, conditionOrder)

figure('Color', [1 1 1], 'Position',[1 1 1200 700], 'name', 'MetRelZscores')
iPlot   = 1; 
nMethod = 2;

for zMethod = 1:nMethod
    for iGroup = 1:length(groups)
                    
        group       = groups{iGroup};
        if exist('EEG.FFT')
            subjects    = EEG.FFT.(group).subjects;
        elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
            subjects    = fieldnames(AnalyzedEEG_Musicians_blsnrFFT.(condName));
        elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
            subjects    = fieldnames(AnalyzedEEG_NonMusicians_blsnrFFT.(condName));
        end
        idx2keep        = find(contains(subjects,'sub'));
        subjects        = subjects(idx2keep);

        % find the subjectIdx corresponding to each of the condition order
        i=1; j=1;
        for iSubject = 1:length(subjects)
            if strcmp(Log.(group).(subjects{iSubject}).condOrder{1},'ABAB1') 
                subIdx2plot(i,1)    = iSubject;    
                condOrderIdx(i,1)   = str2double(strrep(subjects{iSubject},'sub0','')); 
                i = i+1;
            elseif strcmp(Log.(group).(subjects{iSubject}).condOrder{1},'CDEF')   
                subIdx2plot(j,2)    = iSubject;
                condOrderIdx(j,2)   = str2double(strrep(subjects{iSubject},'sub0',''));  
                j = j+1;
            end  
        end
        

        %%%%%%%%%%%
        % scatter %
        %%%%%%%%%%%
        subplot(nMethod,length(groups),iPlot)
        iCondPlot = 1;
        for iCond = [2,1,3,4]
            for iCondOrder = 1:2
                
                condName            = condNames{iCond};
                data2plot           = Zscores.(methodNames{zMethod}).(group).(condName).zMetRel(subIdx2plot(:,iCondOrder),:);
                condOrder_shifts    = [-0.1 0.1];
                colors              = {[0.8500 0.3250 0.0980],[0 0.4470 0.7410]};
                
                scatter(repmat(iCondPlot + condOrder_shifts(iCondOrder),length(data2plot),1), ...
                        data2plot, ...
                        'markerEdgeColor', colors{iCondOrder}, 'markerEdgeAlpha',0.8, ...
                        'markerFaceColor', colors{iCondOrder}, 'markerEdgeAlpha',0.6) %, ... 'jitter','on', 'jitterAmount',0.15);
                         hold on
                 
                 clear data2plot  
            end
            
            iCondPlot = iCondPlot + 1;
            
            % retrieve stim zscores
            if iCond == 1 || iCond == 3;    stimName = stimNames{2};
            else;                           stimName = condName;
            end
            stim2plot(iCond) = Zscores.(methodNames{zMethod}).cochMod.(stimName).zMetRel;            
        end
        
        % plot stim
        scatter([2,1,3,4], stim2plot,'k')
        
        % layout
        set(gca,'XTick',[1:1:4],'XTickLabel',{'AAAA','ABAB1','ABAB2','CDEF'},'TickDir','out','LineWidth',1.25,'fontsize',14,'XLim',[0.6 4.4],'YLim',[-1.3 2.1],'YTick',[-1:0.5:2]) %,'YTick',[-0.5:0.25:0.5],'YLim',[-0.65 0.3] 
        box off
        ylabel('FFT Meter-related Zscores','fontsize',16)
        line([0 10],[0 0],'LineStyle','--','Color','k','LineWidth',0.5)
        title ([methodNames{zMethod},' - ', group])

        iPlot = iPlot + 1;
    end
end


if printFig
%     cd(fullfile(projectPath,'4_Figures'))
%     set(gcf,'PaperPositionMode','auto')
%     print(['Fig1 - Method comparison'],'-djpeg','-r800')
%     cd(fullfile(projectPath,'3_nalysis'))
end


%% Figures Spectral domain visualization

% Main parameters
zMethod         = 1;
fullCondNames   = {'Long Pattern Repetition #1', ...
                   'Medium Pattern Repetition', ...
                   'Long Pattern Repetition #2', ...
                   'No Pattern Repetition'};
fontsize = 17;   %  Cfg.Figure.fontsize ; 
linewidth = 1.3; %Cfg.Figure.linewidth-0.5; 


% Figure preparation
figure('Units','centimeters','Position',[1 1 40 20],'Color',[1 1 1])
plotPos = {[0.30 0.62 0.20 0.23],[0.30 0.36 0.20 0.23],[0.30 0.1 0.20 0.23];...
           [0.05 0.62 0.20 0.23],[0.05 0.36 0.20 0.23],[0.05 0.1 0.20 0.23];...
           [0.55 0.62 0.20 0.23],[0.55 0.36 0.20 0.23],[0.55 0.1 0.20 0.23];...
           [0.80 0.62 0.20 0.23],[0.80 0.36 0.20 0.23],[0.80 0.1 0.20 0.23]};
       
plotPos = {[0.35 0.60 0.18 0.22],[0.35 0.35 0.18 0.22],[0.35 0.1 0.18 0.22];...
           [0.12 0.60 0.18 0.22],[0.12 0.35 0.18 0.22],[0.12 0.1 0.18 0.22];...
           [0.58 0.60 0.18 0.22],[0.58 0.35 0.18 0.22],[0.58 0.1 0.18 0.22];...
           [0.81 0.60 0.18 0.22],[0.81 0.35 0.18 0.22],[0.81 0.1 0.18 0.22]};       

xLimMax = 10.5;

for iCond = [2,1,3,4]
    
    condName = condNames{iCond};
    if iCond == 1 || iCond == 3;    stimName = stimNames{2};
    else;                           stimName = condName;
    end
    
    clear frexIdx metRelIdx metUnrelIdx
       
    % stim 
    % ---- 
        axes('Position',plotPos{iCond,1})
        freqVec     = CochMod.(stimName).F3.AN.freqVec;
        sFFT        = CochMod.(stimName).F3.AN.sFFT;
        locs        = Zscores.(methodNames{zMethod}).cochMod.(stimName).locs;
        amps        = Zscores.(methodNames{zMethod}).cochMod.(stimName).amps;
        metRelIdx   = Zscores.(methodNames{zMethod}).cochMod.(stimName).metRelIdx;
        metUnrelIdx = Zscores.(methodNames{zMethod}).cochMod.(stimName).metUnrelIdx;

        stem1 = stem(freqVec,sFFT,'Color',[0.6 0.6 0.6],'Marker','none');
        hold on
        stem2 = stem(locs(metRelIdx),    amps(metRelIdx),    'Color',[Cfg.Figure.metRelColor],'Marker','none','LineWidth',Cfg.Figure.linewidth+0.25);  % Colors{iCond,3} Cfg.Figure.linewidth+1.5
        stem3 = stem(locs(metUnrelIdx),  amps(metUnrelIdx),  'Color',[Cfg.Figure.metUnrelColor],'Marker','none','LineWidth',Cfg.Figure.linewidth);

        % plot layout
        box off
        set(gca,'tickDir','out','LineWidth',linewidth,'fontsize',fontsize)
        set(gca,'xlim',[0 xLimMax],'xticklabel',{[]})
        line([0 0],[-1 5],'color','k','LineWidth',linewidth)
        
        if iCond == 1
            set(gca,'ylim',[-0.1 3.25],'yTick',[0:0.5:3],'yticklabel',{'0','','','','','','3'}) % set(gca,'yticklabel',{[]})         
        elseif iCond == 2
            ylabel('a.u.','rotation',90,'Position',[-0.45 mean([0 2.75]) 0])
            set(gca,'ylim',[-0.1 2.75],'YTick',[0:0.5:2.5],'yticklabel',{'0','','','','','2.5'})
        elseif iCond == 3
            set(gca,'ylim',[-0.1 3.25],'yTick',[0:0.5:3],'yticklabel',{'0','','','','','','3'})
        elseif iCond == 4
            set(gca,'ylim',[-0.1 1.8],'YTick',[0:0.5:1.5],'yticklabel',{'0','','','1.5'})
        end

        clear frexIdx metRelIdx metUnrelIdx    
    
    % EEG    
    for iGroup = 1:length(groups)

        % subplot(3,4,(iGroup*4)+iStim)
        axes('Position',plotPos{iCond,1+iGroup})
        
        % load magnitude spectrum
        group       = groups{iGroup};
        freqVec     = Log.AllParticipants.SpectralAnalysis.freqVec;
        if exist('EEG.FFT')
            sFFT        = EEG.FFT.(group).(condName).AllParticipants.AveragedTrials.data(electrode,:);
        elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
            sFFT        = AnalyzedEEG_Musicians_blsnrFFT.(condName).AllParticipants.AveragedTrials.data(electrode,:);
        elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
            sFFT        = AnalyzedEEG_NonMusicians_blsnrFFT.(condName).AllParticipants.AveragedTrials.data(electrode,:);
        end  
        frexIdx     = Zscores.(methodNames{zMethod}).(group).(condName).frexIdx(1,:);
        metRelIdx   = Zscores.(methodNames{zMethod}).(group).(condName).metRelIdx(1,:);
        metUnrelIdx = Zscores.(methodNames{zMethod}).(group).(condName).metUnrelIdx(1,:);
        
        stem(freqVec,sFFT,'Color',[0.6 0.6 0.6],'Marker','none');%,'LineWidth',Cfg.Figure.linewidth-0.5); 
        hold on
        stem2 = stem(freqVec(frexIdx(metRelIdx)),   sFFT(frexIdx(metRelIdx)),   'Color',[Cfg.Figure.metRelColor],'Marker','none','LineWidth',linewidth+0.25);
        stem3 = stem(freqVec(frexIdx(metUnrelIdx)), sFFT(frexIdx(metUnrelIdx)), 'Color',[Cfg.Figure.metUnrelColor],'Marker','none','LineWidth',linewidth);

        % plot layout
        box off
        set(gca,'tickDir','out','LineWidth',linewidth,'fontsize',fontsize)
        set(gca,'xlim',[0 xLimMax],'xticklabel',{[]})
        line([0 0],[-1 5],'color','k','LineWidth',linewidth)
        
        % xlabel 
        if iGroup == 2
            set(gca,'xlim',[0 xLimMax],'xTick',[0:2.5:10],'xticklabel',{'0','2.5','5','7.5','10'})
            xlabel('Frequencies (Hz)','fontsize',fontsize)
        else
            set(gca,'xlim',[0 xLimMax])
        end
        
        % ylabel
        if iCond == 4
            set(gca,'ylim',[-0.02 0.17],'yTick',[0:0.05:0.15],'yticklabel',{'0','','','0.15'})
        elseif iCond == 2
            set(gca,'ylim',[-0.02 0.36],'yTick',[0:0.1:0.35],'yticklabel',{'0','','','0.3'})
            ylabel('µV','rotation',90,'Position',[-0.45 mean([0 0.36]) 0])
        else
            set(gca,'ylim',[-0.02 0.36],'yTick',[0:0.1:0.35],'yticklabel',{'0','','','0.3'})
        end
        
        % legend
        if iCond == 4 && iGroup == 2
            legend([stem2 stem3],{['Meter-related' newline 'frequencies'],['Meter-unrelated' newline 'frequencies']},...
                'Position',[0.849393738977072,0.728574530500828,0.165972222222222,0.126102292768959],'fontsize',fontsize-3)
            legend('boxoff')
        end
    end
    
    % condition titles
    annotation('textbox',[plotPos{iCond,1}(1), 0.9, plotPos{iCond,1}(3),0.1] ,'String',fullCondNames{iCond},'EdgeColor','none',...
                     'HorizontalAlignment','center','fontsize',fontsize+3,'BackgroundColor','w','FontWeight','bold') % 
            
    iStim = iStim +1;
    
end

annotation('textbox',[0.01, plotPos{4,1}(2)+plotPos{4,1}(4)/4, 0.07, 0.1] ,'String','Cochlear Model','EdgeColor','none',...
                 'HorizontalAlignment','center','fontsize',fontsize,'BackgroundColor','w','FontWeight','bold') % ,'FontWeight','bold'
annotation('textbox',[0.01, plotPos{4,2}(2)+plotPos{4,2}(4)/4, 0.07, 0.1] ,'String','Musicians','EdgeColor','none',...
                 'HorizontalAlignment','center','fontsize',fontsize,'BackgroundColor','w','FontWeight','bold') % ,'FontWeight','bold'
annotation('textbox',[0.01, plotPos{4,3}(2)+plotPos{4,3}(4)/4, 0.07, 0.1] ,'String','Non- musicians','EdgeColor','none',...
                 'HorizontalAlignment','center','fontsize',fontsize,'BackgroundColor','w','FontWeight','bold') 


if printFig
    cd(fullfile(projectPath,'4_Figures'))
     set(gcf,'PaperPositionMode','auto')
    print(['Fig4 - Spectral Domain visualization'],'-djpeg','-r800')
    cd(fullfile(projectPath,'3_analysis'))
end

%% Figures zscore visualization

% Figure preparation
figure('Units','centimeters','Position',[1 1 40 17],'Color',[1 1 1])
fontsize    = 17;   
linewidth   = 1.3;
markers     = {'o','square'};
axPos       = {[0.1 0.15 0.4 0.75],[0.55 0.15 0.4 0.75]};
Colors      = {[214 115 47 ]./256  ,[246 219 196]./256; ...
               [31  47  55 ]./256  ,[198 216 230]./256; ...
               [80  155 143]./256  ,[185 215 210]./256};

for iGroup = 1:length(groups)
         
    group       = groups{iGroup};
    if exist('EEG.FFT')
        subjects    = EEG.FFT.(group).subjects;
    elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
        subjects    = fieldnames(AnalyzedEEG_Musicians_blsnrFFT.(condName));
    elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
        subjects    = fieldnames(AnalyzedEEG_NonMusicians_blsnrFFT.(condName));
    end
    idx2keep        = find(contains(subjects,'sub'));
    subjects        = subjects(idx2keep);    

    % find the subjectIdx corresponding to each of the condition order
    i=1; j=1;
    for iSubject = 1:length(subjects)
        if strcmp(Log.(group).(subjects{iSubject}).condOrder{1},'ABAB1') 
            subIdx2plot(i,1)    = iSubject;    
            condOrderIdx(i,1)   = str2double(strrep(subjects{iSubject},'sub0','')); 
            i = i+1;
        elseif strcmp(Log.(group).(subjects{iSubject}).condOrder{1},'CDEF')   
            subIdx2plot(j,2)    = iSubject;
            condOrderIdx(j,2)   = str2double(strrep(subjects{iSubject},'sub0',''));  
            j = j+1;
        end  
    end        
       
    axes('Position', axPos{iGroup})
    
    iCondPlot = 1; iBox = 1; ixPos = 1;
    for iCond = [2,1,3,4]
        
        condName = condNames{iCond};
        if iCond == 1 || iCond == 3;    stimName = stimNames{2}; iStim = 2;
        elseif iCond == 2;              stimName = condName;     iStim = 1;
        elseif iCond == 4;              stimName = condName;     iStim = 3;
        end
         
        for iCondOrder = 1:2

            data2plot           = Zscores.(methodNames{zMethod}).(group).(condName).zMetRel(subIdx2plot(:,iCondOrder),:);
            condOrder_shifts    = [-0.13 0.13];

            % scatter to initialise the legend in black
            scatter(1, -5, 'marker',markers{iCondOrder}, 'markerEdgeColor', [0 0 0], 'DisplayName', condOrderNames{iCondOrder}); hold on % 'DisplayName', ['Condition order #', num2str(iCondOrder)]
            
            % scatter   
            scatter(repmat(iCondPlot + condOrder_shifts(iCondOrder),length(data2plot),1), ...
                        data2plot, ...
                        'SizeData',100,...
                        'marker',markers{iCondOrder},...
                        'markerEdgeColor', Colors{iStim,1}, 'markerEdgeAlpha',0.8, ...
                        'markerFaceColor', Colors{iStim,2}, 'markerEdgeAlpha',0.6, ...
                        'jitter','on', 'jitterAmount',0.07,...
                        'DisplayName', condOrderNames{iCondOrder}) 

            % boxplot preparation
            boxZscore(iBox:iBox+length(data2plot)-1) = data2plot;
            boxCond(iBox:iBox+length(data2plot)-1)   = repmat(iCondPlot + condOrder_shifts(iCondOrder),length(data2plot),1);
            xPos(ixPos) = iCondPlot + condOrder_shifts(iCondOrder);

            iBox = iBox+length(data2plot);
            ixPos = ixPos+1;
            clear data2plot 
        end
           
        % retrieve stim zscores
        stim2plot(iCond) = Zscores.(methodNames{zMethod}).cochMod.(stimName).zMetRel;
        line([iCondPlot-0.25 iCondPlot+0.25],[stim2plot(iCond) stim2plot(iCond)], ...
             'LineWidth',linewidth+1.5, 'Color', [0.2 0.2 0.2])
        
        iCondPlot = iCondPlot + 1;
    end
    
    h = boxplot(boxZscore,boxCond,'Colors','k','Widths',0.1,'Positions',xPos,'symbol','','Labels','','Whisker',0, 'symbol', '');    % ,'Notch','on'
        set(h,{'linew'},{linewidth-0.25});
    
    % https://nl.mathworks.com/matlabcentral/answers/101922-how-do-i-create-a-multi-line-tick-label-for-a-figure-using-matlab-7-10-r2010a
    xlabel_row1 = {'Medium Pattern','Long Pattern', 'Long Pattern', 'No Pattern'};
    xlabel_row2 = {'Repetition',' Repetition #1', ' Repetition #2', 'Repetition'};
    labelArray  = [xlabel_row1; xlabel_row2];
    labelArray  = strjust(pad(labelArray),'left'); % 'left'(default)|'right'|'center
    tickLabels  = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    ax          = gca();
    ax.TickLabelInterpreter = 'tex';  % needed for some plots like boxplot.
    ax.XTick        = 1:4; 
    ax.XLim         = [0.6 4.4];
    ax.XTickLabel   = tickLabels; 
    ax.FontSize     = 16;
    % ax.FontWeight   ='bold';
    xtickangle(45)
    
    if iGroup == 2
        lh = legend(legendUnq());
        lh.AutoUpdate   = 'off';
        lh.Box          = 'off';
        lh.Position     = [0.8273,0.8342151675485,0.14365837191358,0.075220458553792];
        annotation('textbox', [0.82,0.9094,0.18,0.04], 'String','Condition order offering:',...
                   'fontsize',lh.FontSize, 'EdgeColor','none', 'FontWeight',lh.FontWeight, 'HorizontalAlignment','left');
    end
    
    set(gca,'XTick',[1:1:4],'TickDir','out','LineWidth',1.25,'fontsize',16,'XLim',[0.6 4.4],'YLim',[-1.3 2.1],'YTick',[-1:0.5:2]) %,'YTick',[-0.5:0.25:0.5],'YLim',[-0.65 0.3] 
    box off
    line([0 10],[0 0],'LineStyle','--','Color','k','LineWidth',0.5)
    title (groupNames{iGroup}) 
    
    if iGroup == 1
        ylabel('Meter-related Zscores','fontsize',16)
    end     
end



if printFig
    cd(fullfile(projectPath,'4_Figures'))
    set(gcf,'PaperPositionMode','auto')
    print(['Fig5 - Zscores'],'-djpeg','-r800')
    cd(fullfile(projectPath,'3_analysis'))
end
        
%% Prepare table for RStudio

% preallocate
i = 1;
Subject=[];

load(fullfile(projectPath,'Log.mat'))

for iGroup = 1:length(groups)
    for iCond = 1:length(condNames)
        
        group       = groups{iGroup};
        condName    = condNames{iCond};
        if exist('EEG.FFT')
            subjects    = EEG.FFT.(group).subjects;
        elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
            subjects    = fieldnames(AnalyzedEEG_Musicians_blsnrFFT.(condName));
        elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
            subjects    = fieldnames(AnalyzedEEG_NonMusicians_blsnrFFT.(condName));
        end
        idx2keep        = find(contains(subjects,'sub'));
        subjects        = subjects(idx2keep);         
        
        
        if iCond == 1 || iCond == 3;    stimName = stimNames{2};
        else;                           stimName = condName;
        end
    
        for iSubject = 1:length(subjects)

                Subject{i,1}        = [group,'_',subjects{iSubject}];
                Group{i,1}          = group;
                Condition{i,1}      = condName;
                
                condOrderName       = Log.(group).(subjects{iSubject}).condOrder{1};
                if strcmp(condOrderName,'ABAB1'); condOrderName = strrep(condOrderName,'1',''); end
                
                CondOrder{i,1}      = [condOrderName,' first'];
                MetRelzScore_MixLargeRange(i,1) ...
                                    = Zscores.(methodNames{1}).(group).(condName).zMetRel(iSubject,:);
                SelectiveEnhancement_MixLargeRange(i,1) ...
                                    = Zscores.(methodNames{1}).(group).(condName).zMetRel(iSubject,:)...
                                    - Zscores.(methodNames{1}).cochMod.(stimName).zMetRel;
                MetRelzScore_MixShortRange(i,1) ...
                                    = Zscores.(methodNames{2}).(group).(condName).zMetRel(iSubject,:);
                SelectiveEnhancement_MixShortRange(i,1) ...
                                    = Zscores.(methodNames{2}).(group).(condName).zMetRel(iSubject,:) ...
                                    - Zscores.(methodNames{2}).cochMod.(stimName).zMetRel;
                MetRelzScore_StimBased(i,1) ...
                                    = Zscores.(methodNames{3}).(group).(condName).zMetRel(iSubject,:);
                SelectiveEnhancement_StimBased(i,1) ...
                                    = Zscores.(methodNames{3}).(group).(condName).zMetRel(iSubject,:) ...
                                    - Zscores.(methodNames{3}).cochMod.(stimName).zMetRel;                                
                OverallGain(i,1)    = sum(Zscores.(methodNames{1}).(group).(condName).amps(iSubject,:));
                i = i+1;
        end
    end
end

T = table(Subject, Group, Condition, CondOrder, MetRelzScore_MixLargeRange, SelectiveEnhancement_MixLargeRange,...
                                                MetRelzScore_MixShortRange, SelectiveEnhancement_MixShortRange,...
                                                MetRelzScore_StimBased,     SelectiveEnhancement_StimBased,...
                                                OverallGain);
T = sortrows(T,{'Group','Subject'});

cd(fullfile(projectPath,'3_analysis'))
writetable(T,'StatisticalAnalysisEEG.xls')



%% visualization of the main effect of condition and condition orders

% Main parameters
zMethod     = 1;
fontsize    = 22;   
linewidth   = 1.8; 
plotEEG     = true;

clear data2plot xPos boxCond boxZscore

% Conditions
% ---------------

% Figure preparation
figure('Units','centimeters','Position',[1 1 40 20],'Color',[1 1 1])
iPlot = 1;      
iBox  = 1;

for iCond = [2,1,3,4]
    
    condName    = condNames{iCond};
    i           = 1;
    if iCond == 1 || iCond == 3;    stimName = stimNames{2};    iStim = 2;
    elseif iCond ==2;               stimName = condName;        iStim = 1;
    elseif iCond ==4;               stimName = condName;        iStim = 3;
    end
    
    clear frexIdx metRelIdx metUnrelIdx
    
    for iGroup = 1:length(groups)
        group       = groups{iGroup};
        if exist('EEG.FFT')
            subjects    = EEG.FFT.(group).subjects;
        elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
            subjects    = fieldnames(AnalyzedEEG_Musicians_blsnrFFT.(condName));
        elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
            subjects    = fieldnames(AnalyzedEEG_NonMusicians_blsnrFFT.(condName));
        end
        idx2keep        = find(contains(subjects,'sub'));
        subjects        = subjects(idx2keep); 

        for iSubject = 1:length(subjects)
                 
            data2plot(i,iCond)  = Zscores.(methodNames{zMethod}).(group).(condName).zMetRel(iSubject);
            i                   = i+1;
        end    
    end

    if plotEEG
        scatter(repmat(iPlot,length(data2plot),1),data2plot(:,iCond), ...
            'jitter','on', 'jitterAmount',0.07,...
            'SizeData',100,...
            'markerEdgeColor', Colors{iStim,1}, 'markerEdgeAlpha',0.8, ...
            'markerFaceColor', Colors{iStim,2}, 'markerEdgeAlpha',0.6); hold on
    end
    
    % stim zscore           
    stim2plot(iCond) = Zscores.(methodNames{zMethod}).cochMod.(stimName).zMetRel; 
    line([iPlot-0.25 iPlot+0.25],[stim2plot(iCond) stim2plot(iCond)],...
        'LineWidth',linewidth+1.5, 'Color', [0.2 0.2 0.2])
    
    % boxplot preparation 
    boxZscore(iBox:iBox+length(data2plot)-1) = data2plot(:,iCond);
    boxCond(iBox:iBox+length(data2plot)-1)   = repmat(iPlot,length(data2plot),1);
    xPos(iPlot)                              = iPlot;
    
    iBox = iBox+length(data2plot);
    iPlot = iPlot + 1;
end
if plotEEG 
    h = boxplot(boxZscore,boxCond,'Colors','k','Widths',0.16,'Positions',xPos,'symbol','','Labels','','Whisker',0, 'symbol', '');    % ,'Notch','on'
        set(h,{'linew'},{linewidth-0.25});
end

% statistical significance
% sigline([1,2],'p = 0.006',[],1.45)
% sigline([1,2],[],1.45)

% layout
set(gca,'XTick',[1:1:4],'TickDir','out','LineWidth',linewidth+0.25,'fontsize',fontsize,'XLim',[0.6 4.4],'YLim',[-1.25 1.6],'YTick',[-1:0.5:1.5]) %,'YTick',[-0.5:0.25:0.5],'YLim',[-0.65 0.3] 
box off
ylabel('Beat-Related Zscores','fontsize',16)
line([0 10],[0 0],'LineStyle','--','Color','k','LineWidth',0.5)

% xlabel 
    % https://nl.mathworks.com/matlabcentral/answers/101922-how-do-i-create-a-multi-line-tick-label-for-a-figure-using-matlab-7-10-r2010a
    xlabel_row1 = {'One Pattern','Two Patterns', 'Two Patterns', 'No Pattern'};
    xlabel_row2 = {'Repetition',' Repetition #1', ' Repetition #2', 'Repetition'};
    labelArray  = [xlabel_row1; xlabel_row2];
    labelArray  = strjust(pad(labelArray),'left'); % 'left'(default)|'right'|'center
    tickLabels  = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    ax          = gca();
    ax.TickLabelInterpreter = 'tex';  % needed for some plots like boxplot.
    ax.XTick        = 1:4; 
    ax.XLim         = [0.6 4.4];
    ax.XTickLabel   = tickLabels; 
    ax.FontSize     = fontsize;
    ax.FontWeight   ='bold';
    xtickangle(45)
    
if printFig
    cd(fullfile(projectPath,'4_Figures'))
    set(gcf,'PaperPositionMode','auto')
    if plotEEG; print(['Fig6 - Main Effect of the Condition'],'-djpeg','-r800')
    else; print(['Fig6 - Main Effect of the Condition (without EEG zscores)'],'-djpeg','-r800'); end; 
    cd(fullfile(projectPath,'3_analysis'))
end


% Condition order 
% ---------------

clear data2plot xPos boxCond boxZscore

% Figure preparation
figure('Units','centimeters','Position',[1 1 20 20],'Color',[1 1 1])     
colors = {[154, 3, 30]./256, [240, 5, 48]./256; ...
          [15, 76, 92]./256, [26, 131, 158]./256};
i=1; j=1; iBox = 1; ixPos = 1;

for iCond = 1:length(condNames)
    for iGroup = 1:length(groups)
        group       = groups{iGroup};
        condName    = condNames{iCond};
        if exist('EEG.FFT')
            subjects    = EEG.FFT.(group).subjects;
        elseif ~exist('EEG.FFT') && strcmp(group,'Musicians')
            subjects    = fieldnames(AnalyzedEEG_Musicians_blsnrFFT.(condName));
        elseif ~exist('EEG.FFT') && strcmp(group,'NonMusicians')
            subjects    = fieldnames(AnalyzedEEG_NonMusicians_blsnrFFT.(condName));
        end
        idx2keep        = find(contains(subjects,'sub'));
        subjects        = subjects(idx2keep); 
        
        for iSubject = 1:length(subjects)

            % find the subjectIdx corresponding to each of the condition order
            if strcmp(Log.(group).(subjects{iSubject}).condOrder{1},'ABAB1') 

                data2plot(i,1) = Zscores.(methodNames{zMethod}).(group).(condName).zMetRel(iSubject);
                i = i+1;

            elseif strcmp(Log.(group).(subjects{iSubject}).condOrder{1},'CDEF')   

                data2plot(j,2) = Zscores.(methodNames{zMethod}).(group).(condName).zMetRel(iSubject);
                j = j+1;
            end 
        end
    end
end

for iCondOrder = 1:2
    scatter(repmat(iCondOrder,length(data2plot),1),data2plot(:,iCondOrder), ...
            'SizeData',100,...
            'markerEdgeColor',colors{iCondOrder,1}, 'markerEdgeAlpha',0.6, ...
            'markerFaceColor',colors{iCondOrder,2}, 'markerFaceAlpha',0.4, ...
            'jitter','on', 'jitterAmount',0.07); hold on

    % boxplot preperation
    boxZscore(iBox:iBox+length(data2plot)-1) = data2plot(:,iCondOrder);
    boxCond(iBox:iBox+length(data2plot)-1)   = repmat(iCondOrder,length(data2plot),1);
    xPos(ixPos) = iCondOrder;

    iBox  = iBox+length(data2plot);
    ixPos = ixPos+1;
        
end
    
h = boxplot(boxZscore,boxCond,'Colors','k','Widths',0.15,'Positions',xPos,'symbol','','Labels','','Whisker',0, 'symbol', '');    % ,'Notch','on'
    set(h,{'linew'},{linewidth-0.25});

% layout
set(gca,'XTick',[1:2],'TickDir','out','LineWidth',linewidth,'fontsize',fontsize,'XLim',[0.6 2.4])%,'YLim',[-1.3 2.1],'YTick',[-1:0.5:2]) %,'YTick',[-0.5:0.25:0.5],'YLim',[-0.65 0.3] 
box off
ylabel('Beat-Related Zscores','fontsize',fontsize)
line([0 10],[0 0],'LineStyle','--','Color','k','LineWidth',0.5)
    
% xlabel
    % https://nl.mathworks.com/matlabcentral/answers/101922-how-do-i-create-a-multi-line-tick-label-for-a-figure-using-matlab-7-10-r2010a
    xlabel_row1 = {'Maximum','Minimum'};
    xlabel_row2 = {'Prior Context','Prior Context'};
    labelArray  = [xlabel_row1; xlabel_row2];
    labelArray  = strjust(pad(labelArray),'left'); % 'left'(default)|'right'|'center
    tickLabels  = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    ax          = gca();
    ax.TickLabelInterpreter = 'tex';  % needed for some plots like boxplot.
    ax.XTick        = 1:2; 
    ax.XLim         = [0.6 2.4];
    ax.XTickLabel   = tickLabels; 
    ax.FontSize     = fontsize;
    ax.FontWeight   ='bold';
    xtickangle(45)
    
    
if printFig
    cd(fullfile(projectPath,'4_Figures'))
    set(gcf,'PaperPositionMode','auto')
    print(['Fig6 - Main Effect of the Condition Order'],'-djpeg','-r800')
    cd(fullfile(projectPath,'3_analysis'))
end