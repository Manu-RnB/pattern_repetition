function [zScoreOut] = Zscore(sFFT,freqVec,varargin)

%%%%%%%%%%%%%%%%%%
% Initial setups %
%%%%%%%%%%%%%%%%%%

global cfg Cfg EEG_zscores

if nargin == 3
    
    Options     = varargin{1};
    optionNames = fieldnames(Options);

    if ~isempty(find(strcmp(optionNames,'gridIOI')))
        cfg.gridIOI = Options.gridIOI;
    end   
    
    if ~isempty(find(strcmp(optionNames,'tripleMeter')))
        tripleMeter = Options.tripleMeter;
    else
        tripleMeter = 'No';
    end
    
    if ~isempty(find(strcmp(optionNames,'iCond')))
        condNames   = Cfg.condNames;
        iCond       = Options.iCond;
    end
    
    if ~isempty(find(strcmp(optionNames,'typeOfSignal')))
        typeOfSignal = Options.typeOfSignal;
    end
    
    if ~isempty(find(strcmp(optionNames,'NonRepProject')))
        NonRep      = true;
    else
        NonRep = false;
    end

end

iZMethod    = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set zscore methods - Classical Methods %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Zscores 
    frex              = 1/(cfg.gridIOI*12) * [1:12*20]; % *20 allowes to select higher harmonics (up to 100Hz)

    % Method 1: Classical method
    methodName{iZMethod,:}   = {'1,3,6,9,12 / 2,4,5,7,8,10,11'};
    metRel{iZMethod,:}       = [1,3,6,9,12];
    metUnrel{iZMethod,:}     = [2,4,5,7,8,10,11];
    iZMethod                 = iZMethod +1;
    
    % Method 2: Classical method minus 1
    methodName{iZMethod,:}   = {'3,6,9,12 / 2,4,5,7,8,10,11'};
    metRel{iZMethod,:}       = [3,6,9,12];
    metUnrel{iZMethod,:}     = [2,4,5,7,8,10,11];
    iZMethod                 = iZMethod +1;
    
    % Method 3: Classical method minus 12
    methodName{iZMethod,:}   = {'1,3,6,9 / 2,4,5,7,8,10,11'};
    metRel{iZMethod,:}       = [1,3,6,9];
    metUnrel{iZMethod,:}     = [2,4,5,7,8,10,11];  
    iZMethod                 = iZMethod +1;
    
    % Method 4: Classical method minus 1 and 12
    methodName{iZMethod,:}   = {'3,6,9 / 2,4,5,7,8,10,11'};
    metRel{iZMethod,:}       = [3,6,9];
    metUnrel{iZMethod,:}     = [2,4,5,7,8,10,11]; 
    iZMethod                 = iZMethod +1;
    
    % Method 5: Classical method minus 1, 6 and 12 (for the triple to
    % dupple comparison)
    methodName{iZMethod,:}   = {'3,9 / 2,4,5,6,7,8,10,11'};
    metRel{iZMethod,:}       = [3,9];
    metUnrel{iZMethod,:}     = [2,4,5,6,7,8,10,11]; 
    iZMethod                 = iZMethod +1;
    
    % Method 6: Classical method minus 1, 2 and 12
    methodName{iZMethod,:}   = {'3,6,9 / 4,5,7,8,10,11'};
    metRel{iZMethod,:}       = [3,6,9];
    metUnrel{iZMethod,:}     = [4,5,7,8,10,11]; 
    iZMethod                 = iZMethod +1;
    
    
    if strcmp(tripleMeter,'Yes')
        
        % Method 7: Classic Triple Meter
        methodName{iZMethod,:}   = {'1,2,4,8,10,12 / 3,5,6,7,9,11'};
        metRel{iZMethod,:}       = [1,2,4,8,10,12];
        metUnrel{iZMethod,:}     = [3,5,6,7,9,11];
        iZMethod                 = iZMethod +1;
        
        % Method 7: Classic Triple Meter minus 1
        methodName{iZMethod,:}   = {'2,4,8,10,12 / 3,5,6,7,9,11'};
        metRel{iZMethod,:}       = [2,4,8,10,12];
        metUnrel{iZMethod,:}     = [3,5,6,7,9,11];
        iZMethod                 = iZMethod +1;
        
        % Method 8: Classic Triple Meter minus 12
        methodName{iZMethod,:}   = {'1,2,4,8,10 / 3,5,6,7,9,11'};
        metRel{iZMethod,:}       = [1,2,4,8,10];
        metUnrel{iZMethod,:}     = [3,5,6,7,9,11];
        iZMethod                 = iZMethod +1;
        
        % Method 9: Classic Triple Meter minus 1 and 12
        methodName{iZMethod,:}   = {'2,4,8,10 / 3,5,6,7,9,11'};
        metRel{iZMethod,:}       = [2,4,8,10];
        metUnrel{iZMethod,:}     = [3,5,6,7,9,11];
        iZMethod                 = iZMethod +1;
    end
    
    
    % selectedIdx & selectedFrex
    for i = 1:size(methodName,1)
        
        selectedIdx{i,:}  = sort([metRel{i,:}, metUnrel{i,:}]);
        selectedFrex{i,:} = sort([frex(metRel{i,:}), frex(metUnrel{i,:})]);
    end
    
    cfg.Zscore.methods = table(methodName,selectedFrex,selectedIdx,metRel,metUnrel);
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Zscore - Classical Methods % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for iZMethod = 1:size(methodName,1)
        
        % compute zscore
        frexIdx                     = dsearchn(freqVec',selectedFrex{iZMethod,:}');
        Amps                        = NaN (1,selectedIdx{iZMethod}(end));
        Amps(selectedIdx{iZMethod}) = sFFT(frexIdx);
        zscores                     = (Amps - nanmean(Amps)) ./ nanstd(Amps);
        zscoresMetRel               = mean(zscores(metRel{iZMethod}));
        zscoresMetUnrel             = mean(zscores(metUnrel{iZMethod}));
        
        % save zScore structure
        zScoreOut(iZMethod).methodName      = methodName{iZMethod};
        zScoreOut(iZMethod).freq            = freqVec(frexIdx);
        zScoreOut(iZMethod).frexIdx         = frexIdx';
        zScoreOut(iZMethod).amps            = Amps;
        zScoreOut(iZMethod).metRelIdx       = metRel{iZMethod,:};
        zScoreOut(iZMethod).metUnrelIdx     = metUnrel{iZMethod,:};
        zScoreOut(iZMethod).zMetRel         = zscoresMetRel;
        zScoreOut(iZMethod).zMetUnrel       = zscoresMetUnrel; 
        
        % add NaNs in the freq and frexIdx of the first frequencies are not taken into account (helps with plotting later) 
        if ~ismember(1,metRel{iZMethod,:})
            zScoreOut(iZMethod).frexIdx     = [NaN, zScoreOut(iZMethod).frexIdx];
            zScoreOut(iZMethod).freq        = [NaN, zScoreOut(iZMethod).freq];
        end
        
        if ~ismember(2,metUnrel{iZMethod,:})
            zScoreOut(iZMethod).frexIdx     = [NaN, zScoreOut(iZMethod).frexIdx];
            zScoreOut(iZMethod).freq        = [NaN, zScoreOut(iZMethod).freq];
        end

        clear frexIdx Amps zscores zscoresMetRel zscoresMetUnrel
        
    end
  
    iZMethod = iZMethod+1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set zscore methods - Based on Tomas' sweeps experiment %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    % Method 10: Groupings in 12 and 16 events (Sweeps Experiment by Tomas)
    
    methodName{iZMethod,:}      = {'3,6,9 / Groupings in 12 and 16'};
    frex12                      = 1/(cfg.gridIOI*12) * [1:12];
    frex16                      = 1/(cfg.gridIOI*16) * [1:16];
    frex12_16                   = sort(unique([frex12,frex16]));
    metRelIdx{iZMethod,:}       = dsearchn(frex12_16',frex([3,6,9,12])')'; % The 9th was in the metRel frequencies ???
    metUnrelIdx{iZMethod,:}     = setdiff (1:length(frex12_16),metRelIdx{iZMethod,:}); 
    metRel{iZMethod,:}          = frex12_16(metRelIdx{iZMethod,:});
    metUnrel{iZMethod,:}        = frex12_16(metUnrelIdx{iZMethod,:});
    
    % compute zscore
    frexIdx                     = dsearchn(freqVec',frex12_16');
    Amps                        = NaN (1,length(frex12_16));
    Amps                        = sFFT(frexIdx);
    zscores                     = (Amps - nanmean(Amps)) ./ nanstd(Amps);
    zscoresMetRel               = mean(zscores(metRelIdx{iZMethod,:}));
    zscoresMetUnrel             = mean(zscores(metUnrelIdx{iZMethod,:}));

    % save zScore structure
    zScoreOut(iZMethod).methodName  = methodName{iZMethod};
    zScoreOut(iZMethod).freq        = freqVec(frexIdx);
    zScoreOut(iZMethod).frexIdx     = frexIdx;
    zScoreOut(iZMethod).amps        = Amps;
    zScoreOut(iZMethod).metRelIdx   = metRelIdx{iZMethod,:} ;
    zScoreOut(iZMethod).metUnrelIdx = metUnrelIdx{iZMethod,:};
    zScoreOut(iZMethod).zMetRel     = zscoresMetRel;
    zScoreOut(iZMethod).zMetUnrel   = zscoresMetUnrel; 

    clear frexIdx Amps zscores zscoresMetRel zscoresMetUnrel metRelIdx metUnrelIdx
    
    iZMethod = iZMethod+1;
    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % zscore methods - takes all possible metUnrel frequencies %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This technique is to be used with the NonRep project. We consider all
% frequencies in the input as frequencies of interest => the twoBars
% condition has 24 frequencies, the FourBars has 48, and the NonRep has
% 336.
    if NonRep

        methodName  = [{'All possible metUnrel frequencies depending on condition (until 30Hz)'};...
                       {'All possible metUnrel frequencies depending on condition (until 5Hz)'}];       
        iDur        = [6,1];
        freq2remove = 5:5:30;

        if strcmp(typeOfSignal,'EEG');      nFreq       = [48,24,48,336];
        elseif strcmp(typeOfSignal,'Stim'); nFreq       = [24,48,336];
        end

        for iMethod = 1:2

            freq = 1/(cfg.gridIOI*nFreq(iCond)) * [1:nFreq(iCond)*iDur(iMethod)]; 

            % remove frequencies lower than 1Hz + 5Hz and harmonics
            freq = freq(freq>1);
            freq (dsearchn(freq',freq2remove')) = [];

            frexIdx                     = dsearchn(freqVec',freq');
            metRelFrex                  = [1.25 : 1.25 : 30];
            metRelFrex (dsearchn(metRelFrex',freq2remove')) = [];
            metRelIdx                   = unique(dsearchn(freq',metRelFrex')); 
            if iMethod==2 % it selects the last metUnrel as metRel
                metRelIdx(end) = []; 
            end                            
            metUnrelIdx                 = setdiff([1:length(freq)],metRelIdx);

            amps                        = sFFT(frexIdx);
            zscores                     = (amps - nanmean(amps)) ./ nanstd(amps);
            zscoresMetRel               = mean(zscores(metRelIdx));
            zscoresMetUnrel             = mean(zscores(metUnrelIdx));   

            zScoreOut(iZMethod).methodName  = methodName(iMethod);
            zScoreOut(iZMethod).freq        = freq;
            zScoreOut(iZMethod).frexIdx     = frexIdx;
            zScoreOut(iZMethod).amps        = amps;
            zScoreOut(iZMethod).metRelIdx   = metRelIdx;
            zScoreOut(iZMethod).metUnrelIdx = metUnrelIdx;
            zScoreOut(iZMethod).zMetRel     = zscoresMetRel;
            zScoreOut(iZMethod).zMetUnrel   = zscoresMetUnrel; 

            iZMethod = iZMethod+1;
        end 
    end
    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % classic zscore method with 24 frequencies of interest %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This technique is to be used with the NonRep project. 
% Instead of looking at the 12 FOI, we take the 24.

    if NonRep
        
        % Method 7: Classical with 24 frequencies without 
        frex                       = 1/(cfg.gridIOI*24) * [1:24]; % *20 allowes to select higher harmonics (up to 100Hz)
        methodName{iZMethod,:}     = {'Classic with 24 FOI (minus < 0.5Hz and 5Hz)'};
        metRel{iZMethod,:}         = [3,6,12,18];%[3,6,9];
        metUnrel{iZMethod,:}       = [4,5,7,8,9,10,11,13,14,15,16,17,19,20,21,22,23]; % don't the frequencies below 1Hz in account
         
        
        selectedIdx{iZMethod,:}     = sort([metRel{iZMethod,:}, metUnrel{iZMethod,:}]);
        selectedFrex{iZMethod,:}    = sort([frex(metRel{iZMethod,:}), frex(metUnrel{iZMethod,:})]);
        
        % compute zscore
        frexIdx                     = dsearchn(freqVec',selectedFrex{iZMethod,:}');
        Amps                        = NaN (1,selectedIdx{iZMethod}(end));
        Amps(selectedIdx{iZMethod}) = sFFT(frexIdx);
        zscores                     = (Amps - nanmean(Amps)) ./ nanstd(Amps);
        zscoresMetRel               = mean(zscores(metRel{iZMethod}));
        zscoresMetUnrel             = mean(zscores(metUnrel{iZMethod}));
        
        % saving
        zScoreOut(iZMethod).methodName      = methodName{iZMethod};
        zScoreOut(iZMethod).freq            = freqVec(frexIdx);
        zScoreOut(iZMethod).frexIdx         = frexIdx';
        zScoreOut(iZMethod).amps            = Amps;
        zScoreOut(iZMethod).metRelIdx       = metRel{iZMethod,:};
        zScoreOut(iZMethod).metUnrelIdx     = metUnrel{iZMethod,:};
        zScoreOut(iZMethod).zMetRel         = zscoresMetRel;
        zScoreOut(iZMethod).zMetUnrel       = zscoresMetUnrel;
        
        iZMethod                    = iZMethod +1; 
        
    end 
end

