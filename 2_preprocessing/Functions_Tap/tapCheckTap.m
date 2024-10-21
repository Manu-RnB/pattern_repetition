function [lwdataset, TappingTaps, TapAnomaly, Log] = tapCheckTap(lwdataset, group, subName, TappingTaps, TapAnomaly, Cfg, Log, Paths)

    checkTapBool = questdlg('Do you want to check the tapping data ?','Taps preprocessing','Yes','No','No');
    checkTapBool = strcmp(checkTapBool, 'Yes');
    if(checkTapBool)
        fprintf ('  -> Checking Taps \n')
        if ~exist(fullfile(Paths.LW,group,'Preprocessing/Tapping/CheckTaps'))
            mkdir(fullfile(Paths.LW,group,'Preprocessing/Tapping/CheckTaps'));
        end
        cd(fullfile(Paths.LW,group,'Preprocessing/Tapping/CheckTaps'));
    
        try
            checkSubExistence = isfield(TappingTaps.AAAA, subName);
        catch %if TappingTaps does not exist
            checkSubExistence = False;
        end
    
        if checkSubExistence
            overwrite         = questdlg('Do you want to overwrite or check?', 'Taps detection','Overwrite','Check', 'Check');
        else
            overwrite         = "Overwrite";
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%
        % Checking Taps %
        %%%%%%%%%%%%%%%%%%%%%%%
    
        gridIOI    = Cfg.gridIOI;
        frex       = 1/(gridIOI*24) * [1:24];
        frex       = frex([2,3,4,6,8,12,24]);
        trialNames = ["Trial1", "Trial2", "Trial3"];
    
        for lwdata = lwdataset
            condName       = extractBetween(lwdata.header.name,"ep_","_tap");
            datasize       = lwdata.header.datasize;
            chanelNames    = string({lwdata.header.chanlocs.labels});
            fs             = 1/lwdata.header.xstep;
            trialDur       = lwdata.header.xstep * datasize(6) - lwdata.header.xstep;
            time           = lwdata.header.xstart:lwdata.header.xstep:trialDur;

            iTrial = 1;
            while iTrial <= datasize(1)
                trialName = trialNames(iTrial);
                if strcmp(overwrite, "Overwrite")
                    tapsIdx        = cell2mat({lwdata.header.events.epoch})==iTrial;
                    tapsIdxtmp     = string({lwdata.header.events.code})=='4';
                    tapsIdx        = tapsIdx & tapsIdxtmp;
                    taps           = cell2mat({lwdata.header.events.latency});
                    taps           = taps(tapsIdx);
                else
                    taps           = TappingTaps.(char(condName)).(subName).(trialName);
                end
                ITI                        = diff(taps); %% Automatically sees double taps
                ITIwdIdx = ITI            >= 0.18; %% Treshold that works great, you can change it
                ITI_wd                     = ITI(ITIwdIdx); %% ITI without double taps
                aimingPeriod               = median(ITI_wd);
                tapswdIdx                  = [1, ITIwdIdx]; %% Indexes of taps that are not double taps, shifted by one because computed with ITI
                taps_wd                    = taps(logical(tapswdIdx)); %% Taps without double taps
                doubleTaps                 = setdiff(taps, taps_wd); %% double taps
                tapRemoved                 = [];
                tapAdded                   = [];
                if strcmp(overwrite, "Check") %% get double taps from structure if only check
                    try %% try because structure may not have be iniated for this subject
                        doubleTaps         = TapAnomaly.(char(condName)).(subName).(trialName).dbTaps;
                        try
                            tapRemoved     = TapAnomaly.(char(condName)).(subName).(trialName).rm;
                            tapAdded       = TapAnomaly.(char(condName)).(subName).(trialName).add;
                        catch
                        end
                    catch
                    end
                end
    
                [~,Idx]                    = min(abs(frex - 1./aimingPeriod));        % index of closest frex
                aimingFrex                 = frex(Idx);
                aimingPeriod               = 1./aimingFrex;
                if sum(contains(chanelNames, "Sound"))>0 %% check if there is sound
                    soundIdx               = contains(chanelNames, "Sound");
                    firstSoundIdx          = find(soundIdx, 1, 'first');
                    soundSignal            = lwdata.data(iTrial,firstSoundIdx,1,1,1,:);
                    soundEnvelope          = envelope(squeeze(soundSignal)');
                    soundEventIdx          = soundEnvelope>(((max(soundEnvelope)-min(soundEnvelope))/2)+min(soundEnvelope)); %% gets where there is sound or not, treshold may be adapted
                    % the idea here is to get the first occurence where the
                    % soundEventIdx goes from 0 to 1. So the first silence
                    % to sound. So we will substract the soundEvent with a
                    % shifted version of itself.
                    shiftRSoundEventIdx      = [soundEventIdx(1), soundEventIdx];
                    soundEventIdx_Extended   = [soundEventIdx, soundEventIdx(end)]; %% extended for the substraciton
                    subSoundEventIdx         = shiftRSoundEventIdx-soundEventIdx_Extended;
                    subSoundEventIdx         = subSoundEventIdx == -1; %% because find() needs a positive scalar
                    idxSecondSoundEvent      = find(subSoundEventIdx, 1, 'first'); %% it will be the index of the onset of the first or second sound event, in case the first one is cutoff.
                    onsetFundamentalGrid     = time(idxSecondSoundEvent);
                    if isempty(onsetFundamentalGrid) %security
                        onsetFundamentalGrid = 0;
                    end
                else
                    onsetFundamentalGrid     = 0;
                end
                fundamentalGrid              = onsetFundamentalGrid:gridIOI:trialDur;
                fundamentalGridIdx           = dsearchn(fundamentalGrid', taps_wd');
                onsetLastGrid                = fundamentalGrid(fundamentalGridIdx(end)); %% get the nearest fundamental grid from the last tap and synchronise the grid on it
                onsetFirstGrid               = (onsetLastGrid - floor(onsetLastGrid/aimingPeriod)*aimingPeriod); %% get the first from the last by dividing by aimingPeriod
                aimingGrid                   = onsetFirstGrid:aimingPeriod:onsetLastGrid; %% the grid the participant might have wanted to tap to
    
                %% plotting
    
                figure('units','normalized','outerposition',[0 0 1 1],'Color',[1 1 1])
    
                %% Pour Manu, pour moi ça ralentit trop le code. Sinon, il faut retirer le plot du grid. 
%                 for gridValue = fundamentalGrid
%                     xline(gridValue, 'color', [0.1, 0.1, 0.1], 'LineStyle','--');
%                     hold on
%                 end
                for gridValue = aimingGrid
                    xline(gridValue, 'color', [.3 .3 .3]);
                    hold on
                end
                %%

%                 xline(fundamentalGrid, 'color', [0.1, 0.1, 0.1], 'LineStyle','--');
%                 hold on
% 
%                 xline(aimingGrid, 'color', [.3 .3 .3]);
%                 hold on

                if sum(contains(chanelNames, "Sound"))>0 %% check if there is sound
                    plot(time,(1/4)*rescale(soundSignal(1,:)), 'Color', [.3 .3 .3]); %% plot sound
                    hold on
                end
    
                if sum(contains(chanelNames, "Force"))>0 %% check if there is a force signal
                    forceIdx               = contains(chanelNames, "Force");
                    firstForceIdx          = find(forceIdx, 1, 'first');
                    forceSignal            = lwdata.data(iTrial,firstForceIdx,1,1,1,:);
                    h1                     = plot(time, rescale(forceSignal(1,:)), 'Color', 'blue');
                    hold on
                end

                if ~(isempty(doubleTaps) && isempty(tapRemoved)) %si aucun des deux n'est empty et qu'on overwrite
                    h2 = stem(sort([doubleTaps,tapRemoved]), ones(length(doubleTaps)+length(tapRemoved)), 'color', [1, 0.5, 0.8]);
                    hold on
                end
    
                h3 = stem(taps_wd,ones(length(taps_wd)),'r');
                hold on

                if ~isempty(tapAdded) % si il y a des taps rajoutés et qu'on overwrite
                    h4 = stem(tapAdded, ones(length(tapAdded)), 'color', [0.6, 0, 0.2]);
                    hold on
                end

                set(gca,'TickDir','out','LineWidth',1.25,'fontsize',14)
                box off
    
                for iXLim = 1 : ceil(trialDur/10)
                    plotrange = [0+(iXLim-1)*10 10+(iXLim-1)*10];
                    xlim(plotrange);
                    title (strcat(subName, ' - ', condName,' - ',trialName),'Position',[5+(iXLim-1)*10 1.03 0])
    
                    % Manually add missing taps
                    wrongTap  = questdlg('Are there taps to add or remove?','Taps detection','Add','Remove','No','No');
    
                    while ~strcmp(wrongTap,'No')
                        % select regions
                        selectedRange    = ginput(2);
                        selectedRangeIdx = [round(selectedRange(1,1)*fs) : round(selectedRange(2,1)*fs)];
    
                        % find if a tap is comprised in the selected region
                        pos = dsearchn(time(selectedRangeIdx)',taps_wd');
                        pos (pos == 1) = [];
                        pos (pos == length(selectedRangeIdx)) = [];
    
                        % remove the tap that was comprised in the selected region
                        if ~isempty(pos)
                            tapIdx2remove                     = dsearchn(taps_wd', time(selectedRangeIdx(pos))'); %tout les taps à retirer
                            [~, ~, ifAddedIdx]                = intersect(taps_wd(tapIdx2remove), tapAdded); % les indices des taps qui sont des rajouts
                            [~, ifRemovedIdx, ~]              = setxor(taps_wd(tapIdx2remove),tapAdded);     % les indices des taps qui ne sont pas des rajouts
                            if ~isempty(ifRemovedIdx)
                                tapRemoved                    = [tapRemoved, taps_wd(ifRemovedIdx)];                       % les taps qui ne sont pas des rajouts sont mis dans les enlevés
                            end
                            if ~isempty(ifAddedIdx)
                                tapAdded(ifAddedIdx)          = [];                                          % les rajouts retirés sont retirés de la liste des rajoutés
                            end
                            taps_wd(tapIdx2remove)            = [];                                          % on enleve les taps
    
                        else % add a new tap
                            if sum(contains(chanelNames, "Force"))>0 %% check if there is Force
                                [~,tapIdx2add]   = max(forceSignal(selectedRangeIdx)); % index within the selected range
                            else
                                anySignal = lwdata.data(iTrial,1,1,1,1,:);
                                [~,tapIdx2add]   = max(anySignal(selectedRangeIdx)); % index within the selected range
                            end
                            tapIdx2add           = tapIdx2add + round(selectedRange(1,1)*fs);
                            taps_wd(end+1)       = time(tapIdx2add);
                            tapAdded             = [tapAdded, time(tapIdx2add)'];
                        end
    
                        taps_wd = sort(taps_wd);
    
                        % update figure
    
                        try
                            delete(h2);
                        catch
                        end
                        delete(h3);
                        try
                            delete(h4);
                        catch
                        end

                        if ~(isempty(doubleTaps) && isempty(tapRemoved)) %%si aucun des deux n'est empty
                            h2 = stem(sort([doubleTaps,tapRemoved]), ones(length(doubleTaps)+length(tapRemoved)), 'color', [1, 0.5, 0.8]);
                            hold on
                        end

                        h3 = stem(taps_wd,ones(length(taps_wd)),'r');
                        hold on
                        
                        if ~isempty(tapAdded) %% si il y a des taps rajoutés
                            h4 = stem(tapAdded, ones(length(tapAdded)), 'color', [0.6, 0, 0.2]);
                            hold on
                        end
    
    
                        % other missing taps?
                        wrongTap       = questdlg('Are there taps to add or remove?','Taps detection','Add','Remove','No','No');
                    end
                end
                wrongTap       = questdlg('Do you want to redo this trial ?','Taps detection','Yes','No','Yes');
                if strcmp(wrongTap,'Yes')
                    close all
                else
                    close all
                    TapAnomaly.(char(condName)).(subName).(trialName).rm     = sort(tapRemoved);
                    TapAnomaly.(char(condName)).(subName).(trialName).add    = sort(tapAdded);
                    TapAnomaly.(char(condName)).(subName).(trialName).dbTaps = sort(doubleTaps);
                    TappingTaps.(char(condName)).(subName).(trialName)       = sort(taps_wd);
                    save('TappingTaps.mat','TappingTaps','TapAnomaly','Paths','Cfg','Log', '-v7.3')
                    iTrial = iTrial+1;
                end
            end
        end
    end
end