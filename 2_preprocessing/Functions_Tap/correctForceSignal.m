function lwdataset = correctForceSignal(lwdataset, group)

fprintf ('  -> Switching force signal with sound if necessary \n')

%%%%%%%%%%%%%%%%%%%%%%%
% Force signal correction %
%%%%%%%%%%%%%%%%%%%%%%%

iTrial = 1;
lwdataTemp = lwdataset(1);
signals = squeeze(lwdataTemp.data(iTrial,:,1,1,1,:));
names = strsplit(lwdataTemp.header.name, " ");
subjectName = string(names(end));
chanels = {lwdataTemp.header.chanlocs.labels};
figure('Color',[1 1 1 1])
ForceIdx = 1;
SoundIdx = 1;
NoiseIdx = 1;
newChan = cell(size(chanels));
for iSign = 1:size(chanels,2)
    for iPlot = 1:size(chanels,2)
        if(size(chanels,2)>1)
            subplot(size(signals,1), 1, iPlot);
            if (iPlot == iSign)
                plot(signals(iPlot, :));
            else
                plot(signals(iPlot, :), "Color", [0.7 0.7 0.7]);
            end
        else
            plot(signals);
        end
        hold on;
    end
    signalType = questdlg('Which type of signal is it ?','Taps preprocessing','Force','Sound','Noise','Noise');
    switch signalType
        case 'Force'
            name = signalType + string(ForceIdx);
            newChan{iSign} = name;
            ForceIdx = ForceIdx + 1;
        case 'Sound'
            name = signalType + string(SoundIdx);
            newChan{iSign} = name;
            SoundIdx = SoundIdx + 1;
        case 'Noise'
            name = signalType + string(NoiseIdx);
            newChan{iSign} = name;
            NoiseIdx = NoiseIdx + 1;
    end
end

option=struct('old_channel',{chanels},'new_channel',{newChan},'suffix','chanlabels','is_save',1);

for iData = 1:length(lwdataset)
    lwdataset(iData) = FLW_electrode_labels.get_lwdata(lwdataset(iData),option);
end

close All;
fprintf ('  ==> Force signal corrected \n')

end