function [lwdataset, Log] = tapEpochSegmentation(lwdata, group, subName, Cfg, Log, Paths)

fprintf ('  -> Epoch Segmentation and correction\n')

if(~exist(fullfile(Paths.LW,group,'Preprocessing/Tapping/Segmentation')))
    mkdir(fullfile(Paths.LW,group,'Preprocessing/Tapping/Segmentation'));
end
cd(fullfile(Paths.LW,group,'Preprocessing/Tapping/Segmentation'));

condNames       = strcat(Cfg.condNames,'_tap');

if isfield(Cfg,'trialDur')  
    trialDur        = Cfg.trialDur;
else
    trialDur        = min(Log.(group).(subName).Design.duration);
end

%%%%%%%%%%%%%%%%%%%%%%%
% Epochs Segmentation %
%%%%%%%%%%%%%%%%%%%%%%%
chanels = {lwdata.header.chanlocs.labels};
Log.(group).(subName).TappingPreprocessing.EpochSegmentation.OldElec = chanels;

if(length(chanels)==8) % New tapping box
    %Segment with 5 seconds before for baseline (force computing)
    option      = struct('event_labels',{condNames}, ...
        'x_start',-5, ...
        'x_end',trialDur, ...
        'x_duration',trialDur+5, ...
        'suffix','ep','is_save',1);
    lwdataset  = FLW_segmentation_separate.get_lwdataset(lwdata,option);

    % Compute force signal 

    for fieldIdx = 1:length(lwdataset)
        lwdataset(fieldIdx) = computeForceSignal(lwdataset(fieldIdx));
        %cropping the 5 seconds 
        option=struct('xcrop_chk',1,'xstart',0,'xend',trialDur,'suffix','crop','is_save',1);
        lwdataset(fieldIdx) = FLW_crop_epochs.get_lwdata(lwdataset(fieldIdx),option);
    end
else %Old tapping box, just matching signal with their names
    option      = struct('event_labels',{condNames}, ...
        'x_start',0, ...
        'x_end',trialDur, ...
        'x_duration',trialDur, ...
        'suffix','ep','is_save',1);
    lwdataset  = FLW_segmentation_separate.get_lwdataset(lwdata,option);

    lwdataset = correctForceSignal(lwdataset, group);
end

chanels = {lwdataset(1).header.chanlocs.labels};
Log.(group).(subName).TappingPreprocessing.EpochSegmentation.NewElec = chanels;

fprintf ('  ==> Epochs segmented and corrected\n')

end

