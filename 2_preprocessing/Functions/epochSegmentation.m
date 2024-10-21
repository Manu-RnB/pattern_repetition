function [lwdata] = epochSegmentation(lwdata, subName, group)

    global Cfg Log Paths
    addpath (genpath (Paths.lw6));
    
    fprintf ('  -> Epoch Segmentation \n')
    
    cd(fullfile(Paths.LW,group,'Preprocessing/EEG/Segmentation'));
    
    condNames       = strcat(Cfg.condNames,'_listen');

    currentDesign   = Log.(group).(subName).Design;
    trialDur        = min(currentDesign.duration);
      
%%%%%%%%%%%%%%%%%%%%%%%   
% Epochs Segmentation %
%%%%%%%%%%%%%%%%%%%%%%%

    option      = struct('event_labels',{condNames}, ...
                         'x_start',0, ...
                         'x_end',trialDur, ...
                         'x_duration',trialDur, ...
                         'suffix','ep','is_save',1);
    % lw7
    % lwdataset   = FLW_segmentation_separate.get_lwdataset(lwdata,option);
    
    % lw6 : [out_datasets,message_string]=RLW_segmentation2(header,data,event_labels,varargin)
    % [lwdataset,message_string] = RLW_segmentation2(lwdata.header,lwdata.data,{condNames},'x_start',0,'x_duration',trialDur);
    for iCond = 1:length(condNames)
        [lwdataSegm.header, lwdataSegm.data] = RLW_segmentation(lwdata.header, lwdata.data, {condNames{iCond}}, 'x_start',0, 'x_duration', trialDur);
        lwdataSegm.header.name = ['ep_',condNames{iCond},' ',lwdataSegm.header.name];
        CLW_save (fullfile(Paths.LW,group,'Preprocessing/EEG/Segmentation'), ...
                  lwdataSegm.header, lwdataSegm.data);
    end
    
    fprintf ('  ==> Epochs segmented \n')
    
end

