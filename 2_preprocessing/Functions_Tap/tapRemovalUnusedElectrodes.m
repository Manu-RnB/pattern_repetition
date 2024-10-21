function [lwdata, Log] = tapRemovalUnusedElectrodes(lwdata, group, subName, ~, Log, Paths)
    
    fprintf ('  -> Removal of unused electrodes \n')
    if(~exist(fullfile(Paths.LW,group,'Preprocessing/Tapping/Removal_Elec')))
        mkdir(fullfile(Paths.LW,group,'Preprocessing/Tapping/Removal_Elec'));
    end
    cd(fullfile(Paths.LW,group,'Preprocessing/Tapping/Removal_Elec'));
      
%%%%%%%%%%%%%%%%%%%%%%%   
%  Electrode removal  %
%%%%%%%%%%%%%%%%%%%%%%%
   chanel_Ana = contains({lwdata.header.chanlocs.labels}, 'Ana');
   chanel_name = {lwdata.header.chanlocs.labels};
   new_chanel_name = chanel_name(chanel_Ana);
   option = struct('type','channel','items',{new_chanel_name},'suffix','sel_chan','is_save',1);
   lwdata = FLW_selection.get_lwdata(lwdata,option);
   Log.(group).(subName).TappingPreprocessing.RemovalElec.originalElec = chanel_name;
   Log.(group).(subName).TappingPreprocessing.RemovalElec.remainingElec = new_chanel_name;
   fprintf ('  ==> Electrodes removed \n')
    
end