function [lwdata] = electrodes(lwdata, subName, group)

    global Cfg Log Paths
    
    fprintf ('  -> Removal of unused electrodes \n')
    
    elecLabels  = Cfg.Preprocessing.elecLabels.labels;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Removal of unused channels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Such as accelerometer channels
    if isfolder(Paths.lw6)
        movefile(Paths.lw6,Paths.lw6OutsideMatlabPath)
    end

    cd(fullfile(Paths.LW,group,'Preprocessing/EEG/ElectrodeRemoval'));
    
    option = struct('type','channel','items',{table2cell(elecLabels)'}, ...
                    'suffix','sel_chan','is_save',1);
                
    lwdata = FLW_selection.get_lwdata(lwdata,option);
    
    
    fprintf ('  ==> Unused channels removed \n')
    
    movefile(Paths.lw6OutsideMatlabPath,Paths.lw6)
    
%% Add electrodes' location

    fprintf ('  -> Add electrode location \n')

    addpath (genpath (Paths.lw6)); 
    
    locs = readlocs(Paths.lw6ElecLoc);
    
    for iLocs = 1:length(locs)
        locs(iLocs).topo_enabled = 1;
    end
    
    [lwdata.header] = RLW_edit_electrodes(lwdata.header,locs);
    
    rmpath(Paths.lw6);
    
    fprintf ('  ==> electrode location added \n')
    
%% Save Log file

    Log.(group).(subName).Preprocessing.ElecLabels = elecLabels;

end

