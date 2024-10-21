function icsRemoval(subName, group)

    global Cfg Log Paths
    
    fprintf ('  -> ICA removal \n')

    condNames   = Cfg.condNames;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Leave some time to look at the ICs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    letswave7
    pause(10)
    waitForMsgBoxCompletion = msgbox('Press Ok when you have had a look at the ICA'); pause(1)
    waitfor (waitForMsgBoxCompletion)  
    
    ICAQuest = questdlg('Would you like to apply an ICA?');
    
    switch ICAQuest
        case 'Yes'
            waitForMsgBoxCompletion = msgbox('Press Ok when you have aplpied the ICA'); pause(1)
            waitfor (waitForMsgBoxCompletion)  
            pause(15)
        case 'No'
    end        
             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Save the directory of the IC filtered files % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     cd(fullfile(Paths.LW,group,'Preprocessing/EEG/ICA')); 
%     
%     % Wait for 10 seconds 
%     pause(10)
%     lwDir = dir (['sp_filter*',subName,'*.lw6']);
%      
%     if size(lwDir,1) < length(condNames)                % If all files haven't been created yet
%         pause(10)
%         lwDir = dir (['sp_filter*',subName,'*.lw6']);
%         
%     elseif size(lwDir,1) > length(condNames)            % If too many versions of the analysis
%         
%         
%         MsgBox ('Too many sp_filter files were found. Load them manually')
%         pause (5)
%         
%         for iCond = 1:size(condNames,2)
%             [Cfg.Preprocessing.ICA.(condNames{iCond}).Dir.name, Cfg.Preprocessing.ICA.(condNames{iCond}).Dir.folder] ...
%                 = uigetfile(['sp_filter*',condNames{iCond},'*.lw6'],['pick the latest ',condNames{iCond},' sp_filter files']);
%         end
%         
%     else                                                % Normal condition    
%         for iCond = 1:size(condNames,2)
%             Cfg.Preprocessing.ICA.(condNames{iCond}).Dir = dir (['sp_filter*',condNames{iCond},'*.lw6']);
%         end        
%     end

    close all force % closes Letswave7
    
    fprintf ('  ==> ICs removed \n')
    
end

