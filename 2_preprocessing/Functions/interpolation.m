function [lwdata] = interpolation(lwdata, subName, group)

    global Log
    
    fprintf ('  -> Interpolation \n')

%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Channels interpolation %
%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Look for previous interpolation
    if isfield(Log.(group).(subName).Preprocessing,'Interpolation')
        
        interpElec  = Log.(group).(subName).Preprocessing.Interpolation.InterpolatedElectrode;
        closestElec = Log.(group).(subName).Preprocessing.Interpolation.ClosestElectrodes;
        
        fprintf ('  -> Previous interpolation \n    Interpolated electrode(s) : %c \n    Closest electrode(s) : %c', interpElec, closestElec )
    end
    

    InterpQuest1 = questdlg ('Would you like to interpolate channels?', ...
                             'Channel Interpolation','Yes','No','No');
                         
                         
    switch InterpQuest1
        case 'Yes'
            
            waitForMsgBoxCompletion = msgbox('Press Ok when ready to start the interpolation');  
            letswave7 
            waitfor (waitForMsgBoxCompletion)
            pause(30)
            
            interpFinished = 'No';
            
            while strcmp(interpFinished,'No') == 1
                
                InterpQuest2 = questdlg('Have you finished the channel interpolation ?', ...
                                        'Channel Interpolation','Yes','No','No');
                    switch InterpQuest2
                        case 'Yes'
                            
                            % Find the latest file in the LW directory
                            lwDir = dir (['chan_interp*',subName,'*.lw6']);

                            for i = 1: size(lwDir,1)
                                time(i,1) = lwDir(i).datenum;            
                            end

                            [~,timeIdx] = max(time);
                            
                            % if the latest file begins with 'chan_interp', load it
                            if strncmp (string (lwDir(timeIdx).name),'chan_interp', 11) == 1
                                option = struct('filename',fullfile(lwDir(timeIdx).name));
                                lwdata = FLW_load.get_lwdata(option);
                                
                                % Complete a dialogbox to indicate the parameters of the interpolation (further used in the ReadMe.txt)
                                chan_interp = inputdlg({'Enter the channel(s) that was/were interpolated', ...
                                                        'Enter the number of closest electrodes selected'},'Interpolation');

                                Log.(group).(subName).Preprocessing.Interpolation.InterpolatedElectrode = chan_interp {1};
                                Log.(group).(subName).Preprocessing.Interpolation.ClosestElectrodes     = chan_interp {2};
                                
                                fprintf ('  ==> Channel(s) correctly interpolated \n')

                            else
                                warning ('The file that was created last is not a chan_interp file')                           
                                Log.(group).(subName).Preprocessing.Interpolation.InterpolatedElectrode = 'None';
                                Log.(group).(subName).Preprocessing.Interpolation.ClosestElectrodes     = 'None';
                                
                                fprintf ('  ==> No channel interpolated \n')
                            end
                            
                            interpFinished = 'Yes'; % stops the while loop
                          
                            
                        case 'No'
                            pause(30)
                    end
            end         

        case 'No'
            fprintf ('  ==> No channel interpolated \n')
            Log.(group).(subName).Preprocessing.Interpolation.InterpolatedElectrode = 'None';
            Log.(group).(subName).Preprocessing.Interpolation.ClosestElectrodes     = 'None';
    end
    
    close all force
end

