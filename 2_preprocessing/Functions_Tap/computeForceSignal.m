function [lwdata] = computeForceSignal(lwdata)
%% 
%
%
%%
    fprintf ('  -> Computing force signal \n')
      
%%%%%%%%%%%%%%%%%%%%%%%   
% Force signal computing %
%%%%%%%%%%%%%%%%%%%%%%%


    % Indicate calibration coefficients
    calib_coeff_fz = [-10.7662351676953 -0.14181979074703 ...
    -10.2649711848043 0.17631695409008 -10.4154237298717 ...
    0.4954806017287];

    forces = zeros(lwdata.header.datasize(1),lwdata.header.datasize(6));
    for iTrial = 1:lwdata.header.datasize(1)
        trialComponents = squeeze(lwdata.data(iTrial,1:6,1,1,1,:));
        for iComponent = 1:size(trialComponents,1)
            mean_baseline = mean(trialComponents(iComponent, 1:(5/lwdata.header.xstep)));
            trialComponents(iComponent,:) = trialComponents(iComponent,:) - mean_baseline;
            trialComponents(iComponent, :) = calib_coeff_fz(iComponent) .* ...
                trialComponents(iComponent, :);
        end
        forces(iTrial,:) = -sum(trialComponents);
        %forces(iTrial,:) = forces(iTrial)/1000000; if you want you force
        %in Volts 
        lwdata.data(iTrial,2,1,1,1,:)=forces(iTrial,:);
        lwdata.data(iTrial,1,1,1,1,:)=lwdata.data(iTrial,8,1,1,1,:);
    end

    chanel_name = {lwdata.header.chanlocs.labels};
    chanel_name = chanel_name(1:2);
    option = struct('type','channel','items',{chanel_name},'suffix','sel_chan','is_save',1);
    lwdata = FLW_selection.get_lwdata(lwdata,option);
    option=struct('old_channel',{chanel_name},'new_channel',{{'Sound1', 'Force1'}},'suffix','chanlabels','is_save',1);
    lwdata = FLW_electrode_labels.get_lwdata(lwdata,option);
   
   fprintf ('  ==> Force signal computed \n')
    
end