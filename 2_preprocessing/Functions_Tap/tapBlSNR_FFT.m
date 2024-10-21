function [lwdatasetFFT, lwdatasetBlSNR, Log] = tapBlSNR_FFT(lwdatasetFFT, group, subName, Cfg, Log, Paths)

    fprintf ('  -> Baseline SNR correction \n')
    
    if ~exist(fullfile(Paths.LW,group,'Preprocessing/Tapping/BlSNR'))
        mkdir(fullfile(Paths.LW,group,'Preprocessing/Tapping/BlSNR'));
    end
    cd(fullfile(Paths.LW,group,'Preprocessing/Tapping/BlSNR'));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % BlSNR %
    %%%%%%%%%%%%%%%%%%%%%%%
    condFindTxt = 'ABAB1 ABAB2 AAAA CDEF';
    option = struct('xstart',Cfg.Preprocessing.Tapping.Bl_snr.lowBin,...
        'xend',Cfg.Preprocessing.Tapping.Bl_snr.highBin,...
        'num_extreme',0,'operation','subtract','suffix','bl_snr','is_save',1);
    for iData = 1:length(lwdatasetFFT)
        lwdatasetBlSNR(iData) = FLW_baseline_SNR.get_lwdata(lwdatasetFFT(iData),option);
    
        freqRes     = lwdatasetBlSNR(iData).header.xstep;
        freqVec     = freqRes : freqRes : lwdatasetBlSNR(iData).header.datasize(6)*freqRes;
        
        %Saving freqRes, freqVec in log
        
        cond = string(intersect(strsplit(condFindTxt),strsplit(lwdatasetBlSNR(iData).header.name, {'_', ' '})));
        Log.(group).(subName).TappingPreprocessing.BlSNR.(cond).freqRes       = freqRes;
        Log.(group).(subName).TappingPreprocessing.BlSNR.(cond).freqVec       = freqVec;
    
    end
    
    fprintf ('  ==> BlSNR saved \n')
end