function [lwdataset, lwdatasetFFT, Log] = tapFFT(lwdataset, group, subName, ~, Log, Paths)

fprintf ('  -> FFT \n')

if ~exist(fullfile(Paths.LW,group,'Preprocessing/Tapping/FFT'))
    mkdir(fullfile(Paths.LW,group,'Preprocessing/Tapping/FFT'));
end
cd(fullfile(Paths.LW,group,'Preprocessing/Tapping/FFT'));

%%%%%%%%%%%%%%%%%%%%%%%
% FFT %
%%%%%%%%%%%%%%%%%%%%%%%
condFindTxt = 'ABAB1 ABAB2 AAAA CDEF';
option = struct('output','amplitude','half_spectrum',1,'suffix','fft','is_save',1);
for iData = 1:length(lwdataset)
    lwdatasetFFT(iData) = FLW_FFT.get_lwdata(lwdataset(iData),option);

    freqRes = lwdatasetFFT(iData).header.xstep;
    freqVec = freqRes : freqRes : lwdatasetFFT(iData).header.datasize(6)*freqRes;

    %Saving freqRes, freqVec in log
    
    cond = string(intersect(strsplit(condFindTxt),strsplit(lwdatasetFFT(iData).header.name, {'_', ' '})));
    Log.(group).(subName).TappingPreprocessing.FFT.(cond).freqRes       = freqRes;
    Log.(group).(subName).TappingPreprocessing.FFT.(cond).freqVec       = freqVec;

end

fprintf ('  ==> FFT saved \n')

end