# pattern_repetition


This is the online repository for the code used in the paper entitled "Repetition of Rhythmic Patterns Fosters Neural Representation of Musical Meter" by E. Coulon, S. Baum, T. Lenc, R. Polak and S. Nozaradan.



Data Description
Data is available on the Zenodo platform  

The raw data (81GB) can be sent upon request.

The output of the cochlear model is available under the name UREAR_BMF[64]_*_F*_noDev.mat for each condition and each pure tone carrier frequency. These files should be put in the folder 1_stimulus.

The preprocessing is already done and the preprocessed dataset is available under the name "PreprocessedEEG_Musicians" for the group of musician participants and "PreprocessedEEG_NonMusicians" for the group of non-musician participants. For the script to run correctly, these two Matlab structures should be put in the 2_preprocessing folder.

The spectral analysis is already done and the analyzed EEG dataset is available under the name "SpectralAnalysis_Musicians" for the group of musician participants and "SpectralAnalysis_NonMusicians" for the group of non-musician participants. For the script to run correctly, these two Matlab structures should be put in the 3_analysis folder.



How to run the code
Please download letswave6 (https://github.com/NOCIONS/letswave6.git) and letswave7 (https://github.com/NOCIONS/letswave7.git) and put them in the lib folder to run the scripts. 

The preprocessing codes are currently available for reviewing only. To run the code correctly, it necessitates the raw data available upon request. 

For the EEG spectral analysis, the output of the fft section (between line 57 and 84) is available in SpectralAnalysis_Musicians and SpectralAnalysis_NonMusicians files in the structure AnalyzedEEG_*_blsnrFFT. To run the code AnalysisEEG.m, please comment this section and the rest should work fine. 



