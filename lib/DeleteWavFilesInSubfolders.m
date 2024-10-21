% This scripts goes through all wavfiles in subfolders and deletes them % 

% used to remove the wavfiles automatically saved when recording the PTB
% scripts

rootdir     = '/Users/emmanuelcoulon/Documents/PTB/NonRepComplexPTB/output/source';
filelist    = dir(fullfile(rootdir, '**/*.wav'));   %get list of files and folders in any subfolder
filelist    = filelist(~[filelist.isdir]);          %remove folders from list

for i = 1:length(filelist)
    delete (fullfile(filelist(i).folder,filelist(i).name))
end