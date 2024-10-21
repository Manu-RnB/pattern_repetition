
clear all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rename the files using the "old" rereferencing names %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('/Users/emmanuelcoulon/Documents/MATLAB/PROJECTS/NonRepComparison/Datasets/LW/NonMusicians/Preprocessing/EEG/ICA/MultiRereference')

for iRerefMethods = 1:3
    
    oldRerefNames = {'reref_mast','reref_allInclMast','reref_allExclMast'};
    newRerefNames = {'rerefMast';'rerefAllInclMast';'rerefAllExclMast'};
    
    
    Files = dir(['*',oldRerefNames{iRerefMethods},'*','sub','*']);

    for iFiles = 1:length(Files)

        oldName = Files(iFiles).name;
        newname = strrep(oldName, ...
                        oldRerefNames{iRerefMethods}, ...
                        newRerefNames{iRerefMethods});
        
        movefile(oldName,newname);
                    
    end
end

