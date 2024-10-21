function [frexIdx,metRelIdx,metUnrelIdx] = freqClassifier(freqVec,locs,frex,varargin)

% determines whether a frequency is meter-related or meter-unrelated

    if nargin > 3

        Options     = varargin{1};
        optionNames = fieldnames(Options);

        if ~isempty(find(strcmp(optionNames,'GroupingBy8inMetRel')))
            GroupingBy8inMetRel = Options.GroupingBy8inMetRel;
        else
            GroupingBy8inMetRel = false;
        end 
    else
        GroupingBy8inMetRel = false;
    end

    iMetRel   = 1;
    iMetUnrel = 1;
    frexIdx   = dsearchn(freqVec',locs');
    
    for iLocs = 1:length(locs)
        
        % if multiple of 0.625Hz
        if GroupingBy8inMetRel
            if mod(round(locs(iLocs),3),frex(3)/2) == 0 
                metRelIdx(iMetRel)     = iLocs;
                iMetRel                = iMetRel +1;                
            else
                metUnrelIdx(iMetUnrel) = iLocs;
                iMetUnrel              = iMetUnrel +1;
            end
            
        % if multiple of 1.25Hz
        else
            if mod(round(locs(iLocs),3),frex(3)) == 0    
                metRelIdx(iMetRel)     = iLocs;
                iMetRel                = iMetRel +1;
            else
                metUnrelIdx(iMetUnrel) = iLocs;
                iMetUnrel              = iMetUnrel +1;
            end
        end
    end   
    
    % create empty output if no metRel or metUnrel index was found
    if ~exist('metRelIdx');     metRelIdx = [];     end
    if ~exist('metUnrelIdx');   metUnrelIdx = [];   end
    
end

