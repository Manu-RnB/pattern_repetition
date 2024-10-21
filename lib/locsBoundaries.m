function [locs,Details,varargout] = locsBoundaries(locs,varargin)

    % 1) Set a lower and upper limit to the detection of peak location in a spectrum. 
    % 2) Remove the 5Hz and harmonics
    
    
%% Initial setup

    if nargin == 2
        
        Options     = varargin{1};
        optionNames = fieldnames(Options);
        
        % remove the selected indeces from the peaks as well
        if ~isempty(find(strcmp(optionNames,'pks')))
            applyOnPks      = true;
            pks             = Options.pks;
        else
            applyOnPks      = false;
        end
        
        % remove locs found lower then a certain frequency
        if ~isempty(find(strcmp(optionNames,'lowerLimit')))
            applyLowerLimit = true;
            lowerLimit      = Options.lowerLimit;
        else
            applyLowerLimit = false;
        end
        
        % remove locs found higher then a certain frequency
        if ~isempty(find(strcmp(optionNames,'upperLimit')))
            applyUpperLimit = true;
            upperLimit      = Options.upperLimit;
        else
            applyUpperLimit = false;
        end  
        
        % remove 5Hz and harmonics
        if ~isempty(find(strcmp(optionNames,'remove5Hz')))
            remove5Hz = Options.remove5Hz;
        else
            remove5Hz = false;
        end
        
    else
        applyOnPks      = false;
        applyLowerLimit = false;
        applyUpperLimit = false;
        remove5Hz       = false;
    end
        
    
%% set range boundaries for locs

    locs = unique(locs);
    
    % remove locs found lower then a certain frequency
    if applyLowerLimit
        idx2remove       = find(locs<lowerLimit);
        locs(idx2remove) = [];
        if applyOnPks;  pks (idx2remove) = [];  end
    end
    clear idx2remove   
    
    % remove locs found higher then a certain frequency
    if applyUpperLimit
        idx2remove       = find(locs>upperLimit);
        locs(idx2remove) = [];
        if applyOnPks;  pks (idx2remove) = [];  end
    end
    clear idx2remove 
    
    % remove 5Hz and harmonics
    if remove5Hz
        i = 1;
        for iLocs = 1:length(locs)
            if mod(round(locs(iLocs),2),5) == 0 % if multiple of 5Hz   
                idx2remove(i) = iLocs;
                i = i+1;
            end
        end
        if exist('idx2remove','var')            
            locs(idx2remove) = [];
            if applyOnPks;  pks(idx2remove) = [];  end
        end
    end
    
%% save details

Details.ApplyOnPks      = applyOnPks;
Details.ApplyLowerLimit = applyLowerLimit;
Details.ApplyUpperLimit = applyUpperLimit;
Details.Remove5Hz       = remove5Hz;

if applyOnPks
    varargout = {pks};
end

clear Options

end

