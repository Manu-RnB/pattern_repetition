function [freqVec,sFFT] = FFT(s,fs,varargin)


% Function to compute the fft 


%%%%%%%%%%%%%%%%%%
% Initial setups % 
%%%%%%%%%%%%%%%%%%

    % Load optional arguments
    if nargin == 3
        Options     = varargin{1};
        optionNames = fieldnames(Options);

        if ~isempty(find(strcmp(optionNames, 'hilb')))
            hilb = Options.hilb;
        else
            hilb = 'No';
        end
        
    elseif nargin == 2
        hilb = 'No';
    end

%%%%%%%%%%%
% Options %  
%%%%%%%%%%%

    % Hilbert transform
    if strcmpi(hilb, 'Yes')
        s = abs(hilbert(s));
    end
    
%%%%%%%%%%%%%%%    
% Compute fft %    
%%%%%%%%%%%%%%% 

    n       = length(s);
    freqVec = linspace(0,fs/2,floor(n/2)+1);

    sFFT     = 2*abs(fft(s))/n;

    % Only keep positive frequencies
    sFFT = sFFT(1:length(freqVec));
    
    % Remove DC
    sFFT(1) = 0;

end


% plot(linspace(0,srate, length(x)),abs(fft(x)) => comprend les fréquences
% négatives
% => set (gca'xlim',[0,srate/2])
