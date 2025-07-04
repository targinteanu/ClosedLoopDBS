function [filtwts, filtorder] = ...
    buildFIRBPF(srate, locutoff, hicutoff, minfac, min_filtorder)
% Build a finite impulse response (FIR) band-pass filter (BPF) for a
% specified signal sample rate and frequency band. Output the filter
% weights, i.e. numerator coefficients that can be used in e.g. 
% x_filt = filter(filtwts, 1, x_unfilt) or filtfilt(filtwts, 1, x_unfilt)

% handle incomplete input args 
if nargin < 5
    min_filtorder = [];
    if nargin < 4
        minfac = [];
    end
end

% filtering bound rules 
if isempty(min_filtorder)
    min_filtorder  = 15;   % minimum filter length
end
if isempty(minfac)
    minfac         = 2;    % this many (lo)cutoff-freq cycles in filter
end

nyq            = srate*0.5;  % Nyquist frequency

% check that specs are proper
if locutoff>0 & hicutoff > 0 & locutoff > hicutoff
    errordlg('locutoff > hicutoff ???', 'filter spec');
    return
end
if locutoff < 0 | hicutoff < 0
    errordlg('locutoff | hicutoff < 0 ???', 'filter spec');
    return
end
if locutoff>nyq
    errordlg('Low cutoff frequency cannot be > srate/2', 'filter spec');
    return
end
if hicutoff>nyq
    errordlg('High cutoff frequency cannot be > srate/2', 'filter spec');
    return
end

% filter order 
if locutoff>0
    filtorder = minfac*fix(srate/locutoff);
elseif hicutoff>0
    filtorder = minfac*fix(srate/hicutoff);
end
if filtorder < min_filtorder
    filtorder = min_filtorder;
end

% build filter 
% usage i.e.: 
% >> filteredSignal = filter(filtwts, 1, unfilteredSignal) 
% -- OR --
% >> filteredSignal = filtfilt(filtwts, 1, unfilteredSignal)
filtwts = fir1(filtorder, [locutoff, hicutoff]./(srate/2));

end