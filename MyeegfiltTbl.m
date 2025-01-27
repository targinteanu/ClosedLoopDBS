function [smoothTbl, filtwts] = MyeegfiltTbl(dataTbl, locutoff, hicutoff, ...
    epochframes,filtorder,revfilt,firtype,causal)
% 
% Performs the function Myeegfilt to bandpass-filter data in a timetable. 
% Necessary inputs: 
%   dataTbl: unfiltered data in a timetable 
%   locutoff: low-frequency cutoff (Hz)
%   hicutoff: high-frequency cutoff (Hz)
% For use of further optional inputs, see Myeegfilt documentation. 
% Outputs: 
%   smoothTbl: filtered data in a timetable 
%   filtwts: fir filter weights, such that filteredData = filter(filtwts, 1, rawData)
% 

data = dataTbl.Variables';
srate = dataTbl.Properties.SampleRate; 

if nargin < 4
    [smoothdata, filtwts] = Myeegfilt(data, srate, locutoff, hicutoff);
elseif nargin < 5
    [smoothdata, filtwts] = Myeegfilt(data, srate, locutoff, hicutoff, ...
        epochframes);
elseif nargin < 6
    [smoothdata, filtwts] = Myeegfilt(data, srate, locutoff, hicutoff, ...
        epochframes, filtorder);
elseif nargin < 7
    [smoothdata, filtwts] = Myeegfilt(data, srate, locutoff, hicutoff, ...
        epochframes, filtorder, revfilt);
elseif nargin < 8
    [smoothdata, filtwts] = Myeegfilt(data, srate, locutoff, hicutoff, ...
        epochframes, filtorder, revfilt, firtype);
else
    [smoothdata, filtwts] = Myeegfilt(data, srate, locutoff, hicutoff, ...
        epochframes, filtorder, revfilt, firtype, causal);
end

smoothTbl = dataTbl;
smoothTbl.Variables = smoothdata';

end