function TTfilt = FilterTimetable(FiltObj, TTunfilt, FiltFun)
%
% Filter time-series data in a timetable 
%
% Inputs: 
%   FiltFun: function [e.g. @(d,x) filtfilt(d,x)] of form 
%            FilteredSignal = FiltFun(FiltObj, UnfilteredSignal) 
%            where FilteredSignal and UnfilteredSignal are matrix/vector 
%            columns over time.
%   FiltObj: digitalFilter object 
%   TTunfilt: unfiltered data in timetable format 
%
% Outputs: 
%   TTfilt: filtered data in timetable format 
%

if nargin < 3
    FiltFun = [];
end
if isempty(FiltFun)
    FiltFun = @(d,x) filter(d,x);
end

%% sample rate check
fs_Filter = FiltObj.SampleRate; 
fs_Signal = TTunfilt.Properties.SampleRate;
Time_Signal = TTunfilt.Time; 
errthresh = .01;
if abs((fs_Filter-fs_Signal)/fs_Signal) > errthresh
    error('Filter and input signal have different sampling rates.')
end

%% filtering 
Xunfilt = table2array(TTunfilt);
Xfilt = FiltFun(FiltObj, Xunfilt);
TTfilt = array2timetable(Xfilt,"RowTimes",Time_Signal,...
    "VariableNames",TTunfilt.Properties.VariableNames); 
TTfilt.Properties.VariableUnits = TTunfilt.Properties.VariableUnits;
TTfilt.Properties.UserData = TTunfilt.Properties.UserData;

end