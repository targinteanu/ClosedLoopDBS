function [TTfilt, filtCond] = FilterTimetable(b,a, TTunfilt, FiltCond, FiltFun)
%
% Filter time-series data in a timetable 
%
% Inputs: 
%   FiltFun: function [e.g. [xFilt, filtCond] = @(b,a,x) filtfilt(b,a,x)] of form 
%            FilteredSignal = FiltFun(b,a, UnfilteredSignal) 
%            where FilteredSignal and UnfilteredSignal are matrix/vector 
%            columns over time.
%   FiltObj: digitalFilter object 
%   TTunfilt: unfiltered data in timetable format 
%
% Outputs: 
%   TTfilt: filtered data in timetable format 
%

if nargin < 5
    FiltFun = @filter;
end

%{
%% sample rate check
fs_Filter = FiltObj.SampleRate; 
fs_Signal = TTunfilt.Properties.SampleRate;
Time_Signal = TTunfilt.Time; 
errthresh = .01;
if abs((fs_Filter-fs_Signal)/fs_Signal) > errthresh
    error('Filter and input signal have different sampling rates.')
end
%}

%% filtering 
Xunfilt = table2array(TTunfilt);
[Xfilt,filtCond] = FiltFun(b, a, Xunfilt, FiltCond);
TTfilt = array2timetable(Xfilt,"RowTimes",TTunfilt.Time,...
    "VariableNames",TTunfilt.Properties.VariableNames); 
TTfilt.Properties.VariableUnits = TTunfilt.Properties.VariableUnits;
TTfilt.Properties.UserData = TTunfilt.Properties.UserData;

end