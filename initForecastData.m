function [emptyPast, emptyFore, emptyBuff] = initForecastData(pastData, foreSamples)
% 
% Input: 
%   cell array of past data buffers as timetables 
%   array of # of forecast samples corresponding to above 
% 
% Output: 
%   correctly-sized buffers of zero data as timetables 
% 

emptyPast = cellfun(@(T) 0.*T, pastData, 'UniformOutput',false);

emptyFore = cell(size(pastData));
emptyBuff = emptyFore;

for ch = 1:length(pastData)
    pastData_ch = pastData{ch};
    bufferSize = height(rawData_ch);
    fs = rawData_ch.Properties.UserData.SampleRate;
    ud.SampleRate = fs;
    emptyFore{ch} = timetable(...
        zeros(bufferSize + foreSamples(ch), width(pastData_ch)), ...
        'SampleRate', fs, ...
        'StartTime', pastData_ch.Time(1), ...
        'VariableNames', pastData_ch.Properties.VariableNames); 
    emptyFore{ch}.Properties.UserData = ud; 
    emptyFore{ch}.Properties.VariableUnits = ...
        pastData_ch.Properties.VariableUnits; 
    emptyBuff{ch} = bufferData(emptyPast{ch}, emptyFore{ch});
end

end