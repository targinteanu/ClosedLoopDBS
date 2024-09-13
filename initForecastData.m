function foreData = initForecastData(pastData, foreSamples)
% 
% Input: 
%   cell array of past data buffers as timetables 
%   array of # of forecast samples corresponding to above 
% 
% Output: 
%   correctly-sized buffers of zero data as timetables 
% 

foreData = cell(size(pastData));

for ch = 1:length(pastData)
    pastData_ch = pastData{ch};
    bufferSize = height(rawData_ch);
    fs = rawData_ch.Properties.UserData.SampleRate;
    ud.SampleRate = fs;
    foreData{ch} = timetable(...
        zeros(buffersize + foreSamples(ch), width(pastData_ch)), ...
        'SampleRate', fs, ...
        'StartTime', pastData_ch.Time(1), ...
        'VariableNames', pastData_ch.Properties.VariableNames); 
    foreData{ch}.Properties.UserData = ud; 
    foreData{ch}.Properties.VariableUnits = ...
        pastData_ch.Properties.VariableUnits; 
end

end