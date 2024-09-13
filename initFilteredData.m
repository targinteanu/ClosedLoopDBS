function filtData = initFilteredData(rawData, filtShift)
% 
% Input: 
%   cell array of unfilt data buffers as timetables 
%   array of # of samples eaten by filter corresponding to above 
% 
% Output: 
%   correctly-sized buffers of zero data as timetables 
% 

filtData = cell(size(rawData));

for ch = 1:length(rawData)
    rawData_ch = rawData{ch};
    bufferSize = height(rawData_ch);
    fs = rawData_ch.Properties.UserData.SampleRate;
    ud.SampleRate = fs;
    if bufferSize > filtShift(ch)
        filtData{ch} = timetable(...
            zeros(bufferSize - filtShift(ch), width(rawData_ch)), ...
            'SampleRate', fs, ...
            'StartTime', rawData_ch.Time(1), ...
            'VariableNames', rawData_ch.Properties.VariableNames);
        filtData{ch}.Properties.UserData = ud;
        filtData{ch}.Properties.VariableUnits = ...
            rawData_ch.Properties.VariableUnits;
    else
        error(['Data buffer is not long enough to be filtered. ',...
               'Try increasing display window or altering filter.'])
    end
end

end