function emptyData = initFilteredData(rawData, filtShift)
% 
% Input: 
%   cell array of unfilt data buffers 
%   array of # of samples eaten by filter corresponding to above 
% 
% Output: 
%   correctly-sized buffers of zero data as timetables 
% 

if length(filtShift) < length(rawData)
    if length(filtShift) == 1
        % assume the one input applies to all channels. 
        filtShift = repmat(filtShift, size(rawData));
    else
        error('Incompatible input dimensions.')
    end
end

emptyData = cell(size(rawData));

for ch = 1:length(rawData)
    rawData_ch = rawData{ch};
    bufferSize = height(rawData_ch);
    if bufferSize > filtShift(ch)
        bufferSize = bufferSize - filtShift(ch);
        t = rawData_ch(1:bufferSize,:);
        emptyData{ch} = [t, ...
            zeros(bufferSize - filtShift(ch), width(rawData_ch)-1)];
    else
        error(['Data buffer is not long enough to be filtered. ',...
               'Try increasing display window or altering filter.'])
    end
end

end