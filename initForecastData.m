function [emptyPast, emptyFore, emptyBuff] = initForecastData(pastData, foreSamples)
% 
% Input: 
%   cell array of past data buffers as timetables 
%   array of # of forecast samples corresponding to above 
% 
% Output: 
%   correctly-sized buffers of zero data as timetables 
% 

if length(foreSamples) < length(pastData)
    if length(foreSamples) == 1
        % assume the one input applies to all channels. 
        foreSamples = repmat(foreSamples, size(pastData));
    else
        error('Incompatible input dimensions.')
    end
end

emptyPast = cellfun(@(D) [1,0].*D, pastData, 'UniformOutput',false);

emptyFore = cell(size(pastData));
emptyBuff = emptyFore;

for ch = 1:length(pastData)
    pastData_ch = pastData{ch};
    emptyFore{ch} = [nan(foreSamples(ch),1), ...
        zeros(foreSamples(ch), width(pastData_ch)-1)]; 
    emptyBuff{ch} = [emptyPast{ch}; emptyFore{ch}];
end

end