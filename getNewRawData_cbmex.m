function [newContData, chname, chnum] = getNewRawData_cbmex(chsel)

% to do: this should also return the channel info for checking 

[spikeEvents, time, continuousData] = cbmex('trialdata',1);

if isempty(continuousData)
    error('No continuous data; ensure that data acquisition has been enabled.')
end

chnum = [continuousData{:,1}]';
%fs = [continuousData{:,2}]';
chname = spikeEvents(:,1);
chname = chname(chnum);
newContDataRaw = continuousData(:,3); 

if isempty(chsel)
    % select all channels
    chsel = chnum;
end

newContData = cell(1,length(chsel));

if width(chsel) == width(chnum)
    % assume order is the same for timing 
    for ch = 1:length(chsel)
        newContData{ch} = ...
            [nan(size(newContDataRaw{ch})), newContDataRaw{ch}];
        newContData{ch}(1,1) = time;
    end
else

% if not same size, can't be same order
for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch)); 
    if ~isempty(chInd)
        newContData{ch} = ...
            [nan(size(newContDataRaw{chInd})), newContDataRaw{chInd}];
        newContData{ch}(1,1) = time;
        %{
        newContData{ch} = timetable(...
            newContDataRaw{chInd}, ...
            'SampleRate', fs(chInd), ...
            'StartTime', seconds(time) + initTime, ...
            'VariableNames', chname(chInd)); 
        %}
    end
end

end

end