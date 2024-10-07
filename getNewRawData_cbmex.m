function [newContData, chname] = getNewRawData_cbmex(chsel)

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

newContData = cell(length(chsel),1);
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