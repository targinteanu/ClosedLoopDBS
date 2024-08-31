function newContData = getNewRawData_cbmex(chsel)

[spikeEvents, time, continuousData] = cbmex('trialdata',1);
chnum = [continuousData{:,1}]';
fs = [continuousData{:,2}]';
chname = spikeEvents(:,1);
chname = chname(chnum);
newContDataRaw = continuousData(:,3); 

newContData = cell(length(chsel),1);
for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch)); 
    if ~isempty(chInd)
        newContData{ch} = timetable(...
            newContDataRaw{chInd}, ...
            'SampleRate', fs(chInd), ...
            'StartTime', time, ...
            'VariableNames', chname(chInd)); 
    end
end

end