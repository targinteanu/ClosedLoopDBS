function newContData = getNewRawData_cbmex(chsel)

[spikeEvents, time, continuousData] = cbmex('trialdata',1);
continuousData = continuousData(chsel,:); % subselect good channels of interest
chnum = [continuousData{:,1}]';
fs = [continuousData{:,2}]';
chname = spikeEvents(:,1);
chname = chname(chnum);
newContDataRaw = continuousData(:,3); 

newContData = cellfun(...
    @(ncd,f,n) timetable(ncd, 'SampleRate',f, 'VariableNames',{n}), ...
    newContDataRaw, fs, chname, 'UniformOutput', false);

end