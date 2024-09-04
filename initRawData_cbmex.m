function contData = initRawData_cbmex(chsel, bufferSize, initTime)

[spikeEvents, time, continuousData] = cbmex('trialdata',1);

if isempty(continuousData)
    error('No continuous data; ensure that data acquisition has been enabled.')
end

chnum = [continuousData{:,1}]';
fs = [continuousData{:,2}]';
chname = spikeEvents(:,1);
chname = chname(chnum);

contData = cell(length(chsel),1);
for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch)); 

    if ~isempty(chInd)

        % Create raw data buffer of zeros of the correct length
        contData{ch} = timetable(...
            zeros(bufferSize,1), ...
            'SampleRate', fs(chInd), ...
            'StartTime', seconds(time - bufferSize/fs(chInd)) + initTime, ...
            'VariableNames', chname(chInd)); 

        % check units 
        config = cbmex('config', chInd);
        unitname_is = lower(config{11,1});
        if contains(unitname_is, 'unit')
            unitname = config{11,2};
        else
            unitname_is = contains(lower(config(:,1)), 'unit');
            unitname_is = find(unitname_is);
            unitname_is = unitname_is(1); 
            unitname = config{unitname_is,2};
        end
        contData{ch}.Properties.VariableUnits = {unitname};

    end
    
end 

end