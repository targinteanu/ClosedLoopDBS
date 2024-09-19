function [emptyData, contData, buffData] = initRawData_cbmex(chsel, bufferSize, initTime)

[spikeEvents, time, continuousData] = cbmex('trialdata',1);

if isempty(continuousData)
    error('No continuous data; ensure that data acquisition has been enabled.')
end

chnum = [continuousData{:,1}]';
fs = [continuousData{:,2}]';
chname = spikeEvents(:,1);
chname = chname(chnum);

if isempty(chsel)
    % select all channels
    chsel = chnum;
end

emptyData = cell(1,length(chsel)); 
contData = emptyData; 
buffData = emptyData;
for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch)); 

    if ~isempty(chInd)

        % Create raw data buffer of zeros of the correct length
        L = length(continuousData{chInd,3}); 
        emptyData{ch} = timetable(...
            zeros(bufferSize,1), ...
            'SampleRate', fs(chInd), ...
            'StartTime', seconds(time - (bufferSize)/fs(chInd)) + initTime, ...
            'VariableNames', chname(chInd)); 
        contData{ch} = timetable(...
            continuousData{chInd,3}, ...
            'SampleRate', fs(chInd), ...
            'StartTime', seconds(time) + initTime, ...
            'VariableNames', chname(chInd));

        % User Data 
        ud.SampleRate = fs(chInd);
        emptyData{ch}.Properties.UserData = ud;
        contData{ch}.Properties.UserData = ud;

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
        emptyData{ch}.Properties.VariableUnits = {unitname};
        contData{ch}.Properties.VariableUnits = {unitname};

        buffData{ch} = bufferData(emptyData{ch}, contData{ch});

    end
    
end 

end