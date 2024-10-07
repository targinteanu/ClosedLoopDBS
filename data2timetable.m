function TTs = data2timetable(datas, chaninfos, initTime)
% Build a cell of timetables using a cell of data in the form [times, data]
% for each channel and an accompanying cell of channel info. 
% Times are nan when there is not a time stamp associated. 
% channel info has fields SampleRate, Name, Unit
% initTime must be a datetime 

TTs = cell(size(datas));

for ch = 1:size(datas,2)
    data = datas{ch};
    chinfo = chaninfos{ch};
    t = data(:,1); 

    timestamps = ~isnan(t);
    timestamps = [1; find(timestamps); length(t)];
    TT = [];
    for tsi = 1:(length(timestamps)-1)
        n1 = timestamps(tsi); n2 = timestamps(tsi+1)-1;
        Di = data(n1:n2,2);
        TTi = timetable(Di, ...
            'SampleRate', chinfo.SampleRate, ...
            'StartTime', seconds(t(n1)) + initTime, ...
            'VariableNames', {chinfo.Name});
        TT = [TT; TTi];
    end

    TT.Properties.VariableUnits = {chinfo.Unit};
    TT.Properties.UserData = chinfo;
    % retime to const sample rate and/or sort rows??
    TTs{ch} = TT;
end

end