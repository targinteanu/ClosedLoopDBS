function TTs = data2timetable(datas, chaninfos, initTime)
% Build a cell of timetables using a cell of data in the form [times, data]
% for each channel and an accompanying cell of channel info. 
% Times are nan when there is not a time stamp associated. 
% channel info has fields SampleRate, Name, Unit
% initTime must be a datetime 

TTs = cell(size(datas));

for ch = 1:size(datas,2)
    data = datas{ch};
    TT = [];
    if numel(data)

        chinfo = chaninfos{ch};
        Unit = chinfo.Unit;
        t = data(:,1); 

        scaleunits = ~isnan(chinfo.Resolution);
        if ~scaleunits
            Unit = [Unit,'[?]']; % indicate uncertainty 
        end
    
        timestamps = ~isnan(t);
        timestamps = [1; find(timestamps); length(t)];
        % *** Is this the cause of x axis jumping around? First block may
        % not have timestamp? ***
        for tsi = 1:(length(timestamps)-1)
            n1 = timestamps(tsi); n2 = timestamps(tsi+1)-1;
            Di = double(data(n1:n2,2));
            if scaleunits
                Di = Di - chinfo.MinDigital; 
                Di = Di .* chinfo.Resolution; 
                Di = Di + chinfo.MinAnalog;
            end
            TTi = timetable(Di, ...
                'SampleRate', chinfo.SampleRate, ...
                'StartTime', seconds(t(n1)) + initTime, ...
                'VariableNames', {chinfo.Name});
            TT = [TT; TTi];
        end
    
        TT.Properties.VariableUnits = {Unit};
        TT.Properties.UserData = chinfo;
        % retime to const sample rate and/or sort rows??

    end
    TTs{ch} = TT;
end

end