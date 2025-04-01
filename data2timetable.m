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
        timestamps = find(timestamps);
        % project/estimate time of t1; although strictly speaking, this is
        % not accurate and could appear misleading
        if isnan(t(1)) % almost always true
            if numel(timestamps)
                n2 = timestamps(1);
                T12 = (n2-1)/chinfo.SampleRate;
                t(1) = t(n2) - T12;
            end
        end
        timestamps = [1; timestamps; length(t)];
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
        clear timestamps
    
        TT.Properties.VariableUnits = {Unit};
        TT.Properties.UserData = chinfo;
        % retime to const sample rate and/or sort rows??

    end
    TTs{ch} = TT;
end

end