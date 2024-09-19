function newBuffer = bufferDataHorizon(tHorizon, oldBuffer, newData)
% Buffer data once it has passed the horizon. If data has not passed the
% horizon, it can be overwritten. 

%% check inputs 

[oldBuffer, norefmt] = reformatBufferData(oldBuffer); 
newData = reformatBufferData(newData);

N = height(newData); W = width(newData); H = height(oldBuffer);
if W~= width(oldBuffer)
    error('Mismatching channels of data to buffer.')
end

%% main function 
newBuffer = oldBuffer;
for c = 1:W
    ob = oldBuffer{c}; nd = newData{c}; 
    obt = oldBuffer.Time; ndt = newData.Time; 
    is_safe_new = ndt <= tHorizon; is_out_old = obt > tHorizon;
    T0 = bufferData(ob(~is_out_old,:), nd(is_safe_new,:));
    n = N - sum(is_safe_new); 
    T1 =  ob(is_out_old,:); T2 =  nd(~is_safe_new,:);
    t1 = obt(is_out_old,:); t2 = ndt(~is_safe_new,:);
    if n > sum(is_out_old)
        error('Too much new data.')
    else
        if n == sum(is_out_old)
            T1 = T2;
        else
            for i2 = 1:length(t2)
                [~,i1] = min(abs(t1 - t2(i2)));
                T1 = [T1(1:(i1-1),:); T2(i2,:); T1((i1+1):end,:)];
            end
        end
    end
    newBuffer{c} = [T0; T1];
end

%% format outputs in the same format as old buffer
if norefmt
    % leave as timetable
    newBuffer = newBuffer{1};
else
    % reformat to array
    nb = seconds(zeros(H,W)); 
    for c = 1:W
        nb(:,c) = newBuffer{c}.Time;
    end
    newBuffer = nb;
end

%% helpers 
    function [cT, refd] = reformatBufferData(bufferData)
        % if not a timetable, assume it is a series of time values (e.g.
        % time of peaks and troughs) in seconds.
        refd = istimetable(bufferData);
        if refd
            cT = {bufferData};
        else
            cT = cell(1,width(bufferData));
            for col = 1:width(cT)
                t = seconds(bufferData(:,col));
                cT{1,col} = array2timetable(true(size(t)), 'RowTimes',t);
            end
        end
    end

end