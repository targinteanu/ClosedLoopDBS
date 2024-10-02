function newBuffer = bufferDataOverwrite(oldBuffer, newData, N)
% Allow the tail end of the old buffer to be overwritten by new data. 
% N = # of points of data that is actually new 

if N <= height(newData)
    if N >= height(oldBuffer)
        newBuffer = newData(end-height(oldBuffer)+1:end, :);
        newBuffer.Properties.VariableUnits = oldBuffer.Properties.VariableUnits;
        newBuffer.Properties.UserData = oldBuffer.Properties.UserData;
    else
        % buffer AND overwrite 
        L = height(newData)-N; % length to overwrite
        newBuffer = [oldBuffer((N+1):(end-L), :); newData];
    end
else
    % there is no data to overwrite; in fact, there is not enough new data
    % nan-pad newData to length N and cycle buffer 
    dT = newData.Properties.TimeStep; 
    nanpad = array2timetable(nan(N-height(newData), width(newData)), ...
        'TimeStep', dT, ...
        'StartTime', newData.Time(end)+dT);
    nanpad.Properties.VariableNames = newData.Properties.VariableNames;
    newData = [newData; nanpad];
    newBuffer = bufferData(oldBuffer, newData);
end

end