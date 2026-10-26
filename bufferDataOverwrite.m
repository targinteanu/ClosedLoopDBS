function newBuffer = bufferDataOverwrite(oldBuffer, newData, N)
% Allow the tail end of the old buffer to be overwritten by new data. 
% N = # of points of data that is actually new 
% TO DO: simplify logic cases since some branches turn out to be the same

isT = istimetable(oldBuffer) || istable(oldBuffer);

if N <= height(newData)
    % there is enough new data 
    if N >= height(oldBuffer)
        % all of old buffer is to be buffered out 
        newBuffer = newData(end-height(oldBuffer)+1:end, :);
        if isT
            newBuffer.Properties.VariableUnits = oldBuffer.Properties.VariableUnits;
            newBuffer.Properties.UserData = oldBuffer.Properties.UserData;
        end
    else
        % some, but not all, of old buffer is to be buffered out 
        L = height(newData)-N; % length to overwrite
        %%{ 
        if L + N >= height(oldBuffer)
            % ALL of old buffer will be replaced, so new data must also be
            % truncated 
            newData = newData(end-height(oldBuffer)+1:end, :);
            % this is the same case as above, so can this be simplified??
        end
        %} 
        % overwrite tail end of old buffer with new data
        newBuffer = [oldBuffer((N+1):(end-L), :); newData];
    end
else
    % there is not enough new data to overwrite old data, 
    % or to buffer without first padding
    % nan-pad newData to length N and cycle buffer 
    nanpad = nan(N-height(newData), width(newData));
    if isT
        dT = newData.Properties.TimeStep; 
        nanpad = array2timetable(nanpad, ...
            'TimeStep', dT, ...
            'StartTime', newData.Time(end)+dT);
        nanpad.Properties.VariableNames = newData.Properties.VariableNames;
    end
    newData = [nanpad; newData];
    newBuffer = bufferData(oldBuffer, newData);
end

end