function [newBuffer, oldData] = bufferData(oldBuffer, newData)

N = height(newData); 
if N >= height(oldBuffer)
    % all data is new
    oldData = oldBuffer;
    newBuffer = newData(end-height(oldBuffer)+1:end, :);
    if istimetable(newBuffer) || istable(newBuffer)
        newBuffer.Properties.VariableUnits = oldBuffer.Properties.VariableUnits;
        newBuffer.Properties.UserData = oldBuffer.Properties.UserData;
    end
else
    % only tail is new 
    oldData = oldBuffer(1:N, :);
    newBuffer = [oldBuffer(N+1:end, :); newData];
end

end