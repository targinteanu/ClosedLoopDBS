function newBuffer = bufferData(oldBuffer, newData)

N = height(newData); 
if N >= height(oldBuffer)
    % all data is new
    newBuffer = newData(end-length(oldBuffer)+1:end, :);
    newBuffer.Properties.VariableUnits = oldBuffer.Properties.VariableUnits;
    newBuffer.Properties.UserData = oldBuffer.Properties.UserData;
else
    % only tail is new 
    newBuffer = [oldBuffer(N+1:end, :); newData];
end

end