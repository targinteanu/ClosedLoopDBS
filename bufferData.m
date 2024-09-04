function newBuffer = bufferData(oldBuffer, newData)

N = height(newData); 
if N >= height(oldBuffer)
    newBuffer = newData(end-length(oldBuffer)+1:end, :);
    newBuffer.Properties.VariableUnits = oldBuffer.Properties.VariableUnits;
else
    newBuffer = [oldBuffer(N+1:end, :); newData];
end

end