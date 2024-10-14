function [storageFull, curStorage, curRow, oldStorage] = ...
    bufferStorage(curStorage, curRow, newData)
% Populate the current storage buffer until it is full; then, begin
% populating the next buffer. 
% curRow indicates current row of the storage. 
% If storageFull is true, storage has been filled, and the old storage to
% be saved is returned as oldStorage. Otherwise, oldStorage is empty. 

N = height(newData); 
storageFull = curRow+N-1 > height(curStorage);
if storageFull
    % storage 1 is now full 
    oldStorage = curStorage;
    curStorage = nan(size(curStorage));
    if N > height(curStorage)
        warning('Data overloaded save buffer; some data may not be saved.')
        N = height(curStorage);
        newData = newData(1:N, :);
    end
    curStorage(1:N, :) = newData; 
    curRow = N+1;
else
    % normal buffering
    curStorage(curRow:(curRow+N-1), :) = newData;
    curRow = curRow+N;
    oldStorage = [];
end

end