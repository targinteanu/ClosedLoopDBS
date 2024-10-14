function [storage1, p1, storage2, p2] = ...
    bufferStorage(storage1, p1, storage2, newData)
% Populate a storage buffer storage1 until it is full; then, begin
% populating the next buffer storage2. 
% storage2 is assumed to be empty at function call. 
% p1 and p2 are pointers to the current rows of the respective storages. If
% zero is returned, storage is full. 

N = height(newData); 
if p1+N-1 > height(storage1)
    % storage 1 is now full 
    if N > height(storage2)
        warning('Data overloaded save buffer; some data may not be saved.')
        N = height(storage2);
        newData = newData(1:N, :);
    end
    p1 = 0; 
    storage2(1:N, :) = newData; 
    p2 = N+1;
else
    % normal buffering
    storage1(p1:(p1+N-1), :) = newData;
    p1 = p1+N;
    p2 = [];
end

end