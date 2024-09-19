function [newBuffTime, newBuffVal] = bufferDataHorizon(...
    tHorizon, oldBuffTime, newDatTime, oldBuffVal, newDatVal)
% Buffer data once it has passed the horizon. If data has not passed the
% horizon, it can be overwritten. 

%% check inputs 
if nargin < 5
    newDatVal = [];
    if nargin < 4
        oldBuffVal = [];
    end
end
if isempty(oldBuffVal)
    oldBuffVal = oldBuffTime;
end
if isempty(newDatVal)
    newDatVal = newDatTime;
end

N = height(newDatTime); W = width(newDatTime); H = height(oldBuffTime);

if (H ~= height(oldBuffVal))
    error('Mismatching buffer time and value.')
end
if (N ~= height(newDatVal))
    error('Mismatching data time and value.')
end

timerep = width(oldBuffTime) == 1;
if (width(oldBuffTime) ~= W)
    if timerep
        oldBuffTime = repmat(oldBuffTime, 1, width(newDatTime));
    else
        error('Mismatching buffer time sizes.')
    end
end

%% main function 
for c = 1:W
    obt = oldBuffTime(:,c); ndt = newDatTime(:,c);
    is_ndt_safe = ndt <= tHorizon; 
    n = N - sum(is_ndt_safe);
    obt = bufferData(obt, ndt(is_ndt_safe));
end

%% edit outputs 
if timerep
    newBuffTime = newBuffTime(:,1);
end

end