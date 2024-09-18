function [newBuffer, newTail, newAll] = ...
    bufferAndRetime(oldBuffer, oldTail, newData, bufferFunc)
% Assume start time of newData is ground truth. Retime oldTail if
% necessary. 

if nargin < 4
    bufferFunc = @bufferData; 
end

t1 = oldTail.Time(end); t2 = newData.Time(1); 
if t2 <= t1
    %seconds(t2-t1)
    T = t2 - oldTail.Time(1);
    %{
    oldTail = retime(oldTail,'regular','nearest', ...
        'TimeStep',T/height(oldTail));
    %}
    oldTail.Properties.TimeStep = T/height(oldTail);
    %seconds(t2-oldTail.Time(end))
end

newBuffer = bufferFunc(oldBuffer, oldTail); 
newTail = newData; 
newAll = bufferFunc(newBuffer, newTail); 
newAll = retime(newAll,'regular','nearest',...
    'SampleRate',newAll.Properties.UserData.SampleRate);

end