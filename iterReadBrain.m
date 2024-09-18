function [rawData, fltData, fltArgs] = iterReadBrain(...
    rawData, daqFun, ...
    fltData, fltFun, fltArgs)
% One iteration of brain reading, to be used as a main loop or timer
% function. 
% 
% Inputs: 
%   xxxData is a cell with 4 rows: name, head, tail, complete 
%   calling daqFun() gets new tails 
%   [fltTails, fltArgs] = fltFun(fltArgs, rawTails)

rawNames = rawData(1,:);

newTails = daqFun(); 
for ch = 1:size(newTails,2)
    newTail = newTails{ch};
    tailName = newTail.Properties.VariableNames{1};
    CH = find(strcmp(tailName, rawNames)); 
    if isempty(CH)
        error(['Unrecognized new raw data label: ',tailName]);
    end
    if length(CH) > 1
        error(['Raw data label ',tailName,' is not unique.']);
    end

    [rawData{2,CH}, rawData{3,CH}, rawData{4,CH}] = ...
        bufferAndRetime(rawData{2,CH},rawData{3,CH},newTail);
end

rawTails = rawData(3,:); 
[fltTails, fltArgs] = fltFun(fltArgs, rawTails);
for CH = 1:size(fltData,2)
    [fltData{2,CH}, fltData{3,CH}, fltData{4,CH}] = ...
        bufferAndRetime(fltData{2,CH},fltData{3,CH},fltTails{CH});
end

end