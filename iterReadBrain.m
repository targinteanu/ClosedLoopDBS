function [...
    timeBuffs, rawData, ...
    artRemData, artRemArgs, ...
    fltData, fltArgs, ...
    forBuffs, forData, forArgs] = ...
    iterReadBrain(...
        timeBuffs, rawData, daqFun, ...
        selRaw2Flt, selRaw2For, selFlt2For, ...
        artRemData, stimRef, artRemFun, artRemArgs, ...
        fltData, fltFun, fltArgs, ...
        forBuffs, forData, forFun, forArgs)
% One iteration of brain reading, to be used as a main loop or timer
% function. 
% 
% Input Details: 
%   timeBuffs columns correspond to columns of rawData 
%   xxxData is a cell with 4 rows: (1) name, (2) head, (3) tail, (4) complete 
%   selxxx2xxx is a horizontal array of selected columns 
%   calling daqFun() gets new raw tails 
%   stimRef is reference signal for new stimuli that would cause artifact
%   [artRemData, artRemArgs] = artRemFun(artRemArgs, artRemHeads, rawTails)
%   [fltTails, fltArgs] = fltFun(fltArgs, rawTails)
%   [forTails, forBuffsAdd, forArgs] = forFun(forArgs, inData)

%% handle/check inputs 
tic

doDAQ = (~isempty(rawData)) && (~isempty(daqFun));
doArt = doDAQ && (~isempty(artRemData)) && (~isempty(stimRef)) && (~isempty(artRemFun));
doFlt = doDAQ && (~isempty(fltFun)); 
doFor = doDAQ && (~isempty(forFun));

if (~doFlt) && (~isempty(selFlt2For))
    error('Filter required to forecast, but filter is not specified.');
end
if (max(selFlt2For) > size(fltData,2)) || (min(selFlt2For) < 1)
    error('Improper selection of filtered channels.');
end

setup_time = toc
%% DAQ 
tic
if doDAQ 

chInfo = rawData(1,:);
rawNames = cellfun(@(c) c.Name, chInfo, 'UniformOutput',false);
% to do: for speed, use channel ID# instead of names 

[newTails, tailNames] = daqFun(); % to do: daqfun should also get chan info
for ch = 1:size(newTails,2)
    newTail = newTails{ch};
    tailName = tailNames{ch};
    tailProcTime = newTail(end,1);
    CH = find(strcmp(tailName, rawNames)); 
    if isempty(CH)
        error(['Unrecognized new raw data label: ',tailName]);
    end
    if length(CH) > 1
        error(['Raw data label ',tailName,' is not unique.']);
    end

    [rawData{2,CH}, rawData{3,CH}, rawData{4,CH}] = ...
        bufferjuggle(rawData{2,CH},rawData{3,CH},newTail,@bufferData);
    timeBuffs{1,CH} = bufferData(timeBuffs{1,CH}, tailProcTime);
end

% selection, etc
rawTails = rawData(3,:); rawAllData = rawData(4,:);
try
    rawTails = rawTails([selRaw2Flt, selRaw2For]); 
    rawAllData = rawAllData([selRaw2Flt, selRaw2For]);
catch ME
    if strcmp(ME.identifier, 'MATLAB:badsubscript')
        error('Requested selected channel(s) that do(es) not exist in raw data.');
    else
        rethrow(ME);
    end
end
selRaw2Flt = 1:length(selRaw2Flt); 
selRaw2For = length(selRaw2Flt) + (1:length(selRaw2For));
lenRaw = cellfun(@height, rawTails); lenFlt = lenRaw(selRaw2Flt);

end

DAQ_time = toc
%% Artifact Removal 
tic
if doArt

artRemHeads = artRemData(2,:);
[artRemTails, artRemArgs] = artRemFun(artRemArgs, artRemHeads, rawTails);
if ~(size(artRemTails,2) == size(artRemData,2))
    error('Artifact removal channels are inconsistent.');
end
for CH = 1:size(artRemData, 2)
    [artRemData{2,CH}, artRemData{3,CH}, artRemData{4,CH}] = ...
        bufferjuggle(artRemData{2,CH},artRemData{3,CH},artRemTails{CH},@bufferData);
    % does this buffering need to have an overwrite for future extended
    % data instead? 
end
rawTails = artRemData(3,:); rawAllData = rawData(4,:);

end
artifact_time = toc
%% Filter
tic
if doFlt
 
[fltTails, fltArgs] = fltFun(fltArgs, rawTails(selRaw2Flt));
if ~(size(fltTails,2) == size(fltData,2))
    error('Filtered channels are inconsistent.');
end
for CH = 1:size(fltData,2)
    [fltData{2,CH}, fltData{3,CH}, fltData{4,CH}] = ...
        bufferjuggle(fltData{2,CH},fltData{3,CH},fltTails{CH},@bufferData);
end
fltTails = fltData(3,:); fltAllData = fltData(4,:);

else
    fltTails = {}; fltAllData = {};
end

filter_time = toc
%% Forecast 
tic
if doFor

[forTails, forBuffsAdd, forArgs] = forFun(forArgs, ...
    [rawAllData(selRaw2For), fltAllData(selFlt2For)]);
lenFor = [lenRaw(selRaw2For), lenFlt(selFlt2For)];
if ~(size(forTails,2) == size(forData,2))
    error('Forecast channels are inconsistent.');
end
for CH = 1:size(forData,2)
    [forData{2,CH}, forData{3,CH}, forData{4,CH}] = ...
        bufferjuggle(forData{2,CH},forData{3,CH},forTails{CH}, ...
        @(old, new) bufferDataOverwrite(old, new, lenFor(CH)));
    forBuffs{1,CH} = bufferData(forBuffs{1,CH}, forBuffsAdd{1,CH});
end

end

forecast_time = toc

%% helper 
    function [newBuffer, newTail, newAll] = ...
        bufferjuggle(oldBuffer, oldTail, newData, bufferFunc)
        newBuffer = bufferFunc(oldBuffer, oldTail); 
        newTail = newData; 
        newAll = bufferFunc(newBuffer, newTail); 
    end

end