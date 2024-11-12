function [...
    timeBuffs, rawData, ...
    artRemData, artRemArgs, ...
    fltData, fltArgs, ...
    forBuffs, forBuffRow, forBuffedOut, forData, forArgs] = ...
    iterReadBrain(...
        timeBuffs, rawData, daqFun, ...
        selRaw2Flt, selRaw2For, selFlt2For, ...
        artRemData, stimRef, artRemFun, artRemArgs, ...
        fltData, fltFun, fltArgs, ...
        forBuffs, forBuffRow, forData, forFun, forArgs)
% One iteration of brain reading, to be used as a main loop or timer
% function. 
% 
% Input Details: 
%   timeBuffs columns correspond to columns of rawData 
%   xxxData is a cell with 4 rows: (1) name, (2) head, (3) tail, (4) complete 
%   selxxx2xxx is a horizontal array of selected column indexes
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
if numel(selFlt2For)
    if (max(selFlt2For) > size(fltData,2)) || (min(selFlt2For) < 1)
        error('Improper selection of filtered channels.');
    end
end

setup_time = toc; disp(['Setup Time = ',num2str(setup_time)])
%% DAQ 
tic
if doDAQ 

chInfo = rawData(1,:);
%rawNames = cellfun(@(c) c.Name, chInfo, 'UniformOutput',false);
rawIDs = cellfun(@(c) c.IDnumber, chInfo);

curTime = nan(1,size(rawData,2));

[newTails, tailNames, tailIDs] = daqFun(); 
for ch = 1:size(newTails,2)
    newTail = newTails{ch};
    tailName = tailNames{ch};
    %CH = find(strcmp(tailName, rawNames)); 
    tailID = tailIDs(ch);
    CH = find(tailID == rawIDs);
    if isempty(CH)
        error(['Unrecognized new raw data label: ',tailName]);
    end
    if length(CH) > 1
        error(['Raw data label ',tailName,' is not unique.']);
    end

    fs = chInfo{CH}.SampleRate;
    tailProcTime = newTail(1,1) + (height(newTail)-1)/fs; 
    curTime(CH) = tailProcTime;

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

DAQ_time = toc; disp(['DAQ time = ',num2str(DAQ_time)])
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

artifact_time = toc; disp(['Artifact time = ',num2str(artifact_time)])
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
    forArgs.TimeStart(CH) = curTime(CH) - fltArgs.TimeShift(CH);
end
fltTails = fltData(3,:); fltAllData = fltData(4,:);

else
    fltTails = {}; fltAllData = {};
end

filter_time = toc; disp(['Filter time = ',num2str(filter_time)])
%% Forecast 
tic
forBuffedOut = cell(size(forBuffs));

if doFor

[forTails, forBuffsAdd, forArgs] = forFun(forArgs, ...
    [rawAllData(selRaw2For), fltAllData(selFlt2For)]);
lenFor = [lenRaw(selRaw2For), lenFlt(selFlt2For)];
rpt = ceil(width(forData)/width(lenFor));
lenFor = repmat(lenFor, [1,rpt]);
if ~(size(forTails,2) == size(forData,2))
    error('Forecast channels are inconsistent.');
end
for CH = 1:size(forData,2)
    [forData{2,CH}, forData{3,CH}, forData{4,CH}] = ...
        bufferjuggle(forData{2,CH},forData{3,CH},forTails{CH}, ...
        @(old, new) bufferDataOverwrite(old, new, lenFor(CH)));
    [~, forBuffs{1,CH}, forBuffRow(1,CH), forBuffedOut{1,CH}] = ...
        bufferStorage(forBuffs{1,CH}, forBuffRow(1,CH), forBuffsAdd{1,CH} + curTime(CH));
end

end

forecast_time = toc; disp(['Forecast time = ',num2str(forecast_time)])

%% helper 
    function [newBuffer, newTail, newAll] = ...
        bufferjuggle(oldBuffer, oldTail, newData, bufferFunc)
        newBuffer = bufferFunc(oldBuffer, oldTail); 
        newTail = newData; 
        newAll = bufferFunc(newBuffer, newTail); 
    end

end