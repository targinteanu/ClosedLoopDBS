function [...
    timeBuffs, rawData, ...
    artRemData, artRemArgs, ...
    fltData, fltArgs, ...
    forBuffs, forData, forArgs, ...
    setup_time, DAQ_time, artifact_time, filter_time, update_time, forecast_time] = ...
    iterReadBrain(...
        timeBuffs, rawData, daqFun, ...
        selRaw2Art, selFor2Art, selRaw2Flt, selRaw2For, selFlt2For, ...
        artRemData, artRemFun, artRemArgs, ...
        fltData, fltFun, fltArgs, ...
        forBuffs, forData, forFun, forArgs)
% One iteration of brain reading, to be used as a main loop or timer
% function. 
% 
% Input Details: 
%   timeBuffs columns correspond to columns of rawData 
%   xxxData is a cell with 4 rows: (1) name, (2) head, (3) tail, (4) complete 
%   selxxx2xxx is a horizontal array of selected column indexes
%   calling daqFun() gets new raw tails 
%   stimRef is reference signal for new stimuli that would cause artifact
%   [artRemTails, artRemArgs] = artRemFun(artRemArgs, rawTails, forTails)
%   [fltTails, fltArgs] = fltFun(fltArgs, rawTails)
%       fltArgs must have an array fltArgs.TimeShift which has each
%       channel's time shift (in seconds) caused by the filter delay 
%   [forTails, forBuffsAdd, forArgs] = forFun(forArgs, inData)
%       forArgs must have an array forArgs.TimeStart which has each
%       channel's starting time (seconds) where the forecast begins 

%% handle/check inputs 
tic

doDAQ = (~isempty(rawData)) && (~isempty(daqFun));
doArt = doDAQ && (~isempty(artRemData)) && (~isempty(artRemFun));
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

setup_time = toc; %disp(['Setup Time = ',num2str(setup_time)])
%% DAQ 
tic
if doDAQ 

chInfo = rawData(1,:);
%rawNames = cellfun(@(c) c.Name, chInfo, 'UniformOutput',false);
rawIDs = cellfun(@(c) c.IDnumber, chInfo);

curTime = nan(1,size(rawData,2));

lenLastRaw = cellfun(@height, rawData(3,:)); lenLastFlt = lenLastRaw(selRaw2Flt);

[newTails, tailNames, tailIDs] = daqFun(); 
if width(rawData) == width(newTails)
    % assume order is the same for timing 
    for CH = 1:width(rawData)
        newTail = newTails{CH};
        [rawData{2,CH}, rawData{3,CH}, rawData{4,CH}] = ...
            bufferjuggle(rawData{2,CH},rawData{3,CH},newTail,@bufferData);
        fs = chInfo{CH}.SampleRate;
        tailProcTime = newTail(1,1) + (height(newTail)-1)/fs; 
        curTime(CH) = tailProcTime;
        timeBuffs{1,CH} = bufferData(timeBuffs{1,CH}, tailProcTime);
    end
else
% assign manually using channel ID
% should the code not get here, since width was set in daqFun() ?
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
end

% selection, etc
rawTails = rawData(3,:); rawAllData = rawData(4,:);
try
    rawTails = rawTails([selRaw2Art, selRaw2Flt, selRaw2For]); 
    rawAllData = rawAllData([selRaw2Art, selRaw2Flt, selRaw2For]);
catch ME
    if strcmp(ME.identifier, 'MATLAB:badsubscript')
        error('Requested selected channel(s) that do(es) not exist in raw data.');
    else
        rethrow(ME);
    end
end
selRaw2Art = 1:length(selRaw2Art);
selRaw2Flt = length(selRaw2Art) + (1:length(selRaw2Flt)); 
selRaw2For = length(selRaw2Art) + length(selRaw2Flt) + (1:length(selRaw2For));
lenRaw = cellfun(@height, rawTails); lenFlt = lenRaw(selRaw2Flt);

end

DAQ_time = toc; %disp(['DAQ time = ',num2str(DAQ_time)])
%% Artifact Removal 
tic
if doArt

nOverlapOld = artRemArgs.nOverlap;
[artRemTails, artRemArgs] = artRemFun(artRemArgs, ...
    rawTails(selRaw2Art), forData(3,selFor2Art));
if ~(size(artRemTails,2) == size(artRemData,2))
    error('Artifact removal channels are inconsistent.');
end
for CH = 1:size(artRemData, 2)
    oldHead = artRemData{2,CH}; oldTail = artRemData{3,CH}; newTail = artRemTails{CH};
    oldTime = nOverlapOld(CH); newTime = artRemArgs.nOverlap(CH); % overwrite times 
    newHead =       bufferDataOverwrite(oldHead, oldTail, oldTime); 
    artRemData{2,CH} = newHead;
    artRemData{3,CH} = newTail; 
    artRemData{4,CH} = bufferDataOverwrite(newHead, newTail, newTime);
end
% overwrite raw data - internally only - ???
rawTails(selRaw2Art) = artRemData(3,:); rawAllData(selRaw2Art) = artRemData(4,:);

end

artifact_time = toc; %disp(['Artifact time = ',num2str(artifact_time)])
%% Filter
tic
if doFlt
 
[fltTails, fltArgs] = fltFun(fltArgs, rawTails(selRaw2Flt));
if ~(size(fltTails,2) == size(fltData,2))
    error('Filtered channels are inconsistent.');
end
curTimeFlt = curTime(selRaw2Flt);
for CH = 1:size(fltData,2)
    [fltData{2,CH}, fltData{3,CH}, fltData{4,CH}] = ...
        bufferjuggle(fltData{2,CH},fltData{3,CH},fltTails{CH},@bufferData);
    CHfor = CH + length(selRaw2For); % corresponding forecast channel
    forArgs.TimeStart(CHfor) = curTimeFlt(CH) - fltArgs.TimeShift(CH);
end
fltTails = fltData(3,:); fltAllData = fltData(4,:);

else
    fltTails = {}; fltAllData = {};
end

filter_time = toc; %disp(['Filter time = ',num2str(filter_time)])

%% Update mdl coeffs 
tic 



update_time = toc; %disp(['Model update time = ',num2str(update_time)])

%% Forecast 
tic

if doFor

[forTails, forBuffsAdd, forArgs] = forFun(forArgs, ...
    [rawAllData(selRaw2For), fltAllData(selFlt2For)]);
lenFor = [lenLastRaw(selRaw2For), lenLastFlt(selFlt2For);
          lenRaw(selRaw2For),     lenFlt(selFlt2For)];
if ~(size(forTails,2) == size(forData,2))
    error('Forecast channels are inconsistent.');
end
curTimeFor = [curTime(selRaw2For), curTimeFlt(selFlt2For)];
for CH = 1:size(forData,2)
    oldHead = forData{2,CH}; oldTail = forData{3,CH}; newTail = forTails{CH};
    oldTime = lenFor(1,CH); newTime = lenFor(2,CH); % overwrite times 
    newHead =       bufferDataOverwrite(oldHead, oldTail, oldTime); 
    forData{2,CH} = newHead;
    forData{3,CH} = newTail; 
    forData{4,CH} = bufferDataOverwrite(newHead, newTail, newTime);
    forBuffCH = forBuffs{1,CH}; forBuffAddCH = forBuffsAdd{1,CH} + curTimeFor(CH);
    for p = 1:width(forBuffCH)
        if forBuffCH(end,p) > curTimeFor(CH) + forArgs.StimulatorLagTime
            % prev point is still in the accessible future; replace it
            forBuffCH(end,p) = forBuffAddCH(end,p);
        else
            % prev point passed; add new point to buffer 
            forBuffCH(:,p) = bufferData(forBuffCH(:,p), forBuffAddCH(end,p));
        end
    end
    forBuffs{1,CH} = forBuffCH;
end

end

forecast_time = toc; %disp(['Forecast time = ',num2str(forecast_time)])

%% helper 
    function [newBuffer, newTail, newAll] = ...
        bufferjuggle(oldBuffer, oldTail, newData, bufferFunc)
        newBuffer = bufferFunc(oldBuffer, oldTail); 
        newTail = newData; 
        newAll = bufferFunc(newBuffer, newTail); 
    end

end