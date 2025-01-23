function [dataReceived, svN, timeBuff, forBuff, stimBuff, ...
    tPltRng, rawPlt, fltPlt, forPlt, artPlt, ...
    rawD1, rawD4, fltD1, fltD4, forD1, forD4, artD1, artD4, ...
    forStore, forP, stimStore, stimP, srlStore, srlP, ...
    srlUserData, srlString] = ...
    pollDataQueue_PhaseDetect_v1(dataQueue, chInd, svname, svN, t0, pollTimeOut, ...
    forStore, forP, stimStore, stimP, srlStore, srlP)
% 
% Poll dataQueue as sent by PhaseDetect function and interpret the results.
% 
% If no data received by pollTimeOut (default 1s), dataReceived will be 
% false and returns will be empty. [This may no longer work!]
% 
% If needed to save, data will be saved as <svname_svN.mat> and svN
% will be incremented. 
% 
% chInd is the index (NOT ID number or name) of the raw channel being
% displayed, filtered, etc. 
% 

if nargin < 6
    pollTimeOut = 1; % s
end
if nargin < 7
    forStore = [];
end
if nargin < 9
    stimStore = [];
end
if nargin < 11
    srlStore = [];
end

if isempty(chInd)
    chInd = 1; % default to first listed channel
end

        tPltRng = []; 
        rawPlt = []; fltPlt = []; forPlt = []; artPlt = [];
        timeBuff = []; forBuff = []; stimBuff = [];
        rawD1 = []; rawD4 = []; 
        fltD1 = []; fltD4 = []; 
        forD1 = []; forD4 = [];
        artD1 = []; artD4 = [];
        srlUserData = []; srlString = '';
        dataReceived = false;

[sentData, dataReceivedNow] = poll(dataQueue, pollTimeOut); 

dopoll = true;
while dopoll
    % poll until Q is empty to get most recent data
    dopoll = dataQueue.QueueLength > 0;
    dataReceived = dataReceived || dataReceivedNow;

    if dataReceivedNow
        if strcmpi(class(sentData), 'MException')
            rethrow(sentData)
        end

        forBuff = sentData{3,3}; stimBuff = sentData{3,2};
        srlBuff = sentData{4,1}; srlUserData = sentData{4,2}; srlString = sentData{4,3};

        if ~isempty(srlStore)
            srlTimeStamp = [srlBuff.TimeStamp];
            srlSel = ~isnan(srlTimeStamp); 
            srlBuff = srlBuff(srlSel); 
            [srlFull, srlStore, srlP, srlBuffSv] = bufferStorage(...
                srlStore, srlP, srlBuff);
            if srlFull
                SerialLog = srlBuffSv; 
                save([svname,num2str(svN),'.mat'], 'SerialLog');
                svN = svN+1;
            end
        end
        
        if ~isempty(stimStore)
            stimBuffAll = [stimStore; stimBuff];
            [~,stimU] = unique(stimBuffAll, 'stable');
            stimU = stimU( stimU > height(stimStore) );
            stimBuffNew = stimBuffAll(stimU, :);
            [stimFull, stimStore, stimP, stimBuffSv] = bufferStorage(...
                stimStore, stimP, removenan(stimBuffNew) );
            if stimFull
                Stim = stimBuffSv;
                save([svname,num2str(svN),'.mat'], 'Stim');
                svN = svN+1;
            end
        end

        if ~isempty(forStore)
            forBuffAll = [forStore; forBuff];
            [~,forU] = unique(forBuffAll, 'stable', 'rows');
            forU = forU( forU > height(forStore) );
            forBuffNew = forBuffAll(forU, :);
            [forFull, forStore, forP, forBuffSv] = bufferStorage(...
                forStore, forP, removenan(forBuffNew) );
            if forFull
                PeakTrough = forBuffSv;
                save([svname,num2str(svN),'.mat'], 'PeakTrough');
                svN = svN+1;
            end
        end

        rawD1 = sentData{1,1}; rawD4 = sentData{2,1};
        fltD1 = sentData(1,2); fltD4 = sentData(2,2);
        forD1 = sentData(1,3); forD4 = sentData(2,3);
        artD1 = sentData(1,4); artD4 = sentData(2,4);
        timeBuffs = sentData{3,1}; 
        timeBuff = timeBuffs{chInd};
    end

    sentData = poll(dataQueue);
    dataReceivedNow = ~isempty(sentData);
end

    if dataReceived
        rawPlt = data2timetable(rawD4(chInd),rawD1(chInd),t0); rawPlt = rawPlt{1};
        fltPlt = data2timetable(fltD4,fltD1,t0); fltPlt = fltPlt{1};
        forPlt = data2timetable(forD4,forD1,t0); forPlt = forPlt{1};
        artPlt = data2timetable(artD4,artD1,t0); artPlt = artPlt{1};
        %tPltRng = [gettimes(rawPlt); gettimes(fltPlt); gettimes(forPlt)];
        tPltRng = gettimes(rawPlt);
        tPltRng = [min(tPltRng), max(tPltRng)];
        tPltRng = tPltRng + [-1,1]*.1*diff(tPltRng);
    end

    function t = gettimes(tbl)
        if numel(tbl)
            t = tbl.Time;
        else
            t = t0 + nan;
        end
    end

    function y = removenan(x)
        b = true(height(x),1);
        for c = 1:width(x)
            b = b & ~isnan(x(:,c));
        end
        y = x(b,:);
    end

end