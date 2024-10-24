function [dataReceived, svN, timeBuff, forBuff, ...
    tPltRng, rawPlt, fltPlt, forPlt, ...
    rawD1, rawD4, fltD1, fltD4, forD1, forD4] = ...
    pollDataQueue_PhaseDetect_v1(dataQueue, chInd, svname, svN, t0, pollTimeOut)
% 
% Poll dataQueue as sent by PhaseDetect function and interpret the results.
% 
% If no data received by pollTimeOut (default 1s), dataReceived will be 
% false and returns will be empty. 
% 
% If there is data to save, it will be saved as <svname_svN.mat> and svN
% will be incremented. 
% 
% chInd is the index (NOT ID number or name) of the raw channel being
% displayed, filtered, etc. 
% 

if nargin < 6
    pollTimeOut = 1; % s
end

if isempty(chInd)
    chInd = 1; % default to first listed channel
end

dataQueue.QueueLength

    [sentData, dataReceived] = poll(dataQueue, pollTimeOut);
    if dataReceived
        if strcmpi(class(sentData), 'MException')
            rethrow(sentData)
        else

        rawD1 = sentData{1,1}; rawD4 = sentData{2,1};
        fltD1 = sentData(1,2); fltD4 = sentData(2,2);
        forD1 = sentData(1,3); forD4 = sentData(2,3);
        timeBuffs = sentData{3,1}; forBuff = sentData{3,3}; 
        forBuffSv = sentData{3,2};

        timeBuff = timeBuffs{chInd};

        if ~isempty(forBuffSv)
            PeakTrough = forBuffSv;
            save([svname,num2str(svN),'.mat'], 'PeakTrough');
            svN = svN+1;
            forBuff = [forBuffSv; forBuff];
        end

        rawPlt = data2timetable(rawD4(chInd),rawD1(chInd),t0); rawPlt = rawPlt{1};
        fltPlt = data2timetable(fltD4,fltD1,t0); fltPlt = fltPlt{1};
        forPlt = data2timetable(forD4,forD1,t0); forPlt = forPlt{1};
        tPltRng = [gettimes(rawPlt); gettimes(fltPlt); gettimes(forPlt)];
        tPltRng = [min(tPltRng), max(tPltRng)];
        tPltRng = tPltRng + [-1,1]*.1*diff(tPltRng);

        end

    else
        tPltRng = []; 
        rawPlt = []; fltPlt = []; forPlt = [];
        timeBuff = []; forBuff = []; 
        rawD1 = []; rawD4 = []; 
        fltD1 = []; fltD4 = []; 
        forD1 = []; forD4 = [];
    end

    function t = gettimes(tbl)
        if numel(tbl)
            t = tbl.Time;
        else
            t = t0 + nan;
        end
    end

end