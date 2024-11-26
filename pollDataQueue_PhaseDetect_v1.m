function [dataReceived, svN, timeBuff, forBuff, stimBuff, ...
    tPltRng, rawPlt, fltPlt, forPlt, ...
    rawD1, rawD4, fltD1, fltD4, forD1, forD4] = ...
    pollDataQueue_PhaseDetect_v1(dataQueue, chInd, svname, svN, t0, pollTimeOut)
% 
% Poll dataQueue as sent by PhaseDetect function and interpret the results.
% 
% If no data received by pollTimeOut (default 1s), dataReceived will be 
% false and returns will be empty. [This may no longer work!]
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

        tPltRng = []; 
        rawPlt = []; fltPlt = []; forPlt = [];
        timeBuff = []; forBuff = []; stimBuff = [];
        rawD1 = []; rawD4 = []; 
        fltD1 = []; fltD4 = []; 
        forD1 = []; forD4 = [];
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

        % This relies on the forecast buffers being 2xN and the stim buffer
        % being 1xN; should be made more robust. 
        ForStim = sentData{3,3}; 
        ForStimSv = sentData{3,2};
        if width(ForStim) > 2
            stimBuff = ForStim(:,3:end); forBuff = ForStim(:,1:2);
        else
            stimBuff = []; forBuff = ForStim;
        end
        if isempty(ForStimSv)
            stimBuffSv = ForStimSv; forBuffSv = ForStimSv; 
        elseif width(ForStimSv) < 2
            stimBuffSv = ForStimSv; forBuffSv = [];
        elseif width(ForStimSv) > 2
            stimBuffSv = ForStimSv(:,3:end); forBuffSv = ForStimSv(:,1:2);
        else
            stimBuffSv = []; forBuffSv = ForStimSv;
        end

        if ~isempty(forBuffSv)
            PeakTrough = forBuffSv;
            save([svname,num2str(svN),'.mat'], 'PeakTrough');
            svN = svN+1;
            forBuff = [forBuffSv; forBuff];
        end
        if ~isempty(stimBuffSv)
            Stim = stimBuffSv;
            save([svname,num2str(svN),'.mat'], 'Stim');
            svN = svN+1;
            stimBuff = [stimBuffSv; stimBuff];
        end

        rawD1 = sentData{1,1}; rawD4 = sentData{2,1};
        fltD1 = sentData(1,2); fltD4 = sentData(2,2);
        forD1 = sentData(1,3); forD4 = sentData(2,3);
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
        tPltRng = [gettimes(rawPlt); gettimes(fltPlt); gettimes(forPlt)];
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

end