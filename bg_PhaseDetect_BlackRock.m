function bgArgOut = bg_PhaseDetect_BlackRock(UserArgs, DQ, StimController)
% 
% Run brain recording with phase detection/prediction for PDS and
% stimulation using BlackRock hardware (cbmex for recording, cerestim for
% stimulation). 
% 
% UserArgs can be e.g. a handles object from GUIDE or an app object. 
% 
% DQ = data and queue (respectively), i.e. DQ = parallel.pool.(Pollable)DataQueue  
% 
% StimController takes in UserArgs, Past Data, and Time to Peak/Trough 
% and outputs time to next stimulus. Meant so that a different function can
% be passed for the motor protocol vs memory protocol. 
% 

bgArgOut = [];

looptime = .01; % starting estimate loop time (s)
guitime = .4; % estimate of gui update time (s)

dT = .001; % s between data requests 

try
%% import filter and model details, etc from user

if ~UserArgs.DAQstatus
    error('Data Acquisition has not been enabled.')
end

cont_loop_2 = UserArgs.DAQstatus && UserArgs.RunMainLoop; 
    % if false, loop should only run once

FilterSetUp = UserArgs.FilterSetUp; % t/f
MdlSetUp = FilterSetUp && UserArgs.MdlSetUp; % t/f
if FilterSetUp
    filtOrds = [UserArgs.FilterOrder]; % array with chans as cols
    filtA = 1; filtB = UserArgs.BPF;
    hico = UserArgs.hicutoff; loco = UserArgs.locutoff; % Hz 
else
    filtOrds = [];
end
if MdlSetUp
    mdls = UserArgs.Mdl;
end
rawIDs = UserArgs.allChannelIDs;
chInd = UserArgs.channelIndex; 
    % index (NOT ID NUMBER) of the recording channel/column; 
    % may be empty if unselected on startup
    if isempty(chInd)
        chInd = 1; % default to first listed channel
    end
forecastwin = UserArgs.PDSwin1; % # samples ahead to forecast
forecastpad = UserArgs.PDSwin2; % # of above to use to pad hilbert transform
buffSize = UserArgs.bufferSizeGrid;
PhaseOfInterest = UserArgs.PhaseOfInterest;

doStim = ((~isempty(StimController)) && UserArgs.StimActive) && (FilterSetUp && MdlSetUp);

% stimulator 
if UserArgs.StimActive
    stimulator = stimSetup_cerestim(UserArgs);
end

%% request raw data and get details  

% access cbmex 
connect_cbmex(); 
pause(1);
[spikeEvents, time, continuousData] = cbmex('trialdata',1);
if isempty(continuousData)
    error('No continuous data; ensure that data acquisition has been enabled.')
end

% setup 
chnum = [continuousData{:,1}]';
Fs = [continuousData{:,2}]';
chname = spikeEvents(:,1);
chname = chname(chnum);
if isempty(rawIDs)
    % select all channels
    rawIDs = chnum;
end

% initialize empty (for now)
emptyData = cell(1,length(rawIDs)); 
rawH = emptyData; % "head" 
rawT = emptyData; % "tail" 
rawB = emptyData; % "buffer" 
chanInfo = emptyData;

% handle/check inputs 
if length(buffSize) < length(rawIDs)
    if length(buffSize) == 1
        % assume the one input applies to all channels.
        buffSize = repmat(buffSize, size(rawIDs));
        buffSize(chInd) = UserArgs.bufferSize;
    else
        error('Incompatible input dimensions.')
    end
end


% assign raw data to the structure 

for ch = 1:length(rawIDs)
    indCH = find(chnum == rawIDs(ch)); 

    if ~isempty(indCH)

        % Create raw data buffer of zeros of the correct length
        L = length(continuousData{indCH,3}); 
        rawH{ch} = [nan(buffSize(ch),1), zeros(buffSize(ch),1)];
        rawH{ch}(end,1) = time - 1/Fs(indCH); 
        rawT{ch} = [nan(L,1), continuousData{indCH,3}];
        rawT{ch}(1,1) = time;

        % check units 
        config = cbmex('config', indCH);
        unitname_is = lower(config{11,1});
        if contains(unitname_is, 'unit')
            unitname = config{11,2};
        else
            unitname_is = contains(lower(config(:,1)), 'unit');
            unitname_is = find(unitname_is);
            unitname_is = unitname_is(1); 
            unitname = config{unitname_is,2};
        end

        % Channel Info
        ud.SampleRate = Fs(indCH);
        ud.Name = chname{indCH};
        ud.Unit = unitname;
        ud.IDnumber = chnum(indCH);
        chanInfo{ch} = ud;

        rawB{ch} = bufferData(rawH{ch}, rawT{ch});

    end
    
end


rawD = [chanInfo; rawH; rawT; rawB];
chIDs = cellfun(@(s) s.IDnumber, chanInfo); chID = chIDs(chInd);

rawInds = nan(size(rawIDs));
for ch = 1:length(rawInds)
    rawInd = find(rawIDs == chanInfo{ch}.IDnumber);
    if ~isempty(rawInd)
        rawInds(ch) = rawInd;
    end
end
if sum(isnan(rawInds))
    warning('Some of the requested IDs were not found.')
    rawInds = rawInds(~isnan(rawInds));
end

%% initialize processed data structs 

IndShiftFIR = ceil(filtOrds/2); % samples
emptyOut = {[]; []; []; []};

% filtered data 
if FilterSetUp && ( numel(IndShiftFIR) && numel(chInd) )
    fltH = initFilteredData(rawH(chInd), IndShiftFIR); 
    fltT = cellfun(@(D) [1,0].*D, rawT(chInd), 'UniformOutput',false);
    fltB = cellfun(@(D) [1,0].*D, rawB(chInd), 'UniformOutput',false);
    fltInfo = chanInfo(chInd);
    fltD = [fltInfo; fltH; fltT; fltB];
else
    fltD = emptyOut;
end

% AR-forecast data 
if MdlSetUp && numel(forecastwin)
    [forH, forT, forB] = initForecastData(fltH(1), forecastwin);
    forD = [fltInfo(1); forH; forT; forB];
else
    forD = emptyOut;
end

% timing buffers 
timeBuffs = cell(size(rawH));
for ch = 1:width(rawH)
    t_ch = rawH{ch}(:,1);
    dt_ch = 1/chanInfo{ch}.SampleRate;
    for it = 2:length(t_ch)
        if isnan(t_ch(it))
            t_ch(it) = t_ch(it-1) + dt_ch;
        end
    end
    timeBuffs{ch} = t_ch;
end

% init forecast-output buffers, i.e. times to/of next phase(s) of interest
forBuffs = cellfun(@(X) (nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);
forBuffs = forBuffs(chInd); 
stimBuff = nan(size(timeBuffs{chInd},1),1); 

%% define input args for filtering/forecasting  
Fs = cellfun(@(s) s.SampleRate, chanInfo); fs = Fs(chInd);
if FilterSetUp
    filtInit = arrayfun(@(ord) zeros(ord,1), filtOrds, 'UniformOutput',false);
    TimeShift = IndShiftFIR(1)/fs; 
end

%% loop
% main loop 

cont_loop = true; first_loop = true; looptime_meas1 = tic; loopcount = 0;
stimLastTime = -inf; stimtime = -inf;
while cont_loop

    % timing 
    looptime_meas2 = toc(looptime_meas1); 
    looptime_meas1 = tic; 
    if ~first_loop
        looptime = .9*looptime + .1*looptime_meas2;
    end
    loopcount = loopcount+1;
    loopsendnum = guitime/looptime;
    pause(dT)

    try
%% data aquisition 

curTime = nan(1,size(rawD,2));

[newTails, tailNames, tailIDs] = getNewRawData_cbmex(rawIDs); 
for ch = 1:size(newTails,2)
    newTail = newTails{ch};
    tailName = tailNames{ch};
    tailID = tailIDs(ch);
    CH = find(tailID == rawIDs);
    if isempty(CH)
        error(['Unrecognized new raw data label: ',tailName]);
    end
    if length(CH) > 1
        error(['Raw data label ',tailName,' is not unique.']);
    end

    fs_ch = chanInfo{CH}.SampleRate;
    tailProcTime = newTail(1,1) + (height(newTail)-1)/fs_ch; 
    curTime(CH) = tailProcTime;

    [rawD{2,CH}, rawD{3,CH}, rawD{4,CH}] = ...
        bufferjuggle(rawD{2,CH},rawD{3,CH},newTail,@bufferData);
    timeBuffs{1,CH} = bufferData(timeBuffs{1,CH}, tailProcTime);
end

% selection, etc
rawTails = rawD(3,:); rawAllData = rawD(4,:);
lenRaw = cellfun(@height, rawTails); lenFlt = lenRaw(chInd);
try
    rawTails = rawTails(chInd); 
    rawAllData = rawAllData(chInd);
catch ME
    if strcmp(ME.identifier, 'MATLAB:badsubscript')
        error('Requested selected channel(s) that do(es) not exist in raw data.');
    else
        rethrow(ME);
    end
end

%% Filter
if FilterSetUp

fltTails = fltD(3,:);
for CH = 1:size(rawTails,2)
    [fltTail,filtInit{CH}] = filter(filtB,filtA,rawTails{CH}(:,2:end),filtInit{CH});
    fltTails{CH} = [rawTails{CH}(:,1) - TimeShift(CH), fltTail];
end

if ~(size(fltTails,2) == size(fltD,2))
    error('Filtered channels are inconsistent.');
end
curTimeFlt = curTime(chInd);
TimeStart = nan(1,width(fltD));
for CH = 1:size(fltD,2)
    [fltD{2,CH}, fltD{3,CH}, fltD{4,CH}] = ...
        bufferjuggle(fltD{2,CH},fltD{3,CH},fltTails{CH},@bufferData);
    TimeStart(CH) = curTimeFlt(CH) - TimeShift(CH);
end
fltAllData = fltD(4,:);

else
    fltAllData = {};
end

%% Forecast 
if MdlSetUp

forTails = cell(size(fltAllData)); forBuffsAdd = forTails; 
for CH = 1:size(fltAllData,2)
    fltTail = myFastForecastAR(mdls, fltAllData{CH}(:,2), forecastwin);
    forTails{CH} = [nan(height(fltTail),1), fltTail];
    forTails{CH}(1,1) = TimeStart(CH);

    fltTail = fltTail(1:forecastpad,:); % use limited duration for hilbert padding
    [t2,i2,phi_inst,f_inst] = blockPDS(...
        fltAllData{CH}(:,2), fltTail, fs(CH), PhaseOfInterest, ...
        TimeShift(CH), loco, hico);
    t2 = t2-TimeShift(CH); % [t2peak, t2trough]
    t2 = max(t2,0);
    forBuffsAdd{CH} = t2; 
end

lenFor = lenFlt(1);
if ~(size(forTails,2) == size(forD,2))
    error('Forecast channels are inconsistent.');
end
curTimeFor = curTimeFlt(1);
for CH = 1:size(forD,2)
    [forD{2,CH}, forD{3,CH}, forD{4,CH}] = ...
        bufferjuggle(forD{2,CH},forD{3,CH},forTails{CH}, ...
        @(old, new) bufferDataOverwrite(old, new, lenFor(CH)));
    forBuffCH = forBuffs{1,CH}; forBuffAddCH = forBuffsAdd{1,CH} + curTimeFor(CH);
    for p = 1:width(forBuffCH)
        if forBuffCH(end,p) > curTimeFor(CH)
            % prev point is still in the future; replace it
            forBuffCH(end,p) = forBuffAddCH(end,p);
        else
            % prev point passed; add new point to buffer 
            forBuffCH(:,p) = bufferData(forBuffCH(:,p), forBuffAddCH(end,p));
        end
    end
    forBuffs{1,CH} = forBuffCH;
end

end

%% do stimulus 
forBuff = forBuffs{1}; 
forBuffNew = [max(forBuff(:,1)), max(forBuff(:,2))]; 
forBuffNew = forBuffNew - timeBuffs{chInd}(end,:); % [t2p, t2t]
if doStim
    [t2stim, stim2q] = StimController(UserArgs, fltD{4,1}(:,2), forBuffNew);
    if stim2q && (t2stim < 2*looptime)
        t2stim = .001*floor(1000*t2stim); % round to nearest 1ms 
        % ensure below max freq
        if 1/(t2stim + timeBuffs{chInd}(end,:) - stimLastTime) <= UserArgs.stimMaxFreq
            % check if stimulator is ready here?
            % assume stim should occur before completion of next loop
            pause(t2stim);
            stimtime1 = cbmex('time');
            stimulator.play(1);
            stimtime2 = cbmex('time');
            stimtime = .5*(stimtime1 + stimtime2);
            stimBuff = bufferData(stimBuff, stimtime);
            stimLastTime = stimtime;
        end
    end
end

%% send data
% User should be ready for new data when loopcount = loopsendum
% send data AT LEAST that frequently so the user never waits for data
if first_loop || (loopcount >= .5*loopsendnum)
    loopcount = 0;
    if isempty(rawIDs)
        rawD_ = rawD;
        timeBuffs_ = timeBuffs;
    else
        rawD_ = rawD(:,rawInds);
        timeBuffs_ = timeBuffs(rawInds);
    end
    % This relies on the forecast buffers being 2xN and the stim buffer
    % being 1xN; should be made more robust. 
    send(DQ, [{rawD_(1,:)}, fltD(1,1), forD(1,1); ...
              {rawD_(4,:)}, fltD(4,1), forD(4,1); ...
              {timeBuffs_}, {stimBuff}, {forBuff}]);
end
% bgArgOut = [forBuff, stimBuff];

    %% end loop logic 
    cont_loop = cont_loop && cont_loop_2;
    first_loop = false;

    % If there are any errors in the loop, stop looping and allow the
    % User to handle them. 
    catch ME_loop
        getReport(ME_loop)
        if contains(ME_loop.message, 'No continuous data')
            % do nothing; proceed to next loop iteration to allow more time
            % for continuous data. 
            % TO DO: should there be some limit; enough of these in a row
            % and it stops cont_loop and sends the error? 
        else
            cont_loop = false;
            send(DQ, ME_loop);
        end
    end
end

%% stop 
cbmex('close'); 

% stimulator 
if UserArgs.StimActive
    stimulator.stop();
    stimulator.disconnect();
end

% If there are any errors in the full func, stop looping and allow the User
% to handle them. 
catch ME_fullfunc
    getReport(ME_fullfunc)
    send(DQ, ME_fullfunc)
end

%% helper function def 

    function [newBuffer, newTail, newAll] = ...
        bufferjuggle(oldBuffer, oldTail, newData, bufferFunc)
        newBuffer = bufferFunc(oldBuffer, oldTail); 
        newTail = newData; 
        newAll = bufferFunc(newBuffer, newTail); 
    end

end