function bgArgOut = bg_PhaseDetect(UserArgs, DQ, SQ, ...
    SetupRecording, ShutdownRecording, ...
    SetupStimulator, ShutdownStimulator, PulseStimulator, ...
    GetNewRawData, GetTime)
% 
% Run brain recording with phase detection/prediction for PDS.
% 
% UserArgs can be e.g. a handles object from GUIDE or an app object. 
% 
% DQ and SQ are data and stimulus queues (respectively), i.e. 
% Q = parallel.pool.(Pollable)DataQueue.  
% 
% InitializeRecording takes arguments:
% ( buffer size(s) , filter order(s) , forecast window(s) , 
%   raw channel ID(s) selected , raw channel ID(s) to filter , 
%   raw channel ID(s) to use for forecasting , 
%   filtered channel ID(s) to use for forecasting ) 
% and returns data structures for: 
% [ raw , filtered , forecast, timing buffer ] data and starting tic
% 
% ShutdownRecording takes no arguments.
% 
% SetupStimulator takes UserArgs and initiates/returns StimArgs;
% ShutdownStimulator and PulseStimulator take StimArgs and UserArgs and
% return StimArgs. 
% 
% InitializeRawData takes ( selRaw , BufferSize ) as input arguments and
% returns [ EmptyDataStruct , Tail , Combined , ChannelInfo, StartTic ]
% 
% GetNewRawData takes selRaw as an argument and returns new raw tails. 
%
% GetTime takes Start Tic and returns elapsed time (s); intended to be
% machine time of recording device. 
% 

bgArgOut = [];

looptime = .01; % starting estimate loop time (s)
guitime = .4; % estimate of gui update time (s)

try
%% import filter and model details, etc

if ~UserArgs.DAQstatus
    error('Data Acquisition has not been enabled.')
end

cont_loop_2 = UserArgs.DAQstatus && UserArgs.RunMainLoop; 
    % if false, loop should only run once

StimController = UserArgs.ControllerResult;
StimulatorLagTime = UserArgs.StimulatorLagTime; % s
doArtRem = UserArgs.check_artifact_Value;
FilterSetUp = UserArgs.FilterSetUp; % t/f
MdlSetUp = FilterSetUp && UserArgs.MdlSetUp; % t/f
if FilterSetUp
    filtOrds = [UserArgs.FilterOrder]; % array with chans as cols
    filtObjs = {1; UserArgs.BPF}; % cell with rows {a; b}; chans as cols
    hico = UserArgs.hicutoff; loco = UserArgs.locutoff; % Hz 
else
    filtOrds = [];
end
if MdlSetUp
    mdls = {UserArgs.Mdl};
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
PhaseOfInterest = UserArgs.PhaseOfInterest;

rawN = UserArgs.allChannelInfo;
chID = cellfun(@(s) s.IDnumber, rawN); chID = chID(chInd);

selRaw2Flt = UserArgs.selInds.selRaw2Flt;
selFlt2For = UserArgs.selInds.selFlt2For;
selFor2Art = UserArgs.selInds.selFor2Art;
selRaw2Art = UserArgs.selInds.selRaw2Art;
selRaw2For = UserArgs.selInds.selRaw2For;

rawD = UserArgs.recDataStructs.rawD;
artD = UserArgs.recDataStructs.artD;
fltD = UserArgs.recDataStructs.fltD;
forD = UserArgs.recDataStructs.forD;
timeBuffs = UserArgs.recDataStructs.timeBuffs;
initTic = UserArgs.recDataStructs.initTic;
forBuffs = UserArgs.recDataStructs.forBuffs;
stimBuff = UserArgs.recDataStructs.stimBuff;

%% init 

dT = .001; % s between data requests 

SetupRecording();

rawInds = nan(size(rawIDs));
for ch = 1:length(rawInds)
    rawInd = find(rawIDs == rawN{ch}.IDnumber);
    if ~isempty(rawInd)
        rawInds(ch) = rawInd;
    end
end
if sum(isnan(rawInds))
    warning('Some of the requested IDs were not found.')
    rawInds = rawInds(~isnan(rawInds));
end

bgArgOut = UserArgs;

% define input args for filtering/forecasting funcs 
Fs = cellfun(@(s) s.SampleRate, rawN); 
if FilterSetUp
    fltN = fltD(1,:); 
    fIC = arrayfun(@(ord) zeros(ord,1), filtOrds, 'UniformOutput',false);
    filtArgs.fltInit = fIC; filtArgs.fltObj = filtObjs;
    filtArgs.TimeShift = filtOrds(1)/(2*Fs(selRaw2Flt)); 
    if MdlSetUp
        forN = forD(1,:);
        foreArgs.K = forecastwin; foreArgs.k = forecastpad;
        foreArgs.TimeStart = nan(size(forN));
        foreArgs.TimeShift = [zeros(size(selRaw2For)), filtArgs.TimeShift(selFlt2For)];
        foreArgs.ARmdls = mdls;
        foreArgs.SampleRates = Fs(selRaw2Flt); % !!! needs improvement
        foreArgs.FreqRange = [loco, hico];
        foreArgs.PhaseOfInterest = PhaseOfInterest;
        foreArgs.StimulatorLagTime = StimulatorLagTime;
        artRemArgs.SampleRates = Fs(selRaw2Art);
        artRemArgs.StimDur = .11; % seconds
        artRemArgs.StimTimes = cell(size(artRemArgs.SampleRates));
        artRemArgs.nOverlap = zeros(size(artRemArgs.SampleRates));
    else
        foreArgs = [];
        artRemArgs = [];
    end
else
    filtArgs = [];
    foreArgs = [];
    artRemArgs = [];
end

% stimulator 
if UserArgs.StimActive
    if ~UserArgs.StimTriggerMode
        StimArgs = SetupStimulator(UserArgs.StimSetupArgs);
    else
        StimArgs = [];
    end
    stimScheduler = timer(...
        "TimerFcn",{@stimFunction, StimArgs, UserArgs, initTic}, ...
        "StartDelay",inf);
end

% enable defined filtering/forecasting funcs only if ready
if FilterSetUp
    filtfun = @filtFun;
else
    filtfun = [];
end
if MdlSetUp
    forefun = @foreFun;
else
    forefun = [];
end
if doArtRem && MdlSetUp
    artremfun = @artRemFun;
else
    artremfun = [];
end
doStim = ((~isempty(StimController)) && UserArgs.StimActive) && (FilterSetUp && MdlSetUp);

%% loop 
cont_loop = true; first_loop = true; looptime_meas1 = tic; loopcount = 0;
stimLastTime = -inf; stimtime = -inf;
while cont_loop
    looptime_meas2 = toc(looptime_meas1); 
    looptime_meas1 = tic; 
    if ~first_loop
        looptime = .9*looptime + .1*looptime_meas2;
    end
    loopcount = loopcount+1;
    loopsendnum = guitime/looptime;
    pause(dT)

    try
    % main iteration 
    [...
    timeBuffs, rawD, ...
    artD, artRemArgs, ...
    fltD, filtArgs, ...
    forBuffs, forD, foreArgs] = ...
    iterReadBrain(...
        timeBuffs, rawD, @() GetNewRawData(rawIDs), ...
        selRaw2Art, selFor2Art, selRaw2Flt, selRaw2For, selFlt2For, ...
        artD, artremfun, artRemArgs, ...
        fltD, filtfun, filtArgs, ...
        forBuffs, forD, forefun, foreArgs);

    % do stimulus 
    forBuff = forBuffs{1}; 
    forBuffNew = [max(forBuff(:,1)), max(forBuff(:,2))]; 
    forBuffNew = forBuffNew - timeBuffs{chInd}(end,:); % [t2p, t2t]
    if doStim && (StimController > 0)
        t2stim = forBuffNew(:,StimController);
        t2stim = t2stim - StimulatorLagTime; % account for hardware delays
        if (t2stim >= 0) % should always be true by construction
            t2stim = .001*floor(1000*t2stim); % round to nearest 1ms 
            if t2stim >= 0 % should the minimum be set any higher?
            % ensure below max freq
            if 1/(t2stim + timeBuffs{chInd}(end,:) - stimLastTime) <= UserArgs.stimMaxFreq
                % check if stimulator is ready here?
                stop(stimScheduler); % if timer hasn't fired, overwrite it
                stimScheduler.StartDelay = t2stim; 
                start(stimScheduler);
            end
            end
        end
        stimtime = stimScheduler.UserData; stimScheduler.UserData = [];
        if numel(stimtime)
            stimBuff = bufferData(stimBuff, stimtime);
            if doArtRem
                artRemArgs.StimTimes = cellfun(...
                    @(T) stimBuff - T(end,:), timeBuffs(selRaw2Art), ...
                    'UniformOutput',false);
            end
            stimLastTime = stimtime(end);
        end
    end

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
        send(DQ, [{rawD_(1,:)}, fltD(1,1), forD(1,1), artD(1,1); ...
                  {rawD_(4,:)}, fltD(4,1), forD(4,1), artD(4,1); ...
                  {timeBuffs_}, {stimBuff}, {forBuff}, {[]} ]);
    end
    % bgArgOut = [forBuff, stimBuff];

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
ShutdownRecording();

% stimulator 
if UserArgs.StimActive && ~UserArgs.StimTriggerMode
    StimArgs = ShutdownStimulator(StimArgs, UserArgs);
end

% If there are any errors in the full func, stop looping and allow the User
% to handle them. 
catch ME_fullfunc
    getReport(ME_fullfunc)
    send(DQ, ME_fullfunc)
end
%end

%% function def 

function [artRemTails, artRemArgs] = ...
        artRemFun(artRemArgs, rawTails, forTails)
% forTails must be the forecast data starting at the same time as rawTails
artRemTails = cell(size(rawTails));
FsArt = artRemArgs.SampleRates;
StimDur = artRemArgs.StimDur; % seconds to remove
StimLen = ceil(StimDur.*FsArt); % #samples to remove 
StimTimesTail = artRemArgs.StimTimes; 
for ch_art = 1:size(rawTails,2)
    tXfor = forTails{ch_art};
    tX = rawTails{ch_art}; % [time, data]
    stimtimes = StimTimesTail{ch_art}; % time to stim (sec)
    stiminds = round(stimtimes * FsArt(ch_art));
    stiminds = stiminds(stiminds > 0);
    for i1 = stiminds'
        i2 = i1 + StimLen(ch_art);
        if i2 > height(tX)
            % artifact will carry over into the next packet 
            nO = i2 - height(tX);
            artRemArgs.nOverlap(ch_art) = nO;
            tX = [tX; nan(nO,2)]; 
        else
            % artifact limited to this packet 
            artRemArgs.nOverlap(ch_art) = 0;
        end
        tX(i1:i2,2) = tXfor(i1:i2,2);
    end
    artRemTails{ch_art} = tX;
end
end

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the (filtered) channel(s) of interest
% foreTails will be the forecast tails starting at the current time minus
% any filtering delay
if isempty(foreArgs)
    error('Attempted forecast function with empty arguments.')
end
K = foreArgs.K; k = foreArgs.k;
ARmdls = foreArgs.ARmdls;
fs = foreArgs.SampleRates;
Ts = foreArgs.TimeShift;
Fco = foreArgs.FreqRange;
phis = foreArgs.PhaseOfInterest;
TstimLag = foreArgs.StimulatorLagTime;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch_fore = 1:size(inData,2)
    armdl = ARmdls{ch_fore};
    FT = myFastForecastAR(armdl, inData{ch_fore}(:,2), K);
    foreTails{ch_fore} = [nan(height(FT),1), FT];
    foreTails{ch_fore}(1,1) = foreArgs.TimeStart(ch_fore);

    FT = FT(1:k,:); % use limited duration for hilbert padding
    %[t2,i2,phi_inst,f_inst] = blockPDS(...
    t2 = blockPDS(...
        inData{ch_fore}(:,2), FT, fs(ch_fore), phis, ...
        Ts(ch_fore)+TstimLag, Fco(1), Fco(2));
    t2 = t2-Ts(ch_fore); % [t2peak, t2trough]
    % t2 = max(t2,0); % Should this be necessary? Will this cause problems?
    t2(t2 < 0) = nan;
    foreBuffsAdd{ch_fore} = t2; 
end
end

function [fltTails, fltArgs] = filtFun(fltArgs, rawTails)
filtObj = fltArgs.fltObj;
filtInit = fltArgs.fltInit;
filtFin = cell(size(filtInit));
fltTails = cell(size(rawTails));
for ch_flt = 1:size(rawTails,2)
    a = filtObj{1,ch_flt}; b = filtObj{2,ch_flt};
    [FT,filtFin{ch_flt}] = filter(b,a,rawTails{ch_flt}(:,2:end),filtInit{ch_flt});
    fltTails{ch_flt} = [rawTails{ch_flt}(:,1) - fltArgs.TimeShift(ch_flt), FT];
end
fltArgs.fltInit = filtFin;
end

    function stimFunction(timerObj, evt, StimArgs, UserArgs, initTic)
        stimtime1 = GetTime(initTic);
        StimArgs = PulseStimulator(StimArgs, UserArgs);
        stimtime2 = GetTime(initTic);
        stimt = .5*(stimtime1 + stimtime2);
        timerObj.UserData = [timerObj.UserData; stimt];
    end

end