function bgArgOut = bg_PhaseDetect(UserArgs, DQ, SQ, ...
    InitializeRecording, ShutdownRecording, ...
    SetupStimulator, ShutdownStimulator, PulseStimulator, ...
    InitializeRawData, GetNewRawData, GetTime, StimController)
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
% StimController takes in UserArgs, Past Data, and Time to Peak/Trough 
% and outputs time to next stimulus. 
% 

SRLA = serialportlist("available");
w = whos;
ww = [w.name,' ',w.class];
f = fields(UserArgs);
ff = [f, repmat({' '},size(f))]';
ff = [ff{:}];

bgArgOut = [];

looptime = .01; % starting estimate loop time (s)
guitime = .4; % estimate of gui update time (s)

cont_fullfunc = true; % run or wait for user input
%while cont_fullfunc
try
%% import filter and model details, etc

if ~UserArgs.DAQstatus
    error('Data Acquisition has not been enabled.')
end

cont_loop_2 = UserArgs.DAQstatus && UserArgs.RunMainLoop; 
    % if false, loop should only run once

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
buffSize = UserArgs.bufferSizeGrid;
PhaseOfInterest = UserArgs.PhaseOfInterest;

%% setup serial comm

SerialArgs = UserArgs.SerialArgs;
noSerialSetup = SerialArgs.NoSerial;
srlUD = SerialArgs.UserData;
%{
if UserArgs.srlHere
    % serial is already set up on user thread! 
    error('Attempt to start serial on multiple threads.')
end
srlCBFn = SerialArgs.CallbackFcn;
%}

%{
if ~sum(strcmp(SRLA, SerialArgs.PortName))
    error(['Port ',char(SerialArgs.PortName),' was not found. ',...
        'UserArgs fields are: ',ff,' ; ',...
        'workspace variables are: ',ww])
end
%}

%{
if ~noSerialSetup
    try
    receiverSerial = serialport(SerialArgs.PortName, 9600);
    configureCallback(receiverSerial,"terminator", ...
        @(hsrl,evt)srlCBFn(hsrl,evt)); % no GUI object args passed - must handle in DQ
    catch MEsrl
        if strcmpi(MEsrl.identifier, 'serialport:serialport:ConnectionFailed')
            srlavailstr = squeeze(char(SRLA));
            srlavailstr = srlavailstr(:)';
            srlavailstr = ['At time of error, available ports are ',srlavailstr,' | '];
            msg = [srlavailstr, MEsrl.message];
            error(msg)
        else
            rethrow(MEsrl)
        end
    end
end
%}
receiverSerial.UserData = srlUD;
%srlString = 'Serial is connected on a parallel thread.';
srlString = 'Serial was bypassed on a parallel thread.';

% saving 
srlUD.TimeStamp = nan;
srlBuff = repmat(srlUD, [100, 1]);
srlLastMsg = srlUD.ReceivedData;

%% init 

dT = .001; % s between data requests 

% first request for raw data to get details only 
[~,~,~,rawN] = InitializeRawData(rawIDs, buffSize);
chID = cellfun(@(s) s.IDnumber, rawN); chID = chID(chInd);
buffSize = UserArgs.bufferSizeGrid .* ones(size(rawN)); % samples
buffSize(chInd) = UserArgs.bufferSize;
ShutdownRecording();

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

% assign channel indexes (NOT IDs!) to use for filtering and forecasting
selRaw2Flt = []; selFlt2For = []; 
if FilterSetUp
    selRaw2Flt = chInd; 
    if MdlSetUp
        selFlt2For = 1;
    end
end
selRaw2For = []; 

% initialize data structs for real
[rawD, fltD, forD, timeBuffs, initTic] = ...
    InitializeRecording(buffSize, filtOrds, forecastwin, ...
    rawIDs, selRaw2Flt, selRaw2For, selFlt2For);
rawN = rawD(1,:); 
% to do: double width of forD and implement sine wave !!
bgArgOut = UserArgs;

% define input args for filtering/forecasting funcs 
Fs = cellfun(@(s) s.SampleRate, rawN); Fs = Fs(selRaw2Flt);
if FilterSetUp
    fltN = fltD(1,:); 
    fIC = arrayfun(@(ord) zeros(ord,1), filtOrds, 'UniformOutput',false);
    filtArgs.fltInit = fIC; filtArgs.fltObj = filtObjs;
    filtArgs.TimeShift = filtOrds(1)/(2*Fs); 
    if MdlSetUp
        forN = forD(1,:);
        foreArgs.K = forecastwin; foreArgs.k = forecastpad;
        foreArgs.TimeStart = nan(size(forN));
        foreArgs.TimeShift = [zeros(size(selRaw2For)), filtArgs.TimeShift(selFlt2For)];
        foreArgs.ARmdls = mdls;
        foreArgs.SampleRates = Fs;
        foreArgs.FreqRange = [loco, hico];
        foreArgs.PhaseOfInterest = PhaseOfInterest;
    else
        foreArgs = [];
    end
else
    filtArgs = [];
    foreArgs = [];
end

% init forecast-output buffers, i.e. times to/of next phase(s) of interest
buffSize2 = (UserArgs.bufferSize / rawN{chInd}.SampleRate) * .5 * UserArgs.stimMaxFreq;
buffSize2 = ceil(buffSize2);
forBuffs = {nan(buffSize2,2)}; 
stimBuff = nan(buffSize2,1); 

% stimulator 
if UserArgs.StimActive
    StimArgs = SetupStimulator(UserArgs);
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
    ~, ~, ...
    fltD, filtArgs, ...
    forBuffs, forD, foreArgs] = ...
    iterReadBrain(...
        timeBuffs, rawD, @() GetNewRawData(rawIDs), ...
        selRaw2Flt, selRaw2For, selFlt2For, ...
        [], [], [], [], ...
        fltD, filtfun, filtArgs, ...
        forBuffs, forD, forefun, foreArgs);

    % do stimulus 
    forBuff = forBuffs{1}; 
    forBuffNew = [max(forBuff(:,1)), max(forBuff(:,2))]; 
    forBuffNew = forBuffNew - timeBuffs{chInd}(end,:); % [t2p, t2t]
    if doStim
        [t2stim, stim2q] = StimController(receiverSerial, UserArgs, fltD{4,1}(:,2), forBuffNew);
        if stim2q && (t2stim < 1.5*looptime)
            t2stim = .001*floor(1000*t2stim); % round to nearest 1ms 
            % ensure below max freq
            if 1/(t2stim + timeBuffs{chInd}(end,:) - stimLastTime) <= UserArgs.stimMaxFreq
                % check if stimulator is ready here?
                % assume stim should occur before completion of next loop
                pause(t2stim);
                stimtime1 = GetTime(initTic);
                StimArgs = PulseStimulator(StimArgs, UserArgs);
                stimtime2 = GetTime(initTic);
                stimtime = .5*(stimtime1 + stimtime2);
                stimBuff = bufferData(stimBuff, stimtime);
                stimLastTime = stimtime;
            end
        end
    end

    % handle serial comm 
    % TO DO: there should be a better way to do this; serial callback
    % should trigger an event or listener that logs the info 
    ReceivedData = receiverSerial.UserData.ReceivedData;
    if ~isempty(ReceivedData)
        srlString = ['Serial Port ',char(receiverSerial.Port),' Received Message: ',ReceivedData];
    end
    if ~strcmp(ReceivedData, srlLastMsg)
        ud = receiverSerial.UserData; 
        ud.TimeStamp = GetTime(initTic); % should this be last proc time from timeBuffs ???
        srlBuff = bufferData(srlBuff, ud);
    end
    srlLastMsg = ReceivedData;

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
        send(DQ, [{rawD_(1,:)}, fltD(1,1), forD(1,1); ...
                  {rawD_(4,:)}, fltD(4,1), forD(4,1); ...
                  {timeBuffs_}, {stimBuff}, {forBuff}; ...
                  {srlBuff}, {receiverSerial.UserData}, {srlString}]);
        srlBuff = repmat(srlUD, size(srlBuff)); % blank 
        srlString = '';
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
if UserArgs.StimActive
    StimArgs = ShutdownStimulator(StimArgs, UserArgs);
end

% If there are any errors in the full func, stop looping and allow the User
% to handle them. 
catch ME_fullfunc
    cont_fullfunc = false;
    getReport(ME_fullfunc)
    send(DQ, ME_fullfunc)
end
%end

%% function def 

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the filtered channel of interest
if isempty(foreArgs)
    error('Attempted forecast function with empty arguments.')
end
K = foreArgs.K; k = foreArgs.k;
ARmdls = foreArgs.ARmdls;
fs = foreArgs.SampleRates;
Ts = foreArgs.TimeShift;
Fco = foreArgs.FreqRange;
phis = foreArgs.PhaseOfInterest;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch_fore = 1:size(inData,2)
    armdl = ARmdls{ch_fore};
    FT = myFastForecastAR(armdl, inData{ch_fore}(:,2), K);
    foreTails{ch_fore} = [nan(height(FT),1), FT];
    foreTails{ch_fore}(1,1) = foreArgs.TimeStart(ch_fore);

    FT = FT(1:k,:); % use limited duration for hilbert padding
    [t2,i2,phi_inst,f_inst] = blockPDS(...
        inData{ch_fore}(:,2), FT, fs(ch_fore), phis, ...
        Ts(ch_fore), Fco(1), Fco(2));
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

end