function bgArgOut = bg_PhaseDetect(UserArgs, DQ, SQ, ...
    InitializeRecording, ShutdownRecording, ...
    InitializeRawData, GetNewRawData, selRaw)
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
% [ raw , filtered , forecast, timing buffer ] data
% 
% ShutdownRecording takes no arguments.
% 
% InitializeRawData takes ( selRaw , BufferSize ) as input arguments and
% returns [ EmptyDataStruct , Tail , Combined , ChannelInfo ]
% 
% GetNewRawData takes selRaw as an argument and returns new raw tails. 
% 

bgArgOut = [];

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
MdlSetUp = UserArgs.MdlSetUp; % t/f
if FilterSetUp
    filtOrds = [UserArgs.filtorder]; % array with chans as cols
    filtObjs = {1; UserArgs.BPF}; % cell with rows {a; b}; chans as cols
    hico = UserArgs.hicutoff; loco = UserArgs.locutoff; % Hz 
else
    filtOrds = [];
end
if MdlSetUp
    mdls = {UserArgs.Mdl};
end
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

%% init 

dT = .001; % s between data requests 

% first request for raw data to get details only 
[~,~,~,rawN] = InitializeRawData(selRaw, buffSize);
chID = cellfun(@(s) s.IDnumber, rawN); chID = chID(chInd);
buffSize = UserArgs.bufferSizeGrid .* ones(size(rawN)); % samples
buffSize(chInd) = UserArgs.bufferSize;
ShutdownRecording();

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
[rawD, fltD, forD, timeBuffs] = ...
    InitializeRecording(buffSize, filtOrds, forecastwin, ...
    selRaw, selRaw2Flt, selRaw2For, selFlt2For);
rawN = rawD(1,:); 
% to do: double width of forD and implement sine wave !!

% define input args for filtering/forecasting funcs 
Fs = cellfun(@(s) s.SampleRate, rawN); Fs = Fs(selRaw2Flt);
if FilterSetUp
    fltN = fltD(1,:); 
    fIC = arrayfun(@(ord) zeros(ord,1), filtOrds, 'UniformOutput',false);
    filtArgs.fltInit = fIC; filtArgs.fltObj = filtObjs;
    filtArgs.TimeShift = filtOrds(1)/Fs; 
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
forBuffs = cellfun(@(X) (nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);
forBuffs = forBuffs(chInd); forBuffRow = ones(size(forBuffs));

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

%% loop 
cont_loop = true;
while cont_loop
    pause(dT)

    try
    % main iteration 
    [...
    timeBuffs, rawD, ...
    ~, ~, ...
    fltD, filtArgs, ...
    forBuffs, forBuffRow, forBuffedOut, forD, foreArgs] = ...
    iterReadBrain(...
        timeBuffs, rawD, @() GetNewRawData(selRaw), ...
        selRaw2Flt, selRaw2For, selFlt2For, ...
        [], [], [], [], ...
        fltD, filtfun, filtArgs, ...
        forBuffs, forBuffRow, forD, forefun, foreArgs);

    % DataQueue is empty when the User polls it, which means the
    % User is ready for new data. 
    if DQ.QueueLength == 0
        if isempty(selRaw)
            rawD_ = rawD;
            timeBuffs_ = timeBuffs;
        else
            rawD_ = rawD(:,selRaw);
            timeBuffs_ = timeBuffs(selRaw);
        end
        send(DQ, [{rawD_(1,:)}, fltD(1,1), forD(1,1); ...
                  {rawD_(4,:)}, fltD(4,1), forD(4,1); ...
                  {timeBuffs_}, forBuffedOut(1), forBuffs(1)]);
    end

    cont_loop = cont_loop && cont_loop_2;

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
    t2 = max(t2,0);
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