function bgArgOut = bg_PhaseDetect(UQ, DQ, SQ, ...
    InitializeRecording, ShutdownRecording, selRaw)
% Run brain recording with phase detection/prediction for PDS

bgArgOut = [];

cont_fullfunc = true; % run or wait for user input
while cont_fullfunc
try
%% poll queue(s) for data 
poll_timeout_sec = 3600;
[UserArgs, ok] = poll(UQ, poll_timeout_sec);
if ~ok
    error('Background process timed out without user input.')
end

%% import filter and model details, etc
% change names according to front-end handles/app/struct!!
filtOrds = UserArgs.filtOrds; % cell with chans as cols
filtObjs = UserArgs.filtObjs; % cell with rows {a; b}; chans as cols
hico = UserArgs.hico; loco = UserArgs.loco; % Hz 
mdls = UserArgs.mdls;
chInd = UserArgs.chInd;
forecastwin = UserArgs.PDSwin1; % # samples ahead to forecast
forecastpad = UserArgs.PDSwin2; % # of above to use to pad hilbert transform
buffSize = UserArgs.bufferSize; % samples
PhaseOfInterest = UserArgs.PhaseOfInterest;

%% init 

dT = .001; % s between data requests 
TimeShiftFIR = filtorder/(2*srate); % seconds
selRaw2Flt = chInd; selRaw2For = []; selFlt2For = 1;

[rawD, fltD, forD, timeBuffs] = ...
    InitializeRecording(buffSize, filtOrds, forecastwin, ...
    selRaw, selRaw2Flt, selRaw2For, selFlt2For);
rawN = rawD(1,:); fltN = fltD(1,:); forN = forD(1,:);

Fs = cellfun(@(s) s.SampleRate, rawN); Fs = Fs(selRaw2Flt);
foreArgs.K = forecastwin; foreArgs.k = forecastpad;
fIC = arrayfun(@(ord) zeros(ord,1), filtOrds, 'UniformOutput',false);
filtArgs.fltInit = fIC; filtArgs.fltObj = filtObjs;
filtArgs.TimeShift = TimeShiftFIR; 
foreArgs.TimeStart = nan(size(forN));
foreArgs.TimeShift = [zeros(size(selRaw2For)), filtArgs.TimeShift(selFlt2For)];
foreArgs.ARmdls = mdls;
foreArgs.SampleRates = Fs;
foreArgs.FreqRange = [loco, hico];
foreArgs.PhaseOfInterest = PhaseOfInterest;

forBuffs = cellfun(@(X) (nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);
forBuffs = forBuffs(chInd); forBuffRow = ones(size(forBuffs));

%% loop 
cont_loop = true;
while cont_loop
    pause(dT)

    try
    [...
    timeBuffs, rawD, ...
    ~, ~, ...
    fltD, filtArgs, ...
    forBuffs, forBuffRow, forBuffedOut, forD, foreArgs] = ...
    iterReadBrain(...
        timeBuffs, rawD, @() getNewRawData_cbmex([]), ...
        selRaw2Flt, selRaw2For, selFlt2For, ...
        [], [], [], [], ...
        fltD, @filtFun, filtArgs, ...
        forBuffs, forBuffRow, forD, @foreFun, foreArgs);

    % DataQueue is empty when the User polls it, which means the
    % User is ready for new data. 
    if DQ.QueueLength == 0
        send(DQ, [rawD(1,chInd), fltD(1,1), forD(1,1); ...
                  rawD(4,chInd), fltD(4,1), forD(4,1); ...
                  timeBuffs(chInd), forBuffedOut(1), forBuffs(1)]);
    end

    % If the UserQueue is non-empty, there are new UserArgs, and the full
    % func should be restarted. 
    cont_loop = cont_loop && UQ.QueueLength < 1;

    % If there are any errors in the loop, stop looping and allow the
    % User to handle them. 
    catch ME_loop
        cont_loop = false;
        getReport(ME_loop)
        send(DQ, ME_loop);
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
end

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