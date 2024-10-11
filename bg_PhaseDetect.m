function bgArgOut = bg_PhaseDetect(UQ, DQ, SQ, ...
    InitializeRecording, ShutdownRecording, selRaw)
% Run brain recording with phase detection/prediction for PDS

try 

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

%% init 

dT = .001; % s between data requests 

selRaw2Flt = chInd; selRaw2For = []; selFlt2For = 1;

[rawD, fltD, forD, timeBuffs] = ...
    InitializeRecording(buffSize, filtorder, forecastwin, ...
    selRaw, selRaw2Flt, selRaw2For, selFlt2For);
rawN = rawD(1,:); fltN = fltD(1,:); forN = forD(1,:);

Fs = cellfun(@(s) s.SampleRate, rawN); Fs = Fs(selRaw2Flt);
foreArgs.K = forecastwin; foreArgs.k = ceil(.02*foreArgs.K);
filtorder = repmat(filtorder, size(fltN));
fIC = arrayfun(@(ord) zeros(ord,1), filtorder, 'UniformOutput',false);
filtArgs.fltInit = fIC; filtArgs.fltObj = filtwts;
filtArgs.TimeShift = repmat(TimeShiftFIR, size(fltN)); 
foreArgs.TimeStart = nan(size(forN));
foreArgs.TimeShift = [zeros(size(selRaw2For)), filtArgs.TimeShift(selFlt2For)];
foreArgs.ARmdls = mdls;
foreArgs.SampleRates = Fs;
foreArgs.FreqRange = [loco, hico];

forBuffs = cellfun(@(X) (nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);
forBuffs = forBuffs(chInd);

%% loop 
cont = true;
while cont
    pause(dT)

    try
    [...
    timeBuffs, rawD, ...
    ~, ~, ...
    fltD, filtArgs, ...
    forBuffs, forD, foreArgs] = ...
    iterReadBrain(...
        timeBuffs, rawD, @() getNewRawData_cbmex([]), ...
        selRaw2Flt, selRaw2For, selFlt2For, ...
        [], [], [], [], ...
        fltD, @filtFun, filtArgs, ...
        forBuffs, forD, @foreFun, foreArgs);

    if DQ.QueueLength == 0
        send(DQ, [rawD(1,chInd), fltD(1,1), forD(1,1); ...
                  rawD(4,chInd), fltD(4,1), forD(4,1); ...
                  timeBuffs(chInd),  {nan}, forBuffs(1)]);
    end

    catch ME_loop
        cont = false;
        getReport(ME_loop)
        send(DQ, ME_loop);
    end
end

%% stop 
ShutdownRecording();

catch ME_fullfunc
    getReport(ME_fullfunc)
    send(DQ, ME_fullfunc)
end

%% function def 

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the filtered channel of interest
K = foreArgs.K; k = foreArgs.k;
ARmdls = foreArgs.ARmdls;
fs = foreArgs.SampleRates;
Ts = foreArgs.TimeShift;
Fco = foreArgs.FreqRange;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch_ = 1:size(inData,2)
    armdl = ARmdls{ch_};
    FT = myFastForecastAR(armdl, inData{ch_}(:,2), K);
    foreTails{ch_} = [nan(height(FT),1), FT];
    foreTails{ch_}(1,1) = foreArgs.TimeStart(ch_);

    FT = FT(1:k,:); % use limited duration for hilbert padding
    [t2,i2,phi_inst,f_inst] = blockPDS(...
        inData{ch_}(:,2), FT, fs(ch_), [0,pi], ...
        Ts(ch_), Fco(1), Fco(2));
    t2 = t2-Ts(ch_); % [t2peak, t2trough]
    t2 = max(t2,0);
    foreBuffsAdd{ch_} = t2; 
end
end

function [fltTails, fltArgs] = filtFun(fltArgs, rawTails)
filtObj = fltArgs.fltObj;
filtInit = fltArgs.fltInit;
filtFin = cell(size(filtInit));
fltTails = cell(size(rawTails));
for ch__ = 1:size(rawTails,2)
    [FT,filtFin{ch__}] = filter(filtObj,1,rawTails{ch__}(:,2:end),filtInit{ch__});
    fltTails{ch__} = [rawTails{ch__}(:,1) - fltArgs.TimeShift(ch__), FT];
end
fltArgs.fltInit = filtFin;
end

end