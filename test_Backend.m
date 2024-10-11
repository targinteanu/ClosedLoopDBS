function test_Backend(DQ)

try

%% set filter 
minfac         = 1;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length
srate = 1000; 
loco = 13; hico = 30; % Hz 
filtorder = minfac*fix(srate/loco);
TimeShiftFIR = filtorder/(2*srate); % seconds

% build filter 
% usage i.e.: 
% >> filteredSignal = filter(filtwts, 1, unfilteredSignal) 
% -- OR --
% >> filteredSignal = filtfilt(filtwts, 1, unfilteredSignal)
filtwts = fir1(filtorder, [loco, hico]./(srate/2));

%% load AR model 
load("20240829_ARmdl.mat","ARmdl");

%% init 

dT = .001; % s between data requests 
forecastwin = 1000; % # samples ahead to forecast
buffSize = 20000; % samples

Fs = cellfun(@(s) s.SampleRate, rawN);
foreArgs.k = forecastwin;
fIC = zeros(filtorder,1); fIC = repmat({fIC}, size(rawH)); 
filtArgs.fltInit = fIC; filtArgs.fltObj = filtwts;
filtArgs.TimeShift = repmat(TimeShiftFIR, size(rawH)); 
foreArgs.TimeStart = nan(size(rawH));
foreArgs.TimeShift = filtArgs.TimeShift;
foreArgs.ARmdls = {ARmdl};
foreArgs.SampleRates = Fs;
foreArgs.FreqRange = [loco, hico];

chInd = 33;
selRaw2Flt = chInd; selRaw2For = []; selFlt2For = 1;

[rawD, fltD, forD, timeBuffs] = ...
    InitializeRecording_cbmex(buffSize, filtorder, forecastwin, ...
    [], selRaw2Flt, selRaw2For, selFlt2For);

%forBuffs = cellfun(@(X) t0+seconds(nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);
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

    catch ME
        cont = false;
        getReport(ME)
        send(DQ, ME);
    end
end

disconnect_cbmex();

catch ME
    getReport(ME)
    send(DQ, ME);
end

%% function def 

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the filtered channel of interest
k = foreArgs.k;
ARmdls = foreArgs.ARmdls;
fs = foreArgs.SampleRates;
Ts = foreArgs.TimeShift;
Fco = foreArgs.FreqRange;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch_ = 1:size(inData,2)
    armdl = ARmdls{ch_};
    FT = myFastForecastAR(armdl, inData{ch_}(:,2), k);
    foreTails{ch_} = [nan(height(FT),1), FT];
    foreTails{ch_}(1,1) = foreArgs.TimeStart(ch_);

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