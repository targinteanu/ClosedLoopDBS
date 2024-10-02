%% set filter 
minfac         = 1;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length
srate = 1000; 
loco = 13; hico = 30; % Hz 
filtorder = minfac*fix(srate/loco);
TimeShiftFIR = filtorder/(2*srate); % seconds
IndShiftFIR = ceil(filtorder/2); % samples

% build filter 
% usage i.e.: 
% >> filteredSignal = filter(filtwts, 1, unfilteredSignal) 
% -- OR --
% >> filteredSignal = filtfilt(filtwts, 1, unfilteredSignal)
filtwts = fir1(filtorder, [loco, hico]./(srate/2));

%% init 

dT = .01; % s between data requests 
N = round(3/dT); % # of data samples to obtain 
PDSwin = 1000; % # samples ahead to forecast
buffSize = 30000; % samples

connect_cbmex(); 
t0 = datetime - seconds(cbmex('time'));
pause(dT);

[rawH, rawT, rawB] = initRawData_cbmex([], buffSize, t0);
fltH = initFilteredData(rawH, repmat(IndShiftFIR, size(rawH))); 
[forH, forT, forB] = initForecastData(fltH, repmat(PDSwin, size(fltH)));

rawN = cellfun(@(T) T.Properties.VariableNames{1}, rawH, 'UniformOutput',false);
rawD = [rawN; rawH; rawT; rawB]; 
forD = [rawN; forH; forT; forB];

fltT = cellfun(@(T) multTbl(0,T), rawT, 'UniformOutput',false);
fltB = cellfun(@(T) multTbl(0,T), rawB, 'UniformOutput',false);
fltD = [rawN; fltH; fltT; fltB];

foreArgs.k = PDSwin;
fIC = zeros(filtorder,1); fIC = repmat({fIC}, size(rawH)); 
filtArgs.fltInit = fIC; filtArgs.fltObj = filtwts;

timeBuffs = cellfun(@(T) T.Time, rawH, 'UniformOutput',false);
forBuffs = cellfun(@(X) t0+seconds(nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);

selRaw2Flt = 1:length(rawN); selRaw2For = []; selFlt2For = selRaw2Flt;

fig = figure; 

% replace below with "snapshot data" channel selector
chInd = 1;
pltRaw = plot(rawB{chInd}.Time, rawB{chInd}.Variables); 
hold on; grid on; 
xlabel('time'); ylabel(rawB{chInd}.Properties.VariableNames{1});
pltFlt = plot(fltB{chInd}.Time, fltB{chInd}.Variables);
pltFor = plot(forB{chInd}.Time, forB{chInd}.Variables);

%% loop 
cont = isvalid(fig);
while cont
    pause(.001)

    try
    [...
    timeBuffs, rawD, ...
    ~, ~, ...
    fltD, filtArgs, ...
    forBuffs, forD, foreArgs] = ...
    iterReadBrain(...
        timeBuffs, rawD, @() getNewRawData_cbmex([], t0), ...
        selRaw2Flt, selRaw2For, selFlt2For, ...
        [], [], [], [], ...
        fltD, @filtFun, filtArgs, ...
        forBuffs, forD, @foreFun, foreArgs);

    pltRaw.YData = rawD{4,chInd}.Variables; pltRaw.XData = rawD{4,chInd}.Time;
    pltFlt.YData = fltD{4,chInd}.Variables; pltFlt.XData = fltD{4,chInd}.Time;
    pltFor.YData = forD{4,chInd}.Variables; pltFor.XData = forD{4,chInd}.Time;

    cont = isvalid(fig);

    catch ME
        cont = false;
        getReport(ME)
    end
end

%% close 
disconnect_cbmex();

%% function def 

function Yf = mySimpleForecast(Yp, k)
% zero-order interp; to be replaced with real 
yf = repmat(Yp{end,:},k,1);
dt = Yp.Properties.TimeStep;
Yf = timetable(yf, ...
    'TimeStep',dt, ...
    'StartTime', Yp.Time(end)+dt, ...
    'VariableNames',Yp.Properties.VariableNames);
end

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the filtered channel of interest
k = foreArgs.k;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch = 1:size(inData,2)
    foreTails{ch} = mySimpleForecast(inData{ch}, k);
    foreBuffsAdd{ch} = [rand,rand]*2 + foreTails{ch}.Properties.StartTime; % replace with time of peak, trough
end
end

function [fltTails, fltArgs] = filtFun(fltArgs, rawTails)
filtObj = fltArgs.fltObj;
filtInit = fltArgs.fltInit;
filtFin = cell(size(filtInit));
fltTails = cell(size(rawTails));
for ch = 1:size(rawTails,2)
    [fltTails{ch},filtFin{ch}] = FilterTimetable(filtObj,1,rawTails{ch},filtInit{ch});
end
fltArgs.fltInit = filtFin;
end