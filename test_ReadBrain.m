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
forH = initForecastData(fltH, PDSwin);

rawN = cellfun(@(T) T.Properties.VariableNames{1}, rawH, 'UniformOutput',false);
rawD = [rawN; rawH; rawT; rawB]; 

%% function def 

function Yf = mySimpleForecast(Yp, k)
% zero-order interp; to be replaced with real 
yf = repmat(Yp.Variables(end),k,1);
dt = Yp.Properties.TimeStep;
Yf = timetable(yf, ...
    'TimeStep',dt, ...
    'StartTime', Yp.Time(end)+dt);
end

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the filtered channel of interest
k = foreArgs.k;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch = 1:size(inData,2)
    foreTails{ch} = mySimpleForecast(inData{ch}, k);
    foreBuffsAdd = [rand,rand]*2 + foreTails{ch}.StartTime; % replace with time of peak, trough
end
end