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

dT = .001; % s between data requests 
PDSwin = 1000; % # samples ahead to forecast
buffSize = 30000; % samples

connect_cbmex(); 
t0 = datetime - seconds(cbmex('time'));
pause(dT);

try
    [rawH, rawT, rawB, rawN] = initRawData_cbmex([], buffSize);
catch ME
    warning(['Error on first attempt: ',ME.message]);
    pause(1);
    [rawH, rawT, rawB, rawN] = initRawData_cbmex([], buffSize);
end
fltH = initFilteredData(rawH, repmat(IndShiftFIR, size(rawH))); 
[forH, forT, forB] = initForecastData(fltH, repmat(PDSwin, size(fltH)));

rawD = [rawN; rawH; rawT; rawB]; 
forD = [rawN; forH; forT; forB];

fltT = cellfun(@(D) [1,0].*D, rawT, 'UniformOutput',false);
fltB = cellfun(@(D) [1,0].*D, rawB, 'UniformOutput',false);
fltD = [rawN; fltH; fltT; fltB];

foreArgs.k = PDSwin;
fIC = zeros(filtorder,1); fIC = repmat({fIC}, size(rawH)); 
filtArgs.fltInit = fIC; filtArgs.fltObj = filtwts;

timeBuffs = cell(size(rawH));
for ch = 1:length(rawH)
    t_ch = rawH{ch}(:,1);
    if isnan(t_ch(end))
        warning([rawN{ch}.Name,' timestamp 0 has not been assigned!'])
    end
    dt_ch = 1/rawN{ch}.SampleRate;
    for it = 2:length(t_ch)
        if isnan(t_ch(it))
            t_ch(it) = t_ch(it-1) + dt_ch;
        end
    end
    timeBuffs{ch} = t_ch;% + t0;
end
%forBuffs = cellfun(@(X) t0+seconds(nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);
forBuffs = cellfun(@(X) (nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);

chInd = 1;
selRaw2Flt = chInd; selRaw2For = []; selFlt2For = chInd;

fltD = fltD(:,chInd); forD = forD(:,chInd);
forBuffs = forBuffs(chInd);

fig = figure; 

% add plotting code below to be output from data2timetable? 
rawPlt = data2timetable(rawB(chInd),rawN(chInd),t0); rawPlt = rawPlt{1};
fltPlt = data2timetable(fltB(chInd),rawN(chInd),t0); fltPlt = fltPlt{1};
forPlt = data2timetable(forB(chInd),rawN(chInd),t0); forPlt = forPlt{1};
pltRaw = plot(rawPlt.Time, rawPlt.Variables); 
hold on; grid on; 
xlabel('time'); ylabel(rawPlt.Properties.VariableNames{1});
pltFlt = plot(fltPlt.Time, fltPlt.Variables);
pltFor = plot(forPlt.Time, forPlt.Variables);

%% loop 
cont = isvalid(fig);
while cont
    pause(10*dT)

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

    rawPlt = data2timetable(rawD(4,chInd),rawD(1,chInd),t0); rawPlt = rawPlt{1};
    fltPlt = data2timetable(fltD(4,chInd),rawD(1,chInd),t0); fltPlt = fltPlt{1};
    forPlt = data2timetable(forD(4,chInd),rawD(1,chInd),t0); forPlt = forPlt{1};

    cont = isvalid(fig);
    if cont
        pltRaw.YData = rawPlt.Variables; pltRaw.XData = rawPlt.Time;
        pltFlt.YData = fltPlt.Variables; pltFlt.XData = fltPlt.Time;
        pltFor.YData = forPlt.Variables; pltFor.XData = forPlt.Time;
    end

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
Yf = repmat(Yp(end,:),k,1);
end

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the filtered channel of interest
k = foreArgs.k;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch = 1:size(inData,2)
    foreTails{ch} = mySimpleForecast(inData{ch}, k);
    t = inData{ch}(:,1); 
    %{
    indf = find(~isnan(t)); indf = indf(end); 
    tf = t(indf); % last logged time
    L = height(t) - indf + 1; % how many samples between ^ and now
    %}
    foreBuffsAdd{ch} = [rand,rand]*2; % replace with ind of peak, trough
end
end

function [fltTails, fltArgs] = filtFun(fltArgs, rawTails)
filtObj = fltArgs.fltObj;
filtInit = fltArgs.fltInit;
filtFin = cell(size(filtInit));
fltTails = cell(size(rawTails));
for ch = 1:size(rawTails,2)
    [fltTails{ch},filtFin{ch}] = filter(filtObj,1,rawTails{ch},filtInit{ch});
end
fltArgs.fltInit = filtFin;
end