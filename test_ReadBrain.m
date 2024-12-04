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

%% load AR model 
load("20240829_ARmdl.mat","ARmdl");

%% init 

dT = .001; % s between data requests 
PDSwin = 1000; % # samples ahead to forecast
buffSize = 2000; % samples

connect_cbmex(); 
t0 = datetime - seconds(cbmex('time'));
pause(10*dT);

svname = ['Saved Data Test',filesep,mfilename,'_SaveFile_',datestr(t0,'yyyymmdd_HHMMSS')];
svN = 1;

try
    [rawH, rawT, rawB, rawN] = initRawData_cbmex([], buffSize);
catch ME
    warning(['Error on first attempt: ',ME.message]);
    pause(1);
    [rawH, rawT, rawB, rawN] = initRawData_cbmex([], buffSize);
end
Fs = cellfun(@(s) s.SampleRate, rawN);
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
filtArgs.TimeShift = repmat(TimeShiftFIR, size(rawH)); 
foreArgs.TimeStart = nan(size(rawH));
foreArgs.TimeShift = filtArgs.TimeShift;
foreArgs.ARmdls = {ARmdl};
foreArgs.SampleRates = Fs;
foreArgs.FreqRange = [loco, hico];

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

chInd = 33;
selRaw2Flt = chInd; selRaw2For = []; selFlt2For = 1;

fltD = fltD(:,chInd); forD = forD(:,chInd);
forBuffs = forBuffs(chInd); forBuffRow = ones(size(forBuffs));
forBuffStore = nan(10*height(forBuffs{1}), width(forBuffs{1}));

fig = figure; 

ax(2) = subplot(2,1,2); 
tPlt = timeBuffs{chInd}; tPltDisp = tPlt; td0 = tPltDisp(end); tDisp1 = tic;
pltTime = stem(t0 + seconds(tPlt), [nan; diff(tPlt)], '--');
grid on; hold on;
pltTimeDisp = stem(t0 + seconds(tPltDisp), [nan; diff(tPltDisp)]);
ylabel('Dur (s)');

% add plotting code below to be output from data2timetable? 
ax(1) = subplot(2,1,1);
rawPlt = data2timetable(rawB(chInd),rawN(chInd),t0); rawPlt = rawPlt{1};
fltPlt = data2timetable(fltB(chInd),rawN(chInd),t0); fltPlt = fltPlt{1};
forPlt = data2timetable(forB(chInd),rawN(chInd),t0); forPlt = forPlt{1};
pltRaw = plot(rawPlt.Time, rawPlt.Variables); 
hold on; grid on; 
xlabel('time'); ylabel(rawPlt.Properties.VariableNames{1});
pltFlt = plot(fltPlt.Time, fltPlt.Variables);
pltFor = plot(forPlt.Time, forPlt.Variables);

tPk = forBuffs{1}(:,1); tTr = forBuffs{1}(:,2);
pltPk = plot(t0 + seconds(tPk), zeros(size(tPk)), '^m');
pltTr = plot(t0 + seconds(tTr), zeros(size(tTr)), 'vm');

%linkaxes(ax, 'x')

%% loop 
cont = isvalid(fig);
while cont
    pause(dT)
    tDisp2 = td0 + toc(tDisp1); %tDisp1 = tic;
    tPltDisp = bufferData(tPltDisp, tDisp2);

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

    forBuffs1 = forBuffs{1};
    [forBuffFull, forBuffStore, forBuffRow, forBuffedOut] = ...
        bufferStorage(forBuffStore, forBuffRow, forBuffs1);
    if forBuffFull
        PeakTrough = forBuffedOut;
        save([svname,'_',num2str(svN),'.mat'], 'PeakTrough');
        svN = svN+1;
        forBuffs1 = [forBuffedOut; forBuffs1];
    end

    rawPlt = data2timetable(rawD(4,chInd),rawD(1,chInd),t0); rawPlt = rawPlt{1};
    fltPlt = data2timetable(fltD(4,1),fltD(1,1),t0); fltPlt = fltPlt{1};
    forPlt = data2timetable(forD(4,1),forD(1,1),t0); forPlt = forPlt{1};
    tPlt = timeBuffs{chInd};
    tPk = forBuffs1(:,1); tTr = forBuffs1(:,2);
    tPltRng = [rawPlt.Time; fltPlt.Time; forPlt.Time];
    tPltRng = [min(tPltRng), max(tPltRng)];
    tPltRng = tPltRng + [-1,1]*.1*diff(tPltRng);

    cont = isvalid(fig);
    if cont
        figure(fig);
        pltRaw.YData = rawPlt.Variables; pltRaw.XData = rawPlt.Time;
        pltFlt.YData = fltPlt.Variables; pltFlt.XData = fltPlt.Time;
        pltFor.YData = forPlt.Variables; pltFor.XData = forPlt.Time;
        % subplot(2,1,1); tPltRng = xlim();
        pltPk.XData = t0 + seconds(tPk); pltPk.YData = zeros(size(tPk));
        pltTr.XData = t0 + seconds(tTr); pltTr.YData = zeros(size(tTr));
        pltTime.YData = [nan; diff(tPlt)]; pltTime.XData = t0 + seconds(tPlt);
        pltTimeDisp.YData = [nan; diff(tPltDisp)]; pltTimeDisp.XData = t0 + seconds(tPltDisp);
        if sum(~isnat(tPltRng))
            subplot(2,1,1); xlim(tPltRng);
            subplot(2,1,2); xlim(tPltRng);
        end
    end

    catch ME
        cont = false;
        getReport(ME)
    end
end

%% close 
disconnect_cbmex();
PeakTrough = forBuffs{1};
save(svname,"PeakTrough");

%% function def 

function [foreTails, foreBuffsAdd, foreArgs] = foreFun(foreArgs, inData)
% inData should just be the filtered channel of interest
k = foreArgs.k;
ARmdls = foreArgs.ARmdls;
Fs = foreArgs.SampleRates;
Ts = foreArgs.TimeShift;
Fco = foreArgs.FreqRange;
foreTails = cell(size(inData)); foreBuffsAdd = foreTails; 
for ch = 1:size(inData,2)
    ARmdl = ARmdls{ch};
    FT = myFastForecastAR(ARmdl, inData{ch}(:,2), k);
    foreTails{ch} = [nan(height(FT),1), FT];
    foreTails{ch}(1,1) = foreArgs.TimeStart(ch);

    [t2,i2,phi_inst,f_inst] = blockPDS(...
        inData{ch}(:,2), FT, Fs(ch), [0,pi], ...
        Ts(ch), Fco(1), Fco(2));
    t2 = t2-Ts(ch); % [t2peak, t2trough]
    t2 = max(t2,0);
    foreBuffsAdd{ch} = t2; 
end
end

function [fltTails, fltArgs] = filtFun(fltArgs, rawTails)
filtObj = fltArgs.fltObj;
filtInit = fltArgs.fltInit;
filtFin = cell(size(filtInit));
fltTails = cell(size(rawTails));
for ch = 1:size(rawTails,2)
    [FT,filtFin{ch}] = filter(filtObj,1,rawTails{ch}(:,2:end),filtInit{ch});
    fltTails{ch} = [rawTails{ch}(:,1) - fltArgs.TimeShift(ch), FT];
end
fltArgs.fltInit = filtFin;
end