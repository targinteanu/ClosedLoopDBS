%% User selects folder; MATLAB loads all files 
filepath = uigetdir('Saved Data PD'); 
OnlineFiles = dir([filepath,filesep,'OnlineDisplaySavedData*.mat']);
OnlineFile = OnlineFiles(1); 
NSFiles = dir([filepath,filesep,'*.ns*']); 
NSFile = NSFiles(1); 
cd00 = cd; cd(filepath); 
load(OnlineFile.name); ns = openNSx(NSFile.name, 'uV'); 
cd(cd00); 

[dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS(ns);
dataOneChannelWithArtifact = dataOneChannel; 
t0 = t(1);

stimtgtphase = inputdlg('Stimulation target phase (degrees):');
stimtgtphase = (str2double(stimtgtphase));
stimtgtphase = stimtgtphase*pi/180; % rad

%% Get indexes of peaks, troughs, and stimulus pulses 
PeakInd = PeakTime*SamplingFreq; 
TroughInd = TroughTime*SamplingFreq; 
StimInd = StimTime*SamplingFreq; 
PeakInd = round(PeakInd); TroughInd = round(TroughInd); StimInd = round(StimInd); 
trimToSize = @(inds, x) inds((inds>=0) & (inds <= length(x)));
PeakInd = trimToSize(PeakInd,t); TroughInd = trimToSize(TroughInd,t); StimInd = trimToSize(StimInd,t);

%% artifact (outlier) detection 
%{
olStartEnd = [false, isoutlier(diff(dataOneChannel), 'mean')];
olStartEnd = find(olStartEnd); 
olStart = olStartEnd(1:2:end); olEnd = olStartEnd(2:2:end);
%}
dataOneChannel = dataOneChannelWithArtifact;
artExtend = 11; % extend artifact by __ samples 
artBegin = 2; % begin artifact __ samples before stim detection
io = isoutlier(dataOneChannel, 'mean');
artIndAll = StimTrainRec;
%artIndAll = io | StimTrainRec; 
%artIndAll(StimInd) = true;
artIndAll = movsum(artIndAll, artExtend) > 0;
artIndAll = find(artIndAll);
artIndAll = artIndAll + ceil(artExtend/2) - artBegin;
artIndAll = artIndAll(artIndAll > 1); 
artIndAll = artIndAll(artIndAll < length(dataOneChannel));
artIndAll = [1, artIndAll, length(dataOneChannel)];
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); baselineStartInd = artIndAll(baselineStartInd); 
baselineStartInd = 1e6; baselineEndInd = 1.45e6;

%% set baseline to fit model 
if ~isempty(artIndAll)
baselineWinLen = 1000; ARlen = 50; % samples 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
DCOS = mean(dataBaseline); dataBaseline = dataBaseline - DCOS;
%dataBaseline = Myeegfilt(dataBaseline,SamplingFreq,13,30, 0, 1024);
baselineWin = (baselineEndInd-baselineStartInd) + [-1,1]*baselineWinLen; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); baselineWin(2) = min(length(dataBaseline),baselineWin(2));
dataBaseline = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl = ar(iddata(dataBaseline', [], 1/SamplingFreq), ARlen, 'yw');

%% remove artifact 
dataOneChannel = dataOneChannelWithArtifact;
dataOneChannel = dataOneChannel - DCOS; % correct DC offset

for ind = artIndAll
    ind0 = ind - ARlen;
    if ind0 > 0
        % using AR mdl forecast:
        dataOneChannel(ind) = myFastForecastAR(ARmdl, dataOneChannel(ind0:(ind-1))', 1); 
        % linear interp (non-causal): 
        %dataOneChannel(ind) = nan;
    end
end

% linear interp (non-causal)
tRelNoArt = tRel(~isnan(dataOneChannel));
dataOneChannel = dataOneChannel(~isnan(dataOneChannel));
[dataOneChannel] = interp1(tRelNoArt, dataOneChannel, tRel, 'linear', 'extrap');

% plot artifact removal 
figure; 
ax(1) = subplot(211); 
plot(t, dataOneChannelWithArtifact-DCOS, 'k', 'LineWidth',1.5); 
grid on; hold on; 
plot(t, dataOneChannel, 'b'); 
title('Artifact Removal'); ylabel(channelName);
ax(2) = subplot(212); 
plot(t, dataAllChannels(channelIndexStim,:)); 
ylabel('ainp1');
grid on; linkaxes(ax, 'x'); 
end

%% filter 
myFilt = buildFIRBPF(SamplingFreq,13,30, 8);
%dataOneChannel = filtfilt(myFilt,1,dataOneChannel);
dataOneChannel = filter(myFilt,1,dataOneChannel);
myFiltShift = ceil(length(myFilt)/2);
dataOneChannel = [dataOneChannel(myFiltShift:end), zeros(1,myFiltShift-1)];

% detect power threshold 
pwrthresh = sqrt(bandpower(dataBaseline,SamplingFreq,[13,30]));
pwrthresh = .75*pwrthresh;

%% find ideal stim pattern 
% Get inst. freq. and phase 
[dataPhase, dataFreq] = instPhaseFreq(dataOneChannel, SamplingFreq);

dtgt = radfix(dataPhase-stimtgtphase);
[~,~,tgtInd] = zerocrossrate(dtgt, "TransitionEdge","rising");

% apply power threshold 
tgtInd = tgtInd & (envelope(dataOneChannel) > pwrthresh);

% exclude stim off times 
indSel = true(size(StimTrainRec));
StimIndRec = find(StimTrainRec); 
indSel(1:(StimIndRec(1)-1)) = false; % before first stim 
indSel((StimIndRec(end)+1):end) = false; % after last stim 
stimTimeAbs = t(StimIndRec);
iOff = diff(stimTimeAbs) > seconds(5);
for ind = find(iOff)
    % stimTimeAbs(ind+1) is more than 5 seconds from stimTimeAbs(ind)
    ind1 = StimIndRec(ind); ind2 = StimIndRec(ind+1);
    indSel(ind1:ind2) = false;
end
tgtInd = tgtInd & indSel;

tgtTimeAbs = t(tgtInd);
tgtISI = seconds([nan, diff(tgtTimeAbs)]);
stimISI = seconds([nan,diff(stimTimeAbs)]);

tgtTimeAbs = tgtTimeAbs(tgtTimeAbs >= stimTimeAbs(1)); % after stim began
tgtTimeAbs = tgtTimeAbs(tgtTimeAbs <= stimTimeAbs(end)); % before stim ended

%% determine red/yellow/green phases of experiment 
if ~isempty(SerialLog)
expStates = {SerialLog.ParadigmPhase}';
SrlTimes = [SerialLog.TimeStamp]';
expStatesU = unique(expStates);
expStateIT = cell(length(expStatesU),4);
for ESi = 1:length(expStatesU)
    expState = expStatesU{ESi};
    [ind, tStart, tEnd] = findExpState(expState, expStates, SrlTimes, tRel);
    indStim = StimInd(ind(StimInd));
    indStimRec = find(ind & StimTrainRec);
    tStart = seconds(tStart) + t0; 
    tEnd   = seconds(tEnd)   + t0;
    expStateIT{ESi,1} = ind; expStateIT{ESi,4} = indStim;
    expStateIT{ESi,2} = tStart; expStateIT{ESi,3} = tEnd;
    expStateIT{ESi,5} = indStimRec;
    clear expState ind tStart tEnd indStim ESi
end
else
    expStatesU = {};
    expStateIT = {};
end

%% Plot time series 

lgd = {'Data'; 'Peaks'; 'Troughs'; 'Stim Sent'; 'Stim Recd'; 'Tgt Stim'}; 
lgdsel = true(size(lgd));

% select the data to plot 
dataMinMax = [min(dataOneChannel), max(dataOneChannel)];
dataMinMax = dataMinMax + [-1,1]*.01*diff(dataMinMax);
dataMin = dataMinMax(1); dataMax = dataMinMax(2); 

% plot data and indicate recorded peaks, troughs, stimuli 
figure; plot(t,dataOneChannel); grid on; hold on; 
if isempty(PeakInd)
    lgdsel(2) = false;
else
    plot(t(PeakInd),   dataOneChannel(PeakInd),   '^m'); 
end
if isempty(TroughInd)
    lgdsel(3) = false;
else
    plot(t(TroughInd), dataOneChannel(TroughInd), 'vm'); 
end
if isempty(StimInd)
    lgdsel(4) = false;
else
    plot(t(StimInd),   dataOneChannel(StimInd),   '.r');
end
if isempty(StimTrainRec)
    lgdsel(5) = false;
else
    plot(t(StimTrainRec), dataOneChannel(StimTrainRec), 'or');
end
if isempty(tgtInd)
    lgdsel(6) = false;
else
    plot(t(tgtInd), dataOneChannel(tgtInd), '+r');
end
lgd = lgd(lgdsel);

if ~isempty(SerialLog)
% shade plot regions indicating encode and decode state of paradigm 
for ESi = 1:length(expStatesU)
    colr = expStatesU{ESi};
    if strcmp(colr,'gray')
        colr = [.5,.5,.5];
    elseif strcmp(colr,'Started')
        colr = [1,1,1];
    elseif strcmp(colr,'Stopped')
        colr = [0,0,0];
    end

    tStart = expStateIT{ESi,2}; tEnd = expStateIT{ESi,3};
    if isempty(tStart)
        tStart = nan;
    end
    if isempty(tEnd)
        tEnd = nan;
    end
    patch([tStart, tEnd, tEnd, tStart]', ...
        repmat([dataMax; dataMax; dataMin; dataMin],1,length(tStart)), ...
        colr, 'FaceAlpha', .2);

    clear colr tStart tEnd ESi
end
end

% label the plot 
legend([lgd; expStatesU])

%% Plot polar histogram 
figure; 
subplot(121); errHistoPolar(dataPhase(PeakInd),0,18); 
title('Predicted Peaks');
subplot(122); errHistoPolar(dataPhase(TroughInd),pi,18); 
title('Predicted Troughs');

figure; 
subplot(121); errHistoPolar(dataPhase(StimInd),stimtgtphase,18); 
title('Stim Sent');
subplot(122); errHistoPolar(dataPhase(StimTrainRec),stimtgtphase,18); 
title('Stim Received');

%% stim polar histogram
H = 2; W = ceil(length(expStatesU)/H);
figure; 

if isempty(expStateIT)
    % there is at most one stimulus cond 
    subplot(1,2,1);
    errHistoPolar(dataPhase(StimInd),stimtgtphase,18);
    %hold on; polarhistogram(dataPhase(tgtind), 18); 
    %legend('Actual', 'Target');
    title('Stim Sent');
    subplot(1,2,2);
    errHistoPolar(dataPhase(StimTrainRec),stimtgtphase,18);
    %hold on; polarhistogram(dataPhase(tgtind), 18); 
    %legend('Actual', 'Target');
    title('Stim Recd');
else

for ESi = 1:length(expStatesU)
    expState = expStatesU{ESi};
    ESind = expStateIT{ESi,4};
    %ESind = ESind(t(ESind) > datetime(2024,08,29,17,08,11));
    subplot(H*2,W,ESi);
    errHistoPolar(dataPhase(ESind),stimtgtphase,18); 
    %hold on; polarhistogram(dataPhase(tgtind), 18); 
    %legend('Actual', 'Target');
    title(['Stim Sent during ',expState])
    ESind = expStateIT{ESi,5};
    subplot(H*2,W,ESi + H*W);
    errHistoPolar(dataPhase(ESind),stimtgtphase,18); 
    %hold on; polarhistogram(dataPhase(tgtind), 18); 
    %legend('Actual', 'Target');
    title(['Stim Recd during ',expState])
end

end

%% inter-stim interval

if isempty(expStateIT)
    numStimPeriods = 2; % assume 1 peak, 1 trough, or change this
    stimISI = seconds(diff(stimTimeAbs));
    stimprd = 1; 
    while stimprd < numStimPeriods
        [maxISI,maxISIind] = max(stimISI);
        stimISI = stimISI([1:(maxISIind-1),(maxISIind+1):end]); 
        indBefore = tgtTimeAbs <= stimTimeAbs(maxISIind);
        if maxISIind < length(stimTimeAbs)
            indAfter = tgtTimeAbs >= stimTimeAbs(maxISIind+1);
        else
            indAfter = false(size(indBefore));
        end
        tgtISI = tgtISI(indBefore | indAfter);
        stimprd = stimprd+1;
    end
else
    for ESi = 2:height(expStateIT)
        t1 = expStateIT{ESi-1, 3}; t2 = expStateIT{ESI, 2};
        tgtISI = tgtISI( (tgtTimeAbs >= t2) | (tgtTimeAbs <= t1) );
        stimISI = stimISI( (stimTimeAbs >= t2) | (stimTimeAbs <= t1) );
    end
end

figure; histogram(stimISI); hold on; grid on; histogram(tgtISI);
xlabel('Seconds'); ylabel('Stim Count'); title('Inter-Stim Interval');
legend('Actual', 'Ideal', 'Location','best');

disp(['Actual: Mean ',num2str(mean(stimISI,'omitnan')), ...
    ', SD ',num2str(std(stimISI,'omitnan')),' seconds'])
disp(['Ideal: Mean ',num2str(mean(tgtISI,'omitnan')), ...
    ', SD ',num2str(std(tgtISI,'omitnan')),' seconds'])

%% cycle-by-cycle analysis 

[dataCycle, numCycle] = analyzeByCycle(...
    StimTrainRec, tgtInd, dataPhase, stimtgtphase);

% show missing/extra stim 
numMissing = sum(dataCycle(4,:) < dataCycle(3,:));
numExtra = sum(dataCycle(4,:) > dataCycle(3,:));
figure; 
pie([numMissing, numExtra, numCycle-numMissing-numExtra], ...
    {['Missing: ',num2str(100*numMissing/numCycle,2),'%'], ...
     ['Extra: ',num2str(100*numExtra/numCycle,2),'%'], ...
     ['Correct: ',num2str(100*(1-(numMissing+numExtra)/numCycle),2),'%']});
title('Num. Stimulations by Cycle')

% show timing error 
inan = isnan(dataCycle(1,:)) | isnan(dataCycle(2,:));
terr = tRel(dataCycle(2,~inan))-tRel(dataCycle(1,~inan));
figure; errHisto(terr);
title('Stim time error: actual - target');
xlabel('dur (s)'); ylabel('count'); 
grid on;

%% polar hist blocks of time 
% to do: this should pull in data from notes.txt
%%{
winTimes = datetime(2026,7,16,15,0,0) + [...
    minutes(47.7), minutes(49); ...
    minutes(49), minutes(51)];
winTimes.TimeZone = t0.TimeZone;
winNames = {'peak rest', 'trough rest'};
winTgt = [0, pi];

figure; 
for w = 1:height(winTimes)

    % stim sent polar histo
    subplot(4, height(winTimes), w); 
    winInd = (t0+seconds(StimTime) >= winTimes(w,1)) & (t0+seconds(StimTime) < winTimes(w,2));
    errHistoPolar(dataPhase(StimInd(winInd)), winTgt(w), 18);
    title(['Stim Sent - ',winNames{w}]);

    % stim recd polar histo 
    subplot(4, height(winTimes), w+height(winTimes));
    StimIndRec = find(StimTrainRec); StimTimeRec = StimIndRec / SamplingFreq;
    winInd = (t0+seconds(StimTimeRec) >= winTimes(w,1)) & (t0+seconds(StimTimeRec) < winTimes(w,2));
    %winInd = find(winInd)-2;
    errHistoPolar(dataPhase(StimIndRec(winInd)), winTgt(w), 18);
    title(['Stim Recd - ',winNames{w}]);

    % missing/extra 
    subplot(4, height(winTimes), w+2*height(winTimes));
    win = (t >= winTimes(w,1)) & (t < winTimes(w,2));
    dc = analyzeByCycle(StimTrainRec&win, tgtInd&win, dataPhase, winTgt(w));
    numMissing = sum(dc(4,:) < dc(3,:));
    numExtra = sum(dc(4,:) > dc(3,:));
    pie([numMissing, numExtra, numCycle-numMissing-numExtra], ...
        {['Missing: ',num2str(100*numMissing/numCycle,2),'%'], ...
         ['Extra: ',num2str(100*numExtra/numCycle,2),'%'], ...
         ['Correct: ',num2str(100*(1-(numMissing+numExtra)/numCycle),2),'%']});
    title(['Num. Stimulations - ',winNames{w}]);

    % timing histo 
    subplot(4, height(winTimes), w+3*height(winTimes));
    inan = isnan(dc(1,:)) | isnan(dc(2,:));
    terr = tRel(dc(2,~inan))-tRel(dc(1,~inan));
    errHisto(terr);
    title(['Stim time error - ',winNames{w}]);
    xlabel('dur (s)'); ylabel('count'); 
    grid on;

end
%}

%% helper functions 

function errHistoPolar(phAct, phTgt, nbin)
if nargin < 3
    nbin = [];
end
if isempty(nbin)
    nbin = 18; % default
end
x = radfix(phAct-phTgt); % error
x = x*180/pi; % report in degrees 
M = mean(x); % mean
SEM = std(x)/sqrt(length(x)); % standard error
ts = tinv(0.975,length(x)-1); % 95% t statistic
CI = ts*SEM; % 95% CI range
polarhistogram(phAct, nbin);
subtitle(['Err Avg ',num2str(M),'°; 95% CI ±',num2str(CI),'°']);
end

function errHisto(err)
M = mean(err); % mean
SEM = std(err)/sqrt(length(err)); % standard error
ts = tinv(0.975,length(err)-1); % 95% t statistic
CI = ts*SEM; % 95% CI range
histogram(err);
subtitle(['Avg ',num2str(M),'; 95% CI ±',num2str(CI)]);
end


function [dataCycle, numCycle] = analyzeByCycle(...
    StimTrainRec, tgtInd, dataPhase, stimtgtphase)
phCycleStart = radfix(stimtgtphase+pi); % put tgt in middle of cycle
phCycleStart = min(phCycleStart, 3); phCycleStart = max(phCycleStart, -3);
[~,~,iCycleStart] = zerocrossrate(dataPhase - phCycleStart);
iCycleStart = find(iCycleStart);
iCycleStart = iCycleStart(1:2:end);
%iCycleStart = sort(iCycleStart); % necessary?
iCycleEnd = iCycleStart(2:end); iCycleStart = iCycleStart(1:(end-1));
numCycle = length(iCycleStart); 

% cycle-by-cycle data: 
%   1: ind of tgt stim 
%   2: ind of actual stim (nan if none; best if many)
%   3: # tgt stim (should always be 1 unless subthresh?)
%   4: # actual stim 
dataCycle = nan(4, numCycle);
for iCycle = 1:numCycle
    iiCycle = iCycleStart(iCycle):iCycleEnd(iCycle);
    tgtCycle = tgtInd(iiCycle);
    stimCycle = StimTrainRec(iiCycle);
    tgtCycle = iiCycle(1) - 1 + find(tgtCycle);
    stimCycle = iiCycle(1) - 1 + find(stimCycle);
    dataCycle(3,iCycle) = length(tgtCycle);
    dataCycle(4,iCycle) = length(stimCycle);
    if dataCycle(3,iCycle) > 0
        pherr = abs(radfix(dataPhase(tgtCycle) - stimtgtphase));
        if dataCycle(3,iCycle) > 1
            [~,jCycle] = min(pherr);
            tgtCycle = tgtCycle(jCycle);
        end
        dataCycle(1,iCycle) = tgtCycle;
    end
    if dataCycle(4,iCycle) > 0
        pherr = abs(radfix(dataPhase(stimCycle) - stimtgtphase));
        if dataCycle(4,iCycle) > 1
            [~,jCycle] = min(pherr);
            stimCycle = stimCycle(jCycle);
        end
        dataCycle(2,iCycle) = stimCycle;
    end
end
end


function [stateInd, stateStartTime, stateEndTime] = ...
    findExpState(expState, expStates, SrlTimes, tRel)

% find timing of state start and end 
stateOn = strcmp(expStates, expState); % true when on, false when off 
stateChange = [0; diff(stateOn)]; 
stateEnd = stateChange < 0; % true -> false (-1) means end
stateStart = stateChange > 0; % false -> true (+1) means start 
stateStartTime = SrlTimes(stateStart); stateEndTime = SrlTimes(stateEnd);

% match start and end times 
stateStartTime = sort(stateStartTime); stateEndTime = sort(stateEndTime);
timeEnd = nan(size(stateStartTime)); 
for ti = 1:length(stateStartTime)
    ts = stateStartTime(ti); 
    te = stateEndTime >= ts; 
    te = stateEndTime(te);
    if ~isempty(te)
        te = min(te); 
        if ti < length(stateStartTime)
            if te <= stateStartTime(ti+1)
                timeEnd(ti) = te;
            end
        else
            timeEnd(ti) = te;
        end
    end
end
stateStartTime = stateStartTime(~isnan(timeEnd)); 
stateEndTime = timeEnd(~isnan(timeEnd)); 

% get timing of this state 
stateInd = false(size(tRel)); 
for ti = 1:length(stateStartTime)
    ts = stateStartTime(ti); te = stateEndTime(ti);
    stateInd((tRel > ts) & (tRel < te)) = true;
end

end