%% User selects folder; MATLAB loads all files 
filepath = uigetdir('Saved Data PD'); 
OnlineFiles = dir([filepath,filesep,'OnlineDisplaySavedData*.mat']);
OnlineFile = OnlineFiles(1); 
NSFiles = dir([filepath,filesep,'*.ns*']); 
NSFile = NSFiles(1); 
cd00 = cd; cd(filepath); 
load(OnlineFile.name); ns = openNSx(NSFile.name); 
cd(cd00); 

[dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS(ns);
dataOneChannelWithArtifact = dataOneChannel; 

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
artExtend = 10; % extend artifact by __ samples 
io = isoutlier(dataOneChannel, 'mean');
artIndAll = io | StimTrainRec; 
artIndAll(StimInd) = true;
artIndAll = movsum(artIndAll, artExtend) > 0;
artIndAll(1) = true; artIndAll(end) = true;
artIndAll = find(artIndAll);
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); baselineStartInd = artIndAll(baselineStartInd); 

%% set baseline to fit model 
if ~isempty(artIndAll)
baselineWinLen = 1000; ARlen = 10; % samples 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
dataBaseline = Myeegfilt(dataBaseline,SamplingFreq,13,30);
baselineWin = (baselineEndInd-baselineStartInd) + [-1,1]*baselineWinLen; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); baselineWin(2) = min(length(dataBaseline),baselineWin(2));
dataBaseline = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl = ar(iddata(dataBaseline', [], 1/SamplingFreq), ARlen, 'yw');

%% remove artifact 
dataOneChannel = dataOneChannelWithArtifact;

for ind = artIndAll
    ind0 = ind - ARlen;
    if ind0 > 0
        dataOneChannel(ind) = myFastForecastAR(ARmdl, dataOneChannel(ind0:(ind-1))', 1);
    end
end

% plot artifact removal 
figure; 
ax(1) = subplot(211); 
plot(t, dataOneChannelWithArtifact, 'k'); 
grid on; hold on; 
plot(t, dataOneChannel, 'b'); 
title('Artifact Removal'); ylabel(channelName);
ax(2) = subplot(212); 
plot(t, dataAllChannels(channelIndexStim,:)); 
ylabel('ainp1');
grid on; linkaxes(ax, 'x'); 
end

%% filter 
dataOneChannel = Myeegfilt(dataOneChannel,SamplingFreq,13,30);

%% determine red/yellow/green phases of experiment 
expStates = {SerialLog.ParadigmPhase}';
SrlTimes = [SerialLog.TimeStamp]';
expStatesU = unique(expStates);
expStateIT = cell(length(expStatesU),4);
for ESi = 1:length(expStatesU)
    expState = expStatesU{ESi};
    [ind, tStart, tEnd] = findExpState(expState, expStates, SrlTimes, tRel);
    indStim = StimInd(ind(StimInd));
    tStart = seconds(tStart) + t0; 
    tEnd   = seconds(tEnd)   + t0;
    expStateIT{ESi,1} = ind; expStateIT{ESi,4} = indStim;
    expStateIT{ESi,2} = tStart; expStateIT{ESi,3} = tEnd;
    clear expState ind tStart tEnd indStim ESi
end

%% Plot time series 

% select the data to plot 
dataMinMax = [min(dataOneChannel), max(dataOneChannel)];
dataMinMax = dataMinMax + [-1,1]*.01*diff(dataMinMax);
dataMin = dataMinMax(1); dataMax = dataMinMax(2); 

% plot data and indicate recorded peaks, troughs, stimuli 
figure; plot(t,dataOneChannel); grid on; hold on; 
plot(t(PeakInd),   dataOneChannel(PeakInd),   '^m'); 
plot(t(TroughInd), dataOneChannel(TroughInd), 'vm'); 
plot(t(StimInd),   dataOneChannel(StimInd),   '*r');

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

% label the plot 
legend([{'Data'; 'Peaks'; 'Troughs'; 'Stim'}; expStatesU])

%% Get inst. freq. and phase 
[dataPhase, dataFreq] = instPhaseFreq(dataOneChannel, SamplingFreq);

%% Plot polar histogram 
figure; 
subplot(121); polarhistogram(dataPhase(PeakInd),18); 
title('Predicted Peaks');
subplot(122); polarhistogram(dataPhase(TroughInd),18); 
title('Predicted Troughs');

%% stim polar histogram
H = 2; W = ceil(length(expStatesU)/H);
figure; 

if isempty(expStateIT)
    % there is at most one stimulus cond 
    polarhistogram(dataPhase(StimInd),18);
    title('Stim');
else

for ESi = 1:length(expStatesU)
    expState = expStatesU{ESi};
    ESind = expStateIT{ESi,4};
    %ESind = ESind(t(ESind) > datetime(2024,08,29,17,08,11));
    subplot(H,W,ESi);
    polarhistogram(dataPhase(ESind),18); 
    title(['Stim during ',expState])
end

end

%% polar hist blocks of time 
% to do: this should pull in data from notes.txt
%{
winTimes = [...
    datetime(2024,10,01,11,15,00), datetime(2024,10,01,11,17,00); ...
    datetime(2024,10,01,11,18,00), datetime(2024,10,01,11,20,00); ...
    datetime(2024,10,01,11,21,00), datetime(2024,10,01,11,23,00); ...
    datetime(2024,10,01,11,24,00), datetime(2024,10,01,11,26,00); ...
    datetime(2024,10,01,11,28,00), datetime(2024,10,01,11,30,00); ...
    datetime(2024,10,01,11,31,00), datetime(2024,10,01,11,33,00); ...
    datetime(2024,10,01,11,34,00), datetime(2024,10,01,11,36,00); ...
    datetime(2024,10,01,11,37,00), datetime(2024,10,01,11,39,00)];
winTimes = winTimes + hours(4); % convert EST to GMT

winNames = [20, 13, 15, 18, 23, 25, 27, 30];

[winNames, winInd] = sort(winNames); winTimes = winTimes(winInd,:); 
winNames = arrayfun(@(x) [num2str(x),' Hz'], winNames, 'UniformOutput',false);

figure;
for w = 1:height(winTimes)
    subplot(2,height(winTimes), w);
    winInd = (t0+seconds(PeakTime) >= winTimes(w,1)) & (t0+seconds(PeakTime) <= winTimes(w,2));
    polarhistogram(dataPhase(PeakInd(winInd)), 18); 
    title(['Predicted Peaks ',winNames{w}]);

    subplot(2,height(winTimes), w+height(winTimes));
    winInd = (t0+seconds(TroughTime) >= winTimes(w,1)) & (t0+seconds(TroughTime) <= winTimes(w,2));
    polarhistogram(dataPhase(TroughInd(winInd)), 18); 
    title(['Predicted Troughs ',winNames{w}]);
end 
%}

%% helper functions 

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