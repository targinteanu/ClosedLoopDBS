%% User selects folder; MATLAB loads all files 
filepath = uigetdir('Saved Data Memory'); 
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

%% Get indexes of peaks, troughs, and stimulus pulses 
PeakInd = PeakTime*SamplingFreq; 
TroughInd = TroughTime*SamplingFreq; 
StimInd = StimTime*SamplingFreq; 
PeakInd = round(PeakInd); TroughInd = round(TroughInd); StimInd = round(StimInd); 
trimToSize = @(inds, x) inds((inds>=0) & (inds <= length(x)));
PeakInd = trimToSize(PeakInd,t); TroughInd = trimToSize(TroughInd,t); StimInd = trimToSize(StimInd,t);

%% artifact detection
artExtend = 10; % extend artifact by __ samples 
if numel(channelIndexStim)
    artIndAll = dataAllChannels(channelIndexStim,:) > 1e4; % cerestim trigs
    % no other sources of artifact in the memory protocol
else
    warning('Stimulus channel ainp1 was not connected.')
    % assume there are cerestim trigs, but they are not recorded
    artIndAll = isoutlier(dataOneChannel, 'mean');
end 
artIndAll(StimInd) = true;
artIndAll = movsum(artIndAll, artExtend) > 0;
artIndAll_PulseTrain = artIndAll;
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
if numel(channelIndexStim)
    plot(t, dataAllChannels(channelIndexStim,:)); 
    ylabel('ainp1');
else
    plot(t, artIndAll_PulseTrain);
    ylabel('outlier?');
end
grid on; linkaxes(ax, 'x'); 
end

%% filter 
dataOneChannel = Myeegfilt(dataOneChannel,SamplingFreq,13,30);

%% select time of interest (manually)
% TO DO: make this automatic, pulled from notes.txt ?
% selind = true(size(t)); % no selection
selind = t >= datetime(2024,10,9,16,35,0) + hours(4);
tSel = t(selind); tRelSel = tRel(selind);
dataOneChannelSel = dataOneChannel(selind);
selind = find(selind); 
trimToInd = @(inds, x) x((x <= max(inds)) & (x >= min(inds))) - min(inds);
PeakIndSel = trimToInd(selind, PeakInd);
TroughIndSel = trimToInd(selind, TroughInd);
StimIndSel = trimToInd(selind, StimInd);

%% determine encode/decode phases of experiment 
expStates = {SerialLog.ParadigmPhase}';
SrlTimes = [SerialLog.TimeStamp]';
SrlTimesSel = (SrlTimes <= max(tRelSel)) & (SrlTimes >= min(tRelSel));
expStates = expStates(SrlTimesSel);
SrlTimesSel = SrlTimes(SrlTimesSel);
[indEncode, encodeStart, encodeEnd] = ...
    findExpState('ENCODE', expStates, SrlTimesSel, tRelSel);
[indDecode, decodeStart, decodeEnd] = ...
    findExpState('DECODE', expStates, SrlTimesSel, tRelSel);
indNeither = (~indEncode)&(~indDecode);

% convert times to absolute 
varnames = {'encodeStart', 'encodeEnd', 'decodeStart', 'decodeEnd'};
for V = varnames
    v = V{:};
    eval([v,' = seconds(',v,') + t0;'])
end

% is each stimulus during encode or decode? 
StimEncodeInd = StimIndSel(indEncode(StimIndSel)); 
StimDecodeInd = StimIndSel(indDecode(StimIndSel));
StimNeitherInd = StimIndSel(indNeither(StimIndSel)); % there should be none 
disp(' ... ')
disp([num2str(length(StimNeitherInd)),' out of ',num2str(length(StimIndSel)), ...
    ' stimuli during neither encode nor decode.'])

%% Plot time series 

% select the data to plot 
dataMinMax = [min(dataOneChannelSel), max(dataOneChannelSel)];
dataMinMax = dataMinMax + [-1,1]*.01*diff(dataMinMax);
dataMin = dataMinMax(1); dataMax = dataMinMax(2); 

% plot data and indicate recorded peaks, troughs, stimuli 
figure; plot(tSel,dataOneChannelSel); grid on; hold on; lgd = ["data"];
if numel(PeakIndSel)
    plot(tSel(PeakIndSel),   dataOneChannelSel(PeakIndSel),   '^m'); 
    lgd = [lgd, "Peaks"];
end
if numel(TroughIndSel)
    plot(tSel(TroughIndSel), dataOneChannelSel(TroughIndSel), 'vm'); 
    lgd = [lgd, "Troughs"];
end
if numel(StimIndSel)
    plot(tSel(StimIndSel),   dataOneChannelSel(StimIndSel),   '*r');
    lgd = [lgd, "Stimulus"];
end

% shade plot regions indicating encode and decode state of paradigm 
if numel(encodeStart)
patch([encodeStart, encodeEnd, encodeEnd, encodeStart]', ...
    repmat([dataMax; dataMax; dataMin; dataMin],1,length(encodeStart)), ...
    'c', 'FaceAlpha', .2); 
lgd = [lgd, "Encode"];
end
if numel(decodeStart)
patch([decodeStart, decodeEnd, decodeEnd, decodeStart]', ...
    repmat([dataMax; dataMax; dataMin; dataMin],1,length(decodeStart)), ...
    'g', 'FaceAlpha', .2); 
lgd = [lgd, "Decode"];
end

% label the plot 
legend(lgd)

%% Get inst. freq. and phase 
[dataPhase, dataFreq] = instPhaseFreq(dataOneChannelSel, SamplingFreq);

%% Plot polar histogram 
figure; 
subplot(221); polarhistogram(dataPhase(PeakIndSel),18); 
title('Predicted Peaks');
subplot(222); polarhistogram(dataPhase(TroughIndSel),18); 
title('Predicted Troughs');
subplot(223); polarhistogram(dataPhase(StimEncodeInd),18); 
title('Stim during Encode');
subplot(224); polarhistogram(dataPhase(StimDecodeInd),18);
title('Stim during Decode');

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