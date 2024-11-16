%% User selects folder; MATLAB loads all files 
filepath = uigetdir('Saved Data Memory'); 
OnlineFiles = dir([filepath,filesep,'OnlineDisplaySavedData*.mat']);
OnlineFile = OnlineFiles(1); 
NSFiles = dir([filepath,filesep,'*.ns*']); 
NSFile = NSFiles(1); 
cd00 = cd; cd(filepath); 
load(OnlineFile.name); ns = openNSx(NSFile.name); 
cd(cd00); 

%% User selects channel 
channelNames = {ns.ElectrodesInfo.Label};
channelIndex = listdlg('ListString', channelNames);
channelName = channelNames{channelIndex};
channelIndexStim = find(contains(channelNames, 'ainp1'));

%% interpret data from ns structure 
SamplingFreq = ns.MetaTags.SamplingFreq;
dataAllChannels = double(ns.Data); 
dataOneChannel = dataAllChannels(channelIndex,:);
dataOneChannelWithArtifact = dataOneChannel;

%% Get timing data
tRel = linspace(0,ns.MetaTags.DataPointsSec,ns.MetaTags.DataPoints);
t = seconds(tRel);
t0 = datetime(ns.MetaTags.DateTime); 
t = t+t0; 

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

%% determine encode/decode phases of experiment 
expStates = {SerialLog.ParadigmPhase}';
SrlTimes = [SerialLog.TimeStamp]';
[indEncode, encodeStart, encodeEnd] = ...
    findExpState('ENCODE', expStates, SrlTimes, tRel);
[indDecode, decodeStart, decodeEnd] = ...
    findExpState('DECODE', expStates, SrlTimes, tRel);
indNeither = (~indEncode)&(~indDecode);

% convert times to absolute 
varnames = {'encodeStart', 'encodeEnd', 'decodeStart', 'decodeEnd'};
for V = varnames
    v = V{:};
    eval([v,' = seconds(',v,') + t0;'])
end

% is each stimulus during encode or decode? 
StimEncodeInd = StimInd(indEncode(StimInd)); 
StimDecodeInd = StimInd(indDecode(StimInd));
StimNeitherInd = StimInd(indNeither(StimInd)); % there should be none 
disp(' ... ')
disp([num2str(length(StimNeitherInd)),' out of ',num2str(length(StimInd)), ...
    ' stimuli during neither encode nor decode.'])

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
patch([encodeStart, encodeEnd, encodeEnd, encodeStart]', ...
    repmat([dataMax; dataMax; dataMin; dataMin],1,length(encodeStart)), ...
    'c', 'FaceAlpha', .2); 
patch([decodeStart, decodeEnd, decodeEnd, decodeStart]', ...
    repmat([dataMax; dataMax; dataMin; dataMin],1,length(decodeStart)), ...
    'g', 'FaceAlpha', .2); 

% label the plot 
legend('Data', 'Peaks', 'Troughs', 'Stimulus', 'Encode', 'Decode')

%% Get inst. freq. and phase 
[dataPhase, dataFreq] = instPhaseFreq(dataOneChannel, SamplingFreq);

%% Plot polar histogram 
figure; 
subplot(221); polarhistogram(dataPhase(PeakInd),18); 
title('Predicted Peaks');
subplot(222); polarhistogram(dataPhase(TroughInd),18); 
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