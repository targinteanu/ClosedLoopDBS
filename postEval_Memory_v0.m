%% User selects folder; MATLAB loads all files 
filepath00 = uigetdir('Saved Data Memory'); 
OnlineFiles = dir([filepath00,filesep,'OnlineDisplaySavedData*.mat']);
OnlineFile = OnlineFiles(1); 
NSFiles = dir([filepath00,filesep,'*.ns*']); 
cd00 = cd; cd(filepath00); 
load(OnlineFile.name); 

% interpret first NS file (should be lower sample rate, e.g. ns2)
NSFile = NSFiles(1); 
ns = openNSx(NSFile.name, 'uV'); 
cd(cd00); 
[dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames, packetLoss]...
    = getRecordedData_NS(ns);

% interpret second NS file if there is one 
if length(NSFiles) > 1
    cd(filepath00); ns5 = openNSx(NSFiles(2).name, 'uV');
    cd(cd00);
    [dataOneChannel2, StimTrainRec2, dataAllChannels2, SamplingFreq2, t2, tRel2, ...
        channelName2, channelIndex2, channelIndexStim2, channelNames2, packetLoss2]...
        = getRecordedData_NS(ns5);
    tNum = seconds(t - t(1)); tNum2 = seconds(t2 - t(1)); % to numeric

    if isempty(dataOneChannel) && ~isempty(dataOneChannel2)
        % resample everything to SamplingFreq2 
        %dataAllChannels = resample(dataAllChannels', round(SamplingFreq2), round(SamplingFreq))';
        dataAllChannels = interp1(tNum,dataAllChannels',tNum2,'nearest','extrap')';
        %StimTrainRec = resample(StimTrainRec, round(SamplingFreq2), round(SamplingFreq));
        if SamplingFreq2 < SamplingFreq
            StimTrainRec = movmean(single(StimTrainRec), ceil(SamplingFreq/SamplingFreq2));
        end
        StimTrainRec = interp1(tNum,single(StimTrainRec),tNum2,'linear','extrap');
        StimTrainRec = StimTrainRec > 0;
        dataAllChannels = [dataAllChannels2; dataAllChannels];
        SamplingFreq = SamplingFreq2; t = t2; tRel = tRel2;
        channelName = channelName2; channelIndex = channelIndex2; 
        StimTrainRec = StimTrainRec | StimTrainRec2;
        channelIndexStim = [channelIndexStim+length(channelNames2), channelIndexStim2];
        dataOneChannel = dataOneChannel2;
        channelNames = [channelNames2, channelNames];

    else
        % resample everything to SamplingFreq 
        %dataAllChannels2 = resample(dataAllChannels2', round(SamplingFreq), round(SamplingFreq2))';
        dataAllChannels2 = interp1(tNum2,dataAllChannels2',tNum,'nearest','extrap')';
        dataOneChannel2 = dataAllChannels2(channelIndex2,:);
        %StimTrainRec2 = resample(StimTrainRec2, round(SamplingFreq), round(SamplingFreq2));
        if SamplingFreq < SamplingFreq2
            StimTrainRec2 = movmean(single(StimTrainRec2), ceil(SamplingFreq2/SamplingFreq));
        end
        StimTrainRec2 = interp1(tNum2,single(StimTrainRec2),tNum,'linear','extrap');
        StimTrainRec2 = StimTrainRec2 > 0;
        dataAllChannels = [dataAllChannels; dataAllChannels2];
        StimTrainRec = StimTrainRec | StimTrainRec2;
        channelIndexStim = [channelIndexStim2+length(channelNames), channelIndexStim];
        channelNames = [channelNames, channelNames2];
    end
else
    packetLoss2 = false;
end

% handle any packet loss 
if packetLoss || packetLoss2
    tRelInt = tRel(1):(1/SamplingFreq):tRel(end);
    dataAllChannels = interp1(tRel,dataAllChannels',tRelInt,"nearest","extrap")';
    dataOneChannel = dataAllChannels(channelIndex,:);
    StimTrainRec = interp1(tRel,single(StimTrainRec),tRelInt,"linear","extrap");
    StimTrainRec = StimTrainRec > 0;
    t = seconds(tRelInt) + t(1); tRel = tRelInt;
end

t0 = t(1);
dataOneChannelWithArtifact = dataOneChannel; 
if numel(channelIndexStim)
    channelIndexStim = channelIndexStim(1);
end

%% Get indexes of peaks, troughs, and stimulus pulses 
PeakInd = PeakTime*SamplingFreq; 
TroughInd = TroughTime*SamplingFreq; 
StimInd = StimTime*SamplingFreq; 
PeakInd = round(PeakInd); TroughInd = round(TroughInd); StimInd = round(StimInd); 
trimToSize = @(inds, x) inds((inds>=0) & (inds <= length(x)));
PeakInd = trimToSize(PeakInd,t); TroughInd = trimToSize(TroughInd,t); StimInd = trimToSize(StimInd,t);

%% artifact detection
dataOneChannel = dataOneChannelWithArtifact;
artExtend = 25; % extend artifact by __ samples 
if numel(channelIndexStim)
    artIndAll = StimTrainRec; % cerestim trigs
    StimInd = find(StimTrainRec); % replace sent stim with observed
    % no other sources of artifact in the memory protocol
else
    warning('Stimulus trigger channel was not found.')
    % assume there are cerestim trigs, but they are not recorded
    artIndAll = isoutlier(dataOneChannel, 'mean');
    artIndAll(StimInd) = true;
end 
artIndAll = movsum(artIndAll, artExtend) > 0;
artIndAll_PulseTrain = artIndAll;
artIndAll = find(artIndAll);
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); baselineStartInd = artIndAll(baselineStartInd); 

%% set baseline to fit model 
if ~isempty(artIndAll)
baselineWinLen = 1000; ARlen = 10; % samples 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
DCOS = mean(dataBaseline);
dataBaseline = Myeegfilt(dataBaseline,SamplingFreq,4,9, 0, 1024);
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
        dataOneChannel(ind) = myFastForecastAR(ARmdl, dataOneChannel(ind0:(ind-1))', 1);
    end
end

% plot artifact removal 
figure; 
ax(1) = subplot(211); 
plot(t, dataOneChannelWithArtifact-DCOS, 'k', 'LineWidth',1.5); 
grid on; hold on; 
plot(t, dataOneChannel, 'b'); 
title('Artifact Removal'); ylabel(channelName);
ax(2) = subplot(212); 
if numel(channelIndexStim)
    plot(t, dataAllChannels(channelIndexStim,:)); 
    ylabel(channelNames{channelIndexStim});
else
    plot(t, artIndAll_PulseTrain);
    ylabel('outlier?');
end
grid on; linkaxes(ax, 'x'); 
end

%% filter 
myFilt = buildFIRBPF(SamplingFreq,4,9, 8);
dataOneChannel = filtfilt(myFilt,1,dataOneChannel);

%% select time of interest (manually)
% TO DO: make this automatic, pulled from notes.txt ?
% selind = true(size(t)); % no selection
% selind = t >= datetime(2024,10,9,16,35,0) + hours(4);
selind = true(size(t)); % bypass selection - fix this!!!
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