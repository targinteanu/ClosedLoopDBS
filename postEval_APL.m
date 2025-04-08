%% access files 

% i.e. if sitm on ns5 and rec on ns2
yn = questdlg('Are rec and stim on different files?');
if strcmp(yn, 'Yes')
    % open two ns files 
    [dataOneChannel, ~, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, ~, channelNames]...
    = getRecordedData_NS();
    [StimTrainRec, ~, dataAllChannelsStim, SamplingFreqStim, tStim, tRelStim, ...
    channelNameStim, ~, channelIndexStim, channelNamesStim]...
    = getRecordedData_NS();
    resampleStim = true;
elseif strcmp(yn, 'No')
    % open one ns file
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS();
    channelNameStim = channelNames{channelIndexStim};
    for varn = {'dataAllChannels', 'SamplingFreq', 't', 'tRel', 'channelIndex', 'channelNames'}
        v = varn{:};
        eval([v,'Stim = ',v,';']);
    end 
    resampleStim = false;
else
    % end/error here ?
end
dataOneChannelWithArtifact = dataOneChannel;
t0 = t(1);

[t(1) t(end)]
[tStim(1) tStim(end)]

%% output CSV file 
[fn,fp] = uigetfile('*.csv'); 
neuromodulation_output_visualization(fullfile(fp,fn)); % APL internal eval 
tbl = readtable(fullfile(fp,fn));

%% interpret APL data
SamplingFreqAPL = 950; % Hz
dataAPL = (tbl.data)';
tAPL = (tbl.dataTimestamp)'/SamplingFreqAPL; % s
StimIndAPL = (tbl.stimOut > 0);
phaseAPL = (tbl.projPhase)';

%% artifact removal 

% detect artifacts 
artExtend = 10; % extend artifact by __ samples 
artIndAll = StimTrainRec; 
artIndAll = movsum(artIndAll, artExtend) > 0;
artIndAll(1) = true; artIndAll(end) = true;
artIndAll = find(artIndAll);
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); baselineStartInd = artIndAll(baselineStartInd); 

% set baseline to fit model 
if ~isempty(artIndAll)
baselineWinLen = 1000; ARlen = 10; % samples 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
dataBaseline = Myeegfilt(dataBaseline,SamplingFreq,13,30);
baselineWin = (baselineEndInd-baselineStartInd) + [-1,1]*baselineWinLen; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); baselineWin(2) = min(length(dataBaseline),baselineWin(2));
dataBaseline = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl = ar(iddata(dataBaseline', [], 1/SamplingFreq), ARlen, 'yw');

% remove artifact 
dataOneChannel = dataOneChannelWithArtifact;
dataOneChannel = dataOneChannel - mean(dataOneChannel); % correct DC offset

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

%% resample 
dataAPL1 = dataAPL; dataOneChannel1 = dataOneChannel;
StimIndAPL1 = StimIndAPL;
tRel1 = tRel; t1 = t;
tStim1 = tStim; tRelStim1 = tRelStim; StimTrainRec1 = StimTrainRec;
if resampleStim
    % resample stim to match rec 
    StimTrainRec1 = interp1(tStim1, StimTrainRec1, t1);
    tStim1 = t1; 
end
if SamplingFreqAPL ~= SamplingFreq
    % resample ns to match APL
    dataOneChannel1 = resample(dataOneChannel1,SamplingFreqAPL,SamplingFreq);
    tRel1 = ((1:length(dataOneChannel1)) - 1) / SamplingFreq;
    t1 = seconds(tRel1) + t0;
end

%% get stim indexes 
StimInd1 = StimTrainRec1; 
figure; 
plot(tStim1, StimTrainRec1); hold on; grid on; 
plot(tStim1(StimInd1), StimTrainRec1(StimInd1), '*r');

%% filter 
dataOneChannel1 = Myeegfilt(dataOneChannel1,SamplingFreqAPL,13,30);
dataAPL1 = Myeegfilt(dataAPL1,SamplingFreqAPL,13,30);

%% align APL-ns timing 
% assumes APL data is shorter-duration than blackrock recording 
if length(dataAPL1) > length(dataOneChannel1)
    warning('APL data is longer duration than recording; this might not work properly.')
end
[r,l] = xcorr(dataOneChannel1, dataAPL1);
[R,ri] = max(r); L = l(ri);
if L > 0
    dataOneChannel1 = dataOneChannel1(L:end); 
    tRel1 = tRel1(L:end); 
    t1 = t1(L:end); 
    tStim1 = tStim1(L:end);
    StimInd1 = StimInd1(L:end);
    StimTrainRec1 = StimTrainRec1(L:end);
elseif L < 0
    dataAPL1 = dataAPL1(-L:end); 
    StimIndAPL1 = StimIndAPL1(-L:end);
end
minlen = min(length(dataAPL1), length(dataOneChannel1));
for varn = {'dataOneChannel1', 'tRel1', 't1', 'dataAPL1', 'StimIndAPL1', ...
        'tStim1', 'StimInd1', 'StimTrainRec1'}
    v = varn{:};
    eval([v,' = ',v,'(1:minlen);']);
end

figure; 
subplot(3,1,1); plot(l,r); grid on; hold on; plot(L,R,'o'); 
xlabel('lag'); ylabel('corr'); 
subplot(3,1,2); plot(t1, dataOneChannel1); grid on; hold on; 
plot(t1, dataAPL1); 
ylabel('aligned data'); legend('BlackRock', 'APL');
subplot(3,1,3); plot(dataOneChannel1, dataAPL1, '.'); grid on; 
xlabel('BlackRock'); ylabel('APL');

%% plot time series 
figure; plot(t1, dataOneChannel1); grid on; hold on; 
plot(t1, dataAPL1);
plot(t1(StimIndAPL1), dataAPL1(StimIndAPL1), '^m'); 
plot(tStim1(StimInd1), dataOneChannel1(StimInd1), '*r');
legend('Data - BlackRock', 'Data - APL', 'Intended Stim', 'Received Stim')

%% plot polar histogram 
[dataPhase, dataFreq] = instPhaseFreq(dataOneChannel1, SamplingFreqAPL);
[dataPhAPL, dataFrAPL] = instPhaseFreq(dataAPL1, SamplingFreqAPL);
figure; 
subplot(2,2,2); polarhistogram(dataPhase(StimIndAPL1),18); 
title('Intended Stim - BlackRock Data');
subplot(2,2,1); polarhistogram(dataPhAPL(StimIndAPL1),18); 
title('Intended Stim - APL Data');
subplot(2,2,4); polarhistogram(dataPhase(StimInd1),18); 
title('Received Stim - BlackRock Data');
subplot(2,2,3); polarhistogram(dataPhAPL(StimInd1),18); 
title('Received Stim - APL Data');
figure; 
subplot(1,2,1); polarhistogram(dataPhAPL(StimIndAPL1),18); 
title('Intended Stim - APL Data');
subplot(1,2,2); polarhistogram(dataPhase(StimInd1),18); 
title('Received Stim - BlackRock Data');