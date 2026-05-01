%% User selects file 
[fn,fp] = uigetfile('C:\Users\Anderson Lab\Desktop\data\2026-04-30_BenchTopTest\*.ns*');
ns = openNSx(fullfile(fp,fn), 'uV');

[dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS(ns, 65, 64);
dataOneChannelWithArtifact = dataOneChannel; 
t0 = t(1);

stimtgtphase = inputdlg('Stimulation target phase (degrees):');
stimtgtphase = (str2double(stimtgtphase));
stimtgtphase = stimtgtphase*pi/180; % rad

PeakInd = []; TroughInd = []; StimInd = [];

%{ 
%% artifact (outlier) detection 
%{
olStartEnd = [false, isoutlier(diff(dataOneChannel), 'mean')];
olStartEnd = find(olStartEnd); 
olStart = olStartEnd(1:2:end); olEnd = olStartEnd(2:2:end);
%}
dataOneChannel = dataOneChannelWithArtifact;
artExtend = 10; % extend artifact by __ samples 
artBegin = 1; % begin artifact __ samples before stim detection
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
        %dataOneChannel(ind) = myFastForecastAR(ARmdl, dataOneChannel(ind0:(ind-1))', 1); 
        % linear interp (non-causal): 
        dataOneChannel(ind) = nan;
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
sgtitle(fn)
end
%}

%% filter 
myFilt = buildFIRBPF(SamplingFreq,13,30, 8);
dataOneChannel = filtfilt(myFilt,1,dataOneChannel);

%% find ideal stim pattern 
% Get inst. freq. and phase 
[dataPhase, dataFreq] = instPhaseFreq(dataOneChannel, SamplingFreq);

dtgt = radfix(dataPhase-stimtgtphase);
[~,~,tgtInd] = zerocrossrate(dtgt, "TransitionEdge","rising");

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

% label the plot 
legend(lgd)
title(fn)

%% Plot polar histogram 

figure; 
polarhistogram(dataPhase(StimTrainRec),18); 
title({fn,'Stim Received'});
disp(['Mean Error: ',...
    num2str( mean(radfix( dataPhase(StimTrainRec)-stimtgtphase )) ),...
    ' rad']);

%% inter-stim interval

tgtTimeAbs = t(tgtInd);
StimIndRec = find(StimTrainRec); stimTimeAbs = t(StimIndRec);
tgtISI = seconds(diff(tgtTimeAbs));
stimISI = seconds(diff(stimTimeAbs));

tgtTimeAbs = tgtTimeAbs(tgtTimeAbs >= stimTimeAbs(1)); % after stim began
tgtTimeAbs = tgtTimeAbs(tgtTimeAbs <= stimTimeAbs(end)); % before stim ended

figure; histogram(stimISI); hold on; grid on; histogram(tgtISI);
xlabel('Seconds'); ylabel('Stim Count'); title({fn,'Inter-Stim Interval'});
legend('Actual', 'Ideal', 'Location','best');

disp(['Actual: Mean ',num2str(mean(stimISI)),', SD ',num2str(std(stimISI)),' seconds'])
disp(['Ideal: Mean ',num2str(mean(tgtISI)),', SD ',num2str(std(tgtISI)),' seconds'])
