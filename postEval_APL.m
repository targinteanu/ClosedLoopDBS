%% access files 

% i.e. if sitm on ns5 and rec on ns2
yn = questdlg('Are rec and stim on different files?');
if strcmp(yn, 'Yes')
    % open two ns files 
    [dataOneChannel, ~, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, ~, channelNames]...
    = getRecordedData_NS();
    [~, StimTrainRec, dataAllChannelsStim, SamplingFreqStim, tStim, tRelStim, ...
    channelNameStim, ~, channelIndexStim, channelNamesStim]...
    = getRecordedData_NS();
    resampleStim = true;
elseif strcmp(yn, 'No')
    % open one ns file
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS();
    for varn = {'dataAllChannels', 'SamplingFreq', 't', 'tRel', 'channelIndex', 'channelNames'}
        v = varn{:};
        eval([v,'Stim = ',v,';']);
    end
    channelNameStim = 'ainp1';
else
    % end/error here ?
end

t(1)
tStim(1)

%% output CSV file 
[fn,fp] = uigetfile('*.csv'); 
neuromodulation_output_visualization(fullfile(fp,fn)); % APL internal eval 
tbl = readtable(fullfile(fp,fn));

%% interpret APL data
SamplingFreqAPL = 1000; % Hz
dataAPL = (tbl.data)';
tAPL = (tbl.dataTimestamp)'/SamplingFreqAPL; % s
StimIndAPL = (tbl.stimOut > 0);
phaseAPL = (tbl.projPhase)';

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
    tRel1 = resample(tRel1,SamplingFreqAPL,SamplingFreq);
    t1 = seconds(tRel1) + t0;
end

%% get stim indexes 
StimInd1 = diff(StimTrainRec1); 
StimInd1 = max(0, StimInd1);
StimInd1 = isoutlier(StimInd1, 'quartiles');
StimInd1 = [StimInd1, false];
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
elseif L < 0
    dataAPL1 = dataAPL1(-L:end); 
    StimIndAPL1 = StimIndAPL1(-L:end);
end
minlen = min(length(dataAPL1), length(dataOneChannel1));
dataOneChannel1 = dataOneChannel1(1:minlen);
tRel1 = tRel1(1:minlen);
t1 = t1(1:minlen);
dataAPL1 = dataAPL1(1:minlen); 
StimIndAPL1 = StimIndAPL1(1:minlen);

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