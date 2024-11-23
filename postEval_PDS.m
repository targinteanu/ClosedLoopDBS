%% access files 

% ns file
[dataOneChannel, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS();

% output file 
[fn,fp] = uigetfile('*SaveFile*.mat');
load([fp,filesep,fn])
PeakTime = PeakTrough(:,1); TroughTime = PeakTrough(:,2);
trimNan = @(x) x(~isnan(x));
PeakTime = trimNan(PeakTime); TroughTime = trimNan(TroughTime);

%% Get indexes of peaks, troughs, and stimulus pulses 
PeakInd = PeakTime*SamplingFreq; 
TroughInd = TroughTime*SamplingFreq; 
PeakInd = round(PeakInd); TroughInd = round(TroughInd);  
trimToSize = @(inds, x) inds((inds>=0) & (inds <= length(x)));
PeakInd = trimToSize(PeakInd,t); TroughInd = trimToSize(TroughInd,t); 

%% filter 
dataOneChannel = Myeegfilt(dataOneChannel,SamplingFreq,13,30);

%% Plot time series 

% plot data and indicate recorded peaks, troughs, stimuli 
figure; plot(t,dataOneChannel); grid on; hold on; 
plot(t(PeakInd),   dataOneChannel(PeakInd),   '^m'); 
plot(t(TroughInd), dataOneChannel(TroughInd), 'vm'); 

% label the plot 
legend({'Data'; 'Peaks'; 'Troughs'})

%% Get inst. freq. and phase 
[dataPhase, dataFreq] = instPhaseFreq(dataOneChannel, SamplingFreq);

%% Plot polar histogram 
figure; 
subplot(121); polarhistogram(dataPhase(PeakInd),18); 
title('Predicted Peaks');
subplot(122); polarhistogram(dataPhase(TroughInd),18); 
title('Predicted Troughs');