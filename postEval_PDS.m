%% access files 

% ns file
[fn,fp] = uigetfile({'*.ns*'; '*.mat'});
[~,fn,fe] = fileparts(fn);
if strcmpi(fe,'.mat')
    load(fullfile(fp,fn,fe), 'NS2', 'ns2', 'NS5', 'ns5', 'NS', 'ns');
    for vtry = {'NS2', 'ns2', 'NS5', 'ns5', 'NS'}
        if exist(vtry{:})
            ns = eval(vtry{:});
        end
    end
else
    openNSx(fullfile(fp,fn,fe));
    ns = eval(['NS',fe(end)]);
end

% output file 
[fn,fp] = uigetfile('*SaveFile*.mat');
load([fp,filesep,fn])
PeakTime = PeakTrough(:,1); TroughTime = PeakTrough(:,2);
trimNan = @(x) x(~isnan(x));
PeakTime = trimNan(PeakTime); TroughTime = trimNan(TroughTime);

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
try
    tRel = linspace(0,ns.MetaTags.DataPointsSec,ns.MetaTags.DataPoints);
catch
    tRel = linspace(0,ns.MetaTags.DataPoints/SamplingFreq,ns.MetaTags.DataPoints);
end
t = seconds(tRel);
t0 = datetime(ns.MetaTags.DateTime); 
t = t+t0; 

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