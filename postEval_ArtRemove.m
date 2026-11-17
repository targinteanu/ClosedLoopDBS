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
t0 = t(1);

%% reconstruct/align saved data 

% convert
chinfo = {struct('Unit','', 'Resolution',nan, 'SampleRate',SamplingFreq, ...
    'Name',channelName)};
sigSaved = data2timetable({fliplr(double(SignalSaved))}, chinfo);
sigSaved = sigSaved{1};
sigSaved.Time = sigSaved.Time + t0;
dataSaved = sigSaved{:,1};

% retime 
dataOrig = interp1(t, dataOneChannel, sigSaved.Time, "nearest", "extrap");

% rescale 
dataOrig = (dataOrig-mean(dataOrig))/std(dataOrig);
dataSaved = (dataSaved-mean(dataSaved))/std(dataSaved);

%% view two signals and stim
figure; 
plot(sigSaved.Time, dataOrig); grid on; hold on; 
plot(sigSaved.Time, dataSaved, '--'); 
stem(seconds(StimTime)+t0, ones(size(StimTime)), '.');
legend('orig', 'GUI', 'stim');