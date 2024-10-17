%% access files 
[fn,fp] = uigetfile('*.mat');
load([fp,filesep,fn])
ns = NS5;

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

%% filter 
dataOneChannel = Myeegfilt(dataOneChannel,SamplingFreq,13,30);

%% fit and save
dta = iddata(dataOneChannel',[],1/SamplingFreq);
ARmdl = ar(dta, 10, 'yw');
save([fp,filesep,fn,'_ARmdl.mat'], "ARmdl", "dta", "SamplingFreq")