%% access files 
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
    openNSx(fullfile(fp,[fn,fe]));
    ns = eval(['NS',fe(end)]);
end

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