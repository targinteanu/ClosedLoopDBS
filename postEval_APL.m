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
    openNSx(fullfile(fp,[fn,fe]));
    ns = eval(['NS',fe(end)]);
end

% output CSV file 
[fn,fp] = uigetfile('*.csv'); 
neuromodulation_output_visualization(fullfile(fp,fn)); % APL internal eval 
tbl = readtable(fullfile(fp,fn));

%% interpret APL data
SamplingFreqAPL = 1000; % Hz
dataAPL = (tbl.data)';
tAPL = (tbl.dataTimestamp)'/SamplingFreqAPL; % s
StimIndAPL = find(tbl.stimOut > 0);

%% User selects channel
channelNames = {ns.ElectrodesInfo.Label}; 
channelIndex = listdlg('ListString', channelNames);
channelName = channelNames{channelIndex};
channelIndexStim = find(contains(channelNames, 'ainp1'));

%% interpret data from ns structure 
SamplingFreq = ns.MetaTags.SamplingFreq;
dataAllChannels = double(ns.Data); 
dataOneChannel = dataAllChannels(channelIndex,:);

%% Get timing data
try
    tRel = linspace(0,ns.MetaTags.DataPointsSec,ns.MetaTags.DataPoints);
catch
    tRel = linspace(0,ns.MetaTags.DataPoints/SamplingFreq,ns.MetaTags.DataPoints);
end
t = seconds(tRel);
t0 = datetime(ns.MetaTags.DateTime); 
t = t+t0; 

%% filter 
dataOneChannel = Myeegfilt(dataOneChannel,SamplingFreq,13,30);
dataAPL = Myeegfilt(dataAPL,SamplingFreqAPL,13,30);

%% align APL-ns timing 
% assumes APL data is shorter-duration than blackrock recording 
dataAPL1 = dataAPL; dataOneChannel1 = dataOneChannel;
StimInd1 = StimIndAPL;
tRel1 = tRel; t1 = t;
if SamplingFreqAPL ~= SamplingFreq
    %{
    % resample APL to match ns
    dataAPL1 = resample(dataAPL1,SamplingFreq,SamplingFreqAPL);
    StimInd1 = (StimInd1-1) * SamplingFreq/SamplingFreqAPL + 1;
    %}
    % resample ns to match APL
    dataOneChannel1 = resample(dataOneChannel1,SamplingFreqAPL,SamplingFreq);
    tRel1 = resample(tRel1,SamplingFreqAPL,SamplingFreq);
    t1 = seconds(tRel1) + t0;
end
if length(dataAPL1) > length(dataOneChannel1)
    warning('APL data is longer duration than recording; this might not work properly.')
end
[r,l] = xcorr(dataOneChannel1, dataAPL1);
[R,ri] = max(r); L = l(ri);
dataOneChannel1 = dataOneChannel1(L:end); 
dataOneChannel1 = dataOneChannel1(1:length(dataAPL1));
tRel1 = tRel1(L:end); tRel1 = tRel1(1:length(dataAPL1));
t1 = t1(L:end); t1 = t1(1:length(dataAPL1));
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
plot(t1(StimInd1), dataOneChannel1(StimInd1), '^m'); 
legend('Data - BlackRock', 'Data - APL', 'Intended Stim')

%% plot polar histogram 
[dataPhase, dataFreq] = instPhaseFreq(dataOneChannel1, SamplingFreqAPL);
[dataPhAPL, dataFrAPL] = instPhaseFreq(dataAPL1, SamplingFreqAPL);
figure; 
subplot(1,2,2); polarhistogram(dataPhase(StimInd1),18); 
title('Intended Stim - BlackRock Data');
subplot(1,2,1); polarhistogram(dataPhAPL(StimInd1),18); 
title('Intended Stim - APL Data');