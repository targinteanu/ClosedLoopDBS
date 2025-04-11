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

t.TimeZone = 'America/New_York';
tStim.TimeZone = 'America/New_York';
[t(1) t(end)]
[tStim(1) tStim(end)]

%% output CSV file 
[fn,fp] = uigetfile('*.csv'); 
neuromodulation_output_visualization(fullfile(fp,fn)); % APL internal eval 
tbl = readtable(fullfile(fp,fn));

%% interpret APL data
SamplingFreqAPL = 1000; % Hz
dataAPL = (tbl.data)';
tRelAPL = (tbl.dataTimestamp)'/SamplingFreqAPL; % s
tRelAPL = tRelAPL - tRelAPL(1); % relative time; start at 0
StimIndAPL = (tbl.stimOut > 0);
phaseAPL = (tbl.projPhase)';
t0APL = sscanf(fn, 'neuromod_output_%f-%f-%f_%f-%f-%f*.csv');
t0APL = datetime(t0APL', 'TimeZone','America/New_York');
tAPL = seconds(tRelAPL) + t0APL;

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
StimInd1 = StimTrainRec1; 
figure; 
plot(tStim1, StimTrainRec1); hold on; grid on; 
plot(tStim1(StimInd1), StimTrainRec1(StimInd1), '*r');

%% filter 
% should this be done after alignment instead? 
% should artifact removal be done after alignment, so artifacts can be used
% for alignment? 

% The beta-band filter we are using:
h = [-0.00235822679806431,-0.00222172563171287,-0.00204188612994641,-0.00183038495484977,-0.00160034390923488,-0.00136581269659799,-0.00114119804128303,-0.000940660566584552,-0.000777502797457442,-0.000663572746103515,-0.00060870768315786,-0.000620241860603209,-0.000702600139262313,-0.000856996727262349,-0.00108125463821507,-0.00136975714725822,-0.00171353761034857,-0.00210050869551155,-0.00251582655346347,-0.00294237994222154,-0.00336138903478806,-0.00375309379670696,-0.00409750762578007,-0.00437520858421974,-0.00456813818108741,-0.00466037640254251,-0.00463886162140641,-0.00449402518340973,-0.00422031285552714,-0.00381656887412917,-0.00328626294222335,-0.00263754604641679,-0.00188312720632949,-0.00103997000959327,-0.000128814776282121,0.000826460826970075,0.00179962321217186,0.00276298068623895,0.0036883666800297,0.00454817311860771,0.00531639560056796,0.00596965038000671,0.0064881229614292,0.00685640947577641,0.00706421489728421,0.00710687651962473,0.00698568681140144,0.0067079966366596,0.00628708762755459,0.00574181096389299,0.00509599864561013,0.00437766221518625,0.00361800246321179,0.00285026159970675,0.00210845637912515,0.00142603643810278,0.000834516390501735,0.000362132821790785,3.25780897138549e-05,-0.000136138309429324,-0.000132653120598098,4.70052982061992e-05,0.000398925699844695,0.000911104654665073,0.00156353765180404,0.00232861857064927,0.00317184445611504,0.00405281127112782,0.00492647522528875,0.00574464379723933,0.00645765107504029,0.00701616390251436,0.00737305886804188,0.00748530568041385,0.00731579015430977,0.00683501001270845,0.00602257906226249,0.00486847998420728,0.00337401289754879,0.00155239580123291,-0.000571016280403731,-0.00295891447561934,-0.00556262940660626,-0.0083231287473506,-0.0111723991318278,-0.0140351789751382,-0.0168309960244307,-0.0194764519931052,-0.0218876868709724,-0.0239829478127192,-0.0256851821892519,-0.0269245716658481,-0.0276409241857489,-0.027785843533624,-0.0273246016808124,-0.0262376472292206,-0.0245216937382922,-0.022190344220844,-0.0192742222365224,-0.0158205953431237,-0.0118924926872715,-0.0075673346993105,-0.0029351086641472,0.00190386116090241,0.00684148672355834,0.0117647588792982,0.0165587765132936,0.0211098722437602,0.0253087409838873,0.0290534767491559,0.032252425249128,0.0348267649435504,0.0367127372169914,0.0378634568836062,0.0382502470358433,0.0378634568836062,0.0367127372169914,0.0348267649435504,0.032252425249128,0.0290534767491559,0.0253087409838873,0.0211098722437602,0.0165587765132936,0.0117647588792982,0.00684148672355834,0.00190386116090241,-0.0029351086641472,-0.0075673346993105,-0.0118924926872715,-0.0158205953431237,-0.0192742222365224,-0.022190344220844,-0.0245216937382922,-0.0262376472292206,-0.0273246016808124,-0.027785843533624,-0.0276409241857489,-0.0269245716658481,-0.0256851821892519,-0.0239829478127192,-0.0218876868709724,-0.0194764519931052,-0.0168309960244307,-0.0140351789751382,-0.0111723991318278,-0.0083231287473506,-0.00556262940660626,-0.00295891447561934,-0.000571016280403731,0.00155239580123291,0.00337401289754879,0.00486847998420728,0.00602257906226249,0.00683501001270845,0.00731579015430977,0.00748530568041385,0.00737305886804188,0.00701616390251436,0.00645765107504029,0.00574464379723933,0.00492647522528875,0.00405281127112782,0.00317184445611504,0.00232861857064927,0.00156353765180404,0.000911104654665073,0.000398925699844695,4.70052982061992e-05,-0.000132653120598098,-0.000136138309429324,3.25780897138549e-05,0.000362132821790785,0.000834516390501735,0.00142603643810278,0.00210845637912515,0.00285026159970675,0.00361800246321179,0.00437766221518625,0.00509599864561013,0.00574181096389299,0.00628708762755459,0.0067079966366596,0.00698568681140144,0.00710687651962473,0.00706421489728421,0.00685640947577641,0.0064881229614292,0.00596965038000671,0.00531639560056796,0.00454817311860771,0.0036883666800297,0.00276298068623895,0.00179962321217186,0.000826460826970075,-0.000128814776282121,-0.00103997000959327,-0.00188312720632949,-0.00263754604641679,-0.00328626294222335,-0.00381656887412917,-0.00422031285552714,-0.00449402518340973,-0.00463886162140641,-0.00466037640254251,-0.00456813818108741,-0.00437520858421974,-0.00409750762578007,-0.00375309379670696,-0.00336138903478806,-0.00294237994222154,-0.00251582655346347,-0.00210050869551155,-0.00171353761034857,-0.00136975714725822,-0.00108125463821507,-0.000856996727262349,-0.000702600139262313,-0.000620241860603209,-0.00060870768315786,-0.000663572746103515,-0.000777502797457442,-0.000940660566584552,-0.00114119804128303,-0.00136581269659799,-0.00160034390923488,-0.00183038495484977,-0.00204188612994641,-0.00222172563171287,-0.00235822679806431];

dataOneChannel1 = filtfilt(h,1,dataOneChannel1);
dataAPL1 = filtfilt(h,1,dataAPL1);

%% align APL-ns timing 

tMaxDiff = .1; % max allowed clock diff (minutes) 

% constrict NS time search based on csv file name 
tNSsel = ...
    (t1 > t0APL-minutes(tMaxDiff)) & ...
    (t1 < tAPL(end)+minutes(tMaxDiff));
for varn = {'dataOneChannel1', 'tRel1', 't1', ...
        'tStim1', 'StimInd1', 'StimTrainRec1'}
    v = varn{:};
    eval([v,' = ',v,'(tNSsel);']);
end

% constrict APL time search based on NS start/end time
tAPLsel = ...
    (tAPL > t(1)-minutes(tMaxDiff)) & ...
    (tAPL < t(end)+minutes(tMaxDiff));
dataAPL1 = dataAPL1(tAPLsel); 
StimIndAPL1 = StimIndAPL1(tAPLsel);

% align start time using cross-corr 
[r,l] = xcorr(dataOneChannel1, dataAPL1);
[R,ri] = max(r); L = l(ri);
if L > 0
    L = L+1;
    dataOneChannel1 = dataOneChannel1(L:end); 
    tRel1 = tRel1(L:end); 
    t1 = t1(L:end); 
    tStim1 = tStim1(L:end);
    StimInd1 = StimInd1(L:end);
    StimTrainRec1 = StimTrainRec1(L:end);
elseif L < 0
    L = L-1;
    dataAPL1 = dataAPL1(-L:end); 
    StimIndAPL1 = StimIndAPL1(-L:end);
end

% align end time by truncating whichever is now longer 
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