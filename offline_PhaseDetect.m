%% offline_PhaseDetect 
% 
% Toren Arginteanu, Johns Hopkins School of Medicine 
% 
% This script loads recorded brain recording data and performs the phase
% detection algorithm as if it were receiving the data one sample at a time
% in real time. 
% Then, it performs the same phase detection using the full data as a
% ground truth, and compares the results of the simulated "real time"
% against the ground truth. 
% For artifact removal, the script will attempt to determine when the
% actual stimuli in the recorded data occurred, using the trigger channel
% if available and using outlier detection to detect stimulus artifacts.
% Note that for real time purposes, the triggers being scheduled by the
% real-time phase detection algorithm would be used to time the artifact
% removal instead. 
% 
% An overview of the phase detection algorithm is as follows: ...

function [phAll, phEst, frAll, frEst, toStim, phStim] = ...
    offline_PhaseDetect(dataOneChannel, SamplingFreq, StimTrainRec, t, channelName, ...
    PhaseOfInterest, FreqRange, ARwin, ARlen, predWin, artDur, packetLength)

% signal to use default values if any arguments are not passed in 
if nargin < 12
    packetLength = [];
    if nargin < 11
        artDur = [];
        if nargin < 10
            predWin = [];
            if nargin < 9
                ARlen = [];
                if nargin < 8
                    ARwin = [];
                    if nargin < 7
                        FreqRange = [];
                        if nargin < 6
                            PhaseOfInterest = [];
                            if nargin < 5
                                channelName = '';
                                if nargin < 4
                                    t = [];
                                    if nargin < 3
                                        StimTrainRec = [];
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Constants: 


% Simulate phase-dependent stimulation at this phase: 
PhaseOfInterest_default = 0; % radians; i.e. 0 for peak, pi for trough stimulation

% frequency range: 
loco_default = 13; hico_default = 30; % low and high cutoff (Hz); e.g. 13-30 for beta band

% Artifact duration: 
artDur_default = 10; % artifact duration is extended by __ samples 

% AR model parameters: 
ARwin_default = 1000; % #samples of baseline to use to fit the AR model
ARlen_default = 10; % AR model order 

% AR-model prediction duration: 
predWin_default = 20; % #samples ahead to predict at each time step 

% length of packets of incomming data 
packetLength_default = 1; % samples; e.g. 1 for each sample entered individually 


% apply default values if necessary 
if isempty(FreqRange)
    hico = hico_default; loco = loco_default;
else
    hico = FreqRange(2); loco = FreqRange(1);
end
varnames = ["PhaseOfInterest", "artDur", "ARwin", "ARlen", "predWin", "packetLength"];
for v = varnames
    if isempty(eval(v))
        eval(v+" = "+v+"_default;");
    end
end
clear v varnames

%% Part A: Setup 
% Load and preprocess brain recording data from a file. 
% Setup the signals, filter, and AR model that will be used to simulate the
% real-time process. 

if nargin < 1
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
        channelName, channelIndex, channelIndexStim, channelNames]...
        = getRecordedData_NS();
else
    if isempty(StimTrainRec)
        StimTrainRec = false(size(dataOneChannel));
    end
    if isempty(t)
        t = 1:length(dataOneChannel); t = (t-1)/SamplingFreq;
    end
    if isempty(channelName)
        channelName = ' ';
    end
end

dataOneChannel = dataOneChannel - mean(dataOneChannel);
dataOneChannelWithArtifact = dataOneChannel; 

% Get indexes of stimulus: 
% find when artifacts are believed to occur 
isOut = isoutlier(dataOneChannel, 'mean');
isArt = isOut | StimTrainRec; 
if artDur > 0
    isArt = movsum(isArt, artDur) > 0;
end

%% A.2 Identify baseline and fit AR model
% In real time, this would be determined by the research team at some point
% during the procedure, preferably when the electrodes have been inserted
% into a stable location but the stimulus has not yet begun. 

% Find a clean baseline to train the AR model. 
artIndAll = isArt;
artIndAll(1) = true; artIndAll(end) = true;
artIndAll = find(artIndAll);
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); 
baselineStartInd = artIndAll(baselineStartInd); 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
baselineWin = (baselineEndInd-baselineStartInd) + [-3,3]*ARwin; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); 
baselineWin(2) = min(length(dataBaseline),baselineWin(2));

% Train the AR model on unfiltered data. 
dataBaseline1 = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl_unfilt = ar(iddata(dataBaseline1', [], 1/SamplingFreq), ARlen, 'yw');

%% A.3 Band-Pass Filtering setup 
% Get FIR filter weights, filter the signal, and train another AR model on
% filtered data. 

% filtering bound rules 
minfac         = 2;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length

% filter order 
if loco>0
    filtord = minfac*fix(SamplingFreq/loco);
elseif hico>0
    filtord = minfac*fix(SamplingFreq/hico);
end
if filtord < min_filtorder
    filtord = min_filtorder;
end

filtwts = fir1(filtord, [loco, hico]./(SamplingFreq/2));
dataBaseline = filtfilt(filtwts,1,dataBaseline);
%[dataBaseline, filtwts] = Myeegfilt(dataBaseline,SamplingFreq,loco,hico);
filtord = length(filtwts);
filtinit = zeros(filtord-1,1); % FIR filter Initial Condition
filtdelay = ceil(filtord/2); % delay (#samples) caused by FIR filter
dataBaseline1 = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl_filt = ar(iddata(dataBaseline1', [], 1/SamplingFreq), ARlen, 'yw');

%% Part B: Real-Time Algorithm 
% Loop through each sample of data and perform the algorithm as if the data
% were being received in real time. 
toStim = false(size(dataOneChannel)); % intended stimulus trigger pulses 
phEst = nan(size(dataOneChannel)); % instantaneous phase estimate (rad)
frEst = nan(size(dataOneChannel)); % instantaneous frequency estimate (Hz)
i2nextStim_prev = Inf; % #samples to next stim pulse 
dataOneChannelFilt = zeros(size(dataOneChannel));

progTick = .05; prog = 0; % track progress

for tind = packetLength:packetLength:length(dataOneChannel)
    tinds = (tind-packetLength+1):tind;

    % track progress
    prog = prog + packetLength/length(dataOneChannel);
    if prog > progTick
        prog = prog - progTick; 
        disp(['Progress: ',num2str(100*tind/length(dataOneChannel)),'%']);
    end

    % Step 1: Remove Artifact 
    % Replace artifact-corrupted signal with AR model-forecasted data.  
    if artDur > 0
        for ind1 = tinds(isArt(tinds))
            ind0 = ind1 - ARlen;
            if ind0 > 0
                dataOneChannel(ind1) = myFastForecastAR(ARmdl_unfilt, ...
                    dataOneChannel(ind0:(ind1-1))', 1);
            end
        end
    end

    % Step 2: Filter 
    % Extract data within desired frequency band, e.g. beta
    [dataOneChannelFilt(tinds),filtinit] = filter(filtwts,1, ...
        dataOneChannel(tinds), filtinit);

    ind0 = tinds(1) - ARwin;
    if ind0 > 0
        % Step 3: AR-based forecasting 
        % Predict some duration ahead using the AR model; this will be used
        % to pad the Hilbert transform
        dataPast = dataOneChannelFilt(ind0:(tind-1))';
        dataFuture = myFastForecastAR(ARmdl_filt, dataPast, predWin);

        % Step 4: Phase and Frequency Estimation
        % Using the Hilbert transform, estimate the current instantaneous
        % phase and frequency, and use this to calculate when the next
        % PhaseOfInterest will likely occur. 
        [t2nextStim,i2nextStim, phEst(tind),frEst(tind)] = blockPDS(...
            dataPast, dataFuture, SamplingFreq, PhaseOfInterest, ...
            filtdelay/SamplingFreq, ... FIR filter imposes this delay (s)
            loco, hico);

        % Step 5: Send a stimulus pulse when appropriate 
        i2nextStim_prev = i2nextStim_prev - packetLength; % one packet passed
        if i2nextStim_prev < packetLength
            toStim(tinds(1)-1+i2nextStim_prev) = true;
        end
        i2nextStim_prev = i2nextStim - filtdelay;
    end
end

% re-align timing; correct for filter delay 
dataOneChannelFilt = [dataOneChannelFilt(filtdelay:end), zeros(1,filtdelay-1)];
phEst = [phEst(filtdelay:end), nan(1,filtdelay-1)]; 
frEst = [frEst(filtdelay:end), nan(1,filtdelay-1)]; 

%% Part C: Evaluate Real-Time results 
% Compare simulated real-time output with offline-computed ground truth 

% plot artifact removal 
figure; 
ax(1) = subplot(211); 
plot(t, dataOneChannelWithArtifact, 'k'); 
grid on; hold on; 
plot(t, dataOneChannel, 'b'); 
hold on; plot(t, dataOneChannel.*(isArt), '--r');
if artDur > 0
    title('Artifact Removal'); 
else
    title('Artifact Removal Disabled');
end
ylabel(channelName);
ax(2) = subplot(212); 
plot(t, StimTrainRec); 
ylabel('ainp1');
grid on; linkaxes(ax, 'x'); 

% compare instantaneous phase, frequency: estimated vs actual
%dataOneChannelFilt2 = Myeegfilt(dataOneChannel,SamplingFreq,loco,hico);
dataOneChannelFilt2 = filtfilt(filtwts,1,dataOneChannel);
[phAll, frAll] = instPhaseFreq(dataOneChannelFilt2, SamplingFreq);
frAll = min(frAll, hico); frAll = max(frAll, loco);
%phAll2 = sin(phAll); phEst2 = sin(phEst);
phErr = phEst - phAll; % [-2*pi -> 2*pi];
%phErr = mod(phErr + 2*pi, 2*pi); % [0 -> 2*pi];
%phErr = phErr - 2*pi*(phErr >= pi); % [-pi -> pi]
frErr = frEst - frAll;
figure; 
subplot(2,2,1); plot(phAll, phEst, '.'); 
grid on; title('Phase Accuracy'); 
xlabel('Offline Calc. Phase (rad)'); ylabel('RealTime Pred. Phase (rad)'); 
xlim([-pi pi]); ylim([-pi pi]);
subplot(2,2,3); polarhistogram(phErr); title('Phase Error (RealTime-Offline)');
subplot(2,2,2); plot(frAll, frEst, '.'); 
xlim([loco hico]); ylim([loco hico]);
grid on; title('Frequency Accuracy'); 
xlabel('Offline Calc. Freq. (Hz)'); ylabel('RealTime Pred. Freq. (Hz)'); 
subplot(2,2,4); histogram(frErr); title('Freq. Error (RealTime-Offline)');
ylabel('Count'); xlabel('Freq. Error (Hz)'); grid on;

% show time-series of intended stim 
figure; 
plot(t, dataOneChannelFilt2); hold on; grid on; 
stem(t(toStim), dataOneChannelFilt2(toStim));
title('Stimulus Timing');
legend('Data', 'Stim', 'Location','westoutside')

% show true phase of intended stim
figure; sgtitle(['Goal = ',num2str(PhaseOfInterest*180/pi),' degrees'])
subplot(2,2,1); polarhistogram(phAll(toStim), 18); 
title('Actual Phase of Stim'); 
subplot(2,2,2); polarhistogram(phEst(toStim), 18);
title('Est. Phase of Stim');
subplot(2,2,3); polarhistogram(phAll(StimTrainRec), 18);
title('Actual Phase of Recorded Stim');

% process function outputs 
phStim = phAll(toStim);
phAll = phAll(~isnan(phEst)); phEst = phEst(~isnan(phEst));
frAll = frAll(~isnan(frEst)); frEst = frEst(~isnan(frEst));

end
