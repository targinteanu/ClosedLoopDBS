% re-writing the offline PDS algorithm to match more closely with the APL
% algorithm

clear
clc

%% define constants

% Simulate phase-dependent stimulation at this phase: 

PhaseOfInterest = 0; % radians; i.e. 0 for peak, pi for trough stimulation

% frequency range: 

loco = 16; hico = 21; % low and high cutoff (Hz); e.g. 13-30 for beta band

% set packet update rate (every sample, every 5th, every 10th, etc.)

packetLength = 1;

% set length of past data packets (# of samples behind every "processed"
% point)

pastDataLength = 1024; % current APL setting

% AR model parameters: 
ARlen = 10; % AR model order 
predWin = 80; % #samples ahead to predict at each time step

% post stimulus hold
stim_hold_ms = 30;

%% define input

input_type = 'SAR'; % either 'sinusoid','ecog' (load from .NSx file), or 'SAR' (artifact-removed .mat file)

switch input_type
    case 'sinusoid' % creates a 100s, 20 Hz sine wave with an amplitude of 10 units
        fs = 1000;
        t = 0:1/fs:100-(1/fs);
        data = 10*sin(40*pi.*t);

        % define first 5s as baseline
        baselineStartInd = 1;
        baselineEndInd = 5000;
        
    case 'mixed_sinusoid' % creates a 100s sine wave from adding 15, 17, and 20 Hz waves
        fs = 1000;
        t = 0:1/fs:100-(1/fs);
        data = 10*sin(40*pi.*t)+10*sin(30*pi.*t)+10*sin(34*pi.*t);

        % define first 5s as baseline
        baselineStartInd = 1;
        baselineEndInd = 5000;

    case 'ecog' % opens channel 25 of selected file
        openNSx() % PD24N011, side1002.ns2
        fs = NS2.MetaTags.SamplingFreq;
        data = double(NS2.Data(25,:));
        data = data - mean(data);
        t = 1:length(data); t = (t-1)/fs;

        % define a clean 5s baseline window (10-15s for PD24N011)
        baselineStartInd = 10000;
        baselineEndInd = 15000;

    case 'SAR' % current subject = PD25N008
        load('NS2.mat')
        load('Ch39_SAR.mat')
        fs = NS2.MetaTags.SamplingFreq;
        data = Ch39_SAR;
        data = data - mean(data);
        t = 1:length(data); t = (t-1)/fs;

        % define a clean 5s baseline window (18-23s for PD25N008)
        baselineStartInd = 18000;
        baselineEndInd = 23000;
end

%% extract baseline and training data

baseline = data(baselineStartInd:baselineEndInd); 

% bandpass filter baseline data

filtord = 229; % hard coded to match APL algorithm

filtwts = designfilt('bandpassfir', FilterOrder=filtord,...
    CutoffFrequency1 = loco, CutoffFrequency2 = hico,...
    SampleRate = fs);

baselineFilt = filtfilt(filtwts,baseline);

% select portion of baseline to be used to train AR model, default is entire baseline segment
baselineFilt_train = baselineFilt; 

% Train AR model on filtered baseline data
ARmdl_filt = ar(iddata(baselineFilt_train', [], 1/fs), ARlen, 'yw');

%% Simulate real-time phase detection

toStim = false(size(data)); % intended stimulus trigger pulses 

progTick = .05; prog = 0; % track progress

i2t = NaN(size(data)); % samples to next target phase
idxs_to_target = inf;

stim_hold = 0;

for tind = packetLength:packetLength:length(data)-75
   
    % track progress
    prog = prog + packetLength/length(data);
    if prog > progTick
        prog = prog - progTick; 
        disp(['Progress: ',num2str(100*tind/length(data)),'%']);
    end
    
    tinds = (tind-packetLength+1):tind;

    ind0 = tinds(1) - pastDataLength; 

    if ind0 > 0 % make sure enough time has passed to obtain a full packet of data

        % filter incoming data

        dataFilt = filtfilt(filtwts, data(ind0:tinds(end)));

        dataPast = dataFilt(1:end-40)';
        
        % predict future data
        
        dataFuture = myFastForecastAR(ARmdl_filt, dataPast, predWin);

        % extract inst. phase/freq with Hilbert transform
        
        blockData = [dataPast;dataFuture];

        % ----
        
        H_block = hilbert(blockData); 
        phi_block = angle(H_block);

        f_block = gradient(unwrap(phi_block)) *(fs/(2*pi));
        f_block = max(f_block, loco); 
        f_block = min(f_block, hico);

        % plot(f_block)
        % hold on

        % fwinlen = 4; 
        % fwin = -97 + ((-fwinlen):fwinlen); 
        fwin = (-57:-25)-40;
        
        f_inst = mean(f_block(end+fwin));
        phi_inst = phi_block(length(blockData)-65);
     
        % estimate samples to next target phase
        
        rad_to_target = (PhaseOfInterest+(2*pi)) - phi_inst;
        rad_s = (2*pi)*f_inst;
        idxs_to_target = round((rad_to_target/rad_s)*1000)-25;

        % if any(phi_block(25:end)>3) % if there is a target phase predicted soon
        %     idxs_to_target = find(phi_block(25:end)>3);
        % end
        % 
        i2t(tinds) = idxs_to_target;

        if stim_hold > 0
            stim_hold = stim_hold - 1;
        end

        if stim_hold == 0 & idxs_to_target < 25 %&& sum(toStim((tinds(1)+idxs_to_target-30):(tinds(1)+idxs_to_target)))==0
            toStim(tinds(1)+idxs_to_target) = true;
            stim_hold = stim_hold_ms;
        end

    end
end

allDataFilt = filtfilt(filtwts, data);

% stim timing plot
figure()
plot(t,allDataFilt)
hold on
stem(t(toStim(1:length(t))), allDataFilt(toStim(1:length(t))));

phAll = angle(hilbert(allDataFilt));

% stim accuracy histogram
figure()
polarhistogram(phAll(toStim(1:length(t))), 18);