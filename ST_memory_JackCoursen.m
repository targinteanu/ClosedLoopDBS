channelIndices = {1,2,3,4};

thetaPowerResults = struct();


%% User selects folder; MATLAB loads all files 
filepath = uigetdir('Saved Data Memory'); 
OnlineFiles = dir([filepath,filesep,'OnlineDisplaySavedData*.mat']);
OnlineFile = OnlineFiles(1); 
NSFiles = dir([filepath,filesep,'*.ns*']); 
NSFile = NSFiles(1); 
cd00 = cd; cd(filepath); 
load(OnlineFile.name); ns = openNSx(NSFile.name, 'uV'); 
cd(cd00); 

for idx = 1:length(channelIndices)

% [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
% channelName, channelIndex, channelIndexStim, channelNames]...
% = getRecordedData_NS_JC(ns,cell2mat(channelIndices(idx)));

[dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
channelName, channelIndex, channelIndexStim, channelNames]...
= getRecordedData_NS(fullfile(NSFile.folder, NSFile.name));

dataOneChannelWithArtifact = dataOneChannel; 
t0 = t(1);

%% Get indexes of peaks, troughs, and stimulus pulses 
PeakInd = PeakTime*SamplingFreq; 
TroughInd = TroughTime*SamplingFreq; 
StimInd = StimTime*SamplingFreq; 
PeakInd = round(PeakInd); TroughInd = round(TroughInd); StimInd = round(StimInd); 
trimToSize = @(inds, x) inds((inds>=0) & (inds <= length(x)));
PeakInd = trimToSize(PeakInd,t); TroughInd = trimToSize(TroughInd,t); StimInd = trimToSize(StimInd,t);

%% artifact detection
artExtend = 10; % extend artifact by __ samples 
if numel(channelIndexStim)
    artIndAll = StimTrainRec; % cerestim trigs
    % no other sources of artifact in the memory protocol
else
    warning('Stimulus channel ainp1 was not connected.')
    % assume there are cerestim trigs, but they are not recorded
    artIndAll = isoutlier(dataOneChannel, 'mean');
end 
artIndAll(StimInd) = true;
artIndAll = movsum(artIndAll, artExtend) > 0;
artIndAll_PulseTrain = artIndAll;
artIndAll = find(artIndAll);
[~,baselineStartInd] = max(diff(artIndAll));
baselineEndInd = artIndAll(baselineStartInd+1); baselineStartInd = artIndAll(baselineStartInd); 

%% set baseline to fit model 
if ~isempty(artIndAll)
baselineWinLen = 1000; ARlen = 10; % samples 
dataBaseline = dataOneChannel(baselineStartInd:baselineEndInd); 
dataBaseline = Myeegfilt(dataBaseline,SamplingFreq,4,9);
baselineWin = (baselineEndInd-baselineStartInd) + [-1,1]*baselineWinLen; 
baselineWin = baselineWin/2; baselineWin = round(baselineWin); 
baselineWin(1) = max(1,baselineWin(1)); baselineWin(2) = min(length(dataBaseline),baselineWin(2));
dataBaseline = dataBaseline(baselineWin(1):baselineWin(2));
ARmdl = ar(iddata(dataBaseline', [], 1/SamplingFreq), ARlen, 'yw');

%% remove artifact 
dataOneChannel = dataOneChannelWithArtifact;
dataOneChannel = dataOneChannel - mean(dataOneChannel); % correct DC offset

for ind = artIndAll
    ind0 = ind - ARlen;
    if ind0 > 0
        dataOneChannel(ind) = myFastForecastAR(ARmdl, dataOneChannel(ind0:(ind-1))', 1);
    end
end
end

%% filter 
dataOneChannel = Myeegfilt(dataOneChannel,SamplingFreq,4,9);

%% select time of interest (manually)
% TO DO: make this automatic, pulled from notes.txt ?
% selind = true(size(t)); % no selection
% selind = t >= datetime(2024,10,9,16,35,0) + hours(4);
selind = true(size(t)); % bypass selection - fix this!!!
tSel = t(selind); tRelSel = tRel(selind);
dataOneChannelSel = dataOneChannel(selind);
selind = find(selind); 
trimToInd = @(inds, x) x((x <= max(inds)) & (x >= min(inds))) - min(inds);
PeakIndSel = trimToInd(selind, PeakInd);
TroughIndSel = trimToInd(selind, TroughInd);
StimIndSel = trimToInd(selind, StimInd);

%% determine encode/decode phases of experiment 
expStates = {SerialLog.ParadigmPhase}';
SrlTimes = [SerialLog.TimeStamp]';
SrlTimesSel = (SrlTimes <= max(tRelSel)) & (SrlTimes >= min(tRelSel));
expStates = expStates(SrlTimesSel);
SrlTimesSel = SrlTimes(SrlTimesSel);
[indEncode, encodeStart, encodeEnd] = ...
    findExpState('ENCODE', expStates, SrlTimesSel, tRelSel);
[indDecode, decodeStart, decodeEnd] = ...
    findExpState('DECODE', expStates, SrlTimesSel, tRelSel);
indNeither = (~indEncode)&(~indDecode);

% convert times to absolute 
varnames = {'encodeStart', 'encodeEnd', 'decodeStart', 'decodeEnd'};
for V = varnames
    v = V{:};
    eval([v,' = seconds(',v,') + t0;'])
end

% is each stimulus during encode or decode? 
StimEncodeInd = StimIndSel(indEncode(StimIndSel)); 
StimDecodeInd = StimIndSel(indDecode(StimIndSel));
StimNeitherInd = StimIndSel(indNeither(StimIndSel)); % there should be none 
disp(' ... ')
disp([num2str(length(StimNeitherInd)),' out of ',num2str(length(StimIndSel)), ...
    ' stimuli during neither encode nor decode.'])

%% Plot time series 

% select the data to plot 
dataMinMax = [min(dataOneChannelSel), max(dataOneChannelSel)];
dataMinMax = dataMinMax + [-1,1]*.01*diff(dataMinMax);
dataMin = dataMinMax(1); dataMax = dataMinMax(2); 

% plot data and indicate recorded peaks, troughs, stimuli 
figure; 
plot(tSel,dataOneChannelSel);
grid on; hold on; lgd = ["data"];

% shade plot regions indicating encode and decode state of paradigm 
if numel(encodeStart)
patch([encodeStart, encodeEnd, encodeEnd, encodeStart]', ...
    repmat([dataMax; dataMax; dataMin; dataMin],1,length(encodeStart)), ...
    'c', 'FaceAlpha', .2); 
lgd = [lgd, "Encode"];
end
if numel(decodeStart)
patch([decodeStart, decodeEnd, decodeEnd, decodeStart]', ...
    repmat([dataMax; dataMax; dataMin; dataMin],1,length(decodeStart)), ...
    'g', 'FaceAlpha', .2); 
lgd = [lgd, "Decode"];
end

% label the plot 
legend(lgd)

encodeData = getEncodeData(tSel, dataOneChannelSel, encodeStart, encodeEnd);
decodeData = getDecodeData(tSel, dataOneChannelSel, decodeStart, decodeEnd);

% Compute Theta Power for Encoding and Decoding
[avgThetaPowerEncoding, stdThetaPowerEncoding] = computeThetaPower(encodeData, SamplingFreq);
[avgThetaPowerDecoding, stdThetaPowerDecoding] = computeThetaPower(decodeData, SamplingFreq);

disp(['Avg Theta Power Encoding: ', num2str(avgThetaPowerEncoding), ' ± ', num2str(stdThetaPowerEncoding)]);
disp(['Avg Theta Power Decoding: ', num2str(avgThetaPowerDecoding), ' ± ', num2str(stdThetaPowerDecoding)]);



% Concatenate all encoding segments back-to-back
concatenatedEncodeData = [];
for i = 1:length(encodeData)
    concatenatedEncodeData = [concatenatedEncodeData; encodeData{i}; NaN]; % Add NaN to separate segments
end

% Concatenate all decoding segments back-to-back
concatenatedDecodeData = [];
for i = 1:length(decodeData)
    concatenatedDecodeData = [concatenatedDecodeData; decodeData{i}; NaN]; % Add NaN to separate segments
end

% Find the global min and max for the y-axis range
globalMin = min([concatenatedEncodeData; concatenatedDecodeData], [], 'omitnan');
globalMax = max([concatenatedEncodeData; concatenatedDecodeData], [], 'omitnan');

% Create a figure for encoding and decoding data
figure;
tiledlayout(3,1);

% Plot Encoding Data
nexttile;
plot(concatenatedEncodeData, 'b'); % Blue for encoding
title('Encoding Signal');
ylabel('Amplitude');
xlabel('Time Points');
ylim([globalMin globalMax]); % Set same y-axis range

% Plot Decoding Data
nexttile;
plot(concatenatedDecodeData, 'g'); % Green for decoding
title('Decoding Signal');
ylabel('Amplitude');
xlabel('Time Points');
ylim([globalMin globalMax]); % Set same y-axis range

% Plot Decoding and Encoding Data
nexttile;
plot(concatenatedDecodeData, 'g');
title('Encoding and Decoding Signal (Overlay)');
ylabel('Amplitude');
xlabel('Time Points');
hold on;
plot(concatenatedEncodeData, 'b');

ylim([globalMin globalMax]); % Set same y-axis range

thetaPowerResults(idx).channelName = channelNames{channelIndex}; % Store channel name
thetaPowerResults(idx).encodingPower = avgThetaPowerEncoding;
thetaPowerResults(idx).decodingPower = avgThetaPowerDecoding;
thetaPowerResults(idx).encodingError = stdThetaPowerEncoding;
thetaPowerResults(idx).decodingError = stdThetaPowerDecoding;

end

channelNamesList = {thetaPowerResults.channelName};
encodingPowers = [thetaPowerResults.encodingPower]; 
decodingPowers = [thetaPowerResults.decodingPower]; 
encodingErrors = [thetaPowerResults.encodingError]; 
decodingErrors = [thetaPowerResults.decodingError]; 
%%
figure;
barVals = [encodingPowers; decodingPowers]'; % Grouped bar values
bar(barVals, 'grouped'); % Plot bar chart

% Compute X positions for error bars
hold on;
numGroups = size(barVals, 1);
numBars = size(barVals, 2);
groupWidth = min(0.8, numBars / (numBars + 1.5));
x = nan(numBars, numGroups);

for i = 1:numBars
    x(i,:) = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
end

% Add error bars
errorbar(x(1,:), encodingPowers, encodingErrors, 'k', 'linestyle', 'none', 'linewidth', 1);
errorbar(x(2,:), decodingPowers, decodingErrors, 'k', 'linestyle', 'none', 'linewidth', 1);

% Customize plot
set(gca, 'XTick', 1:length(channelNamesList));
set(gca, 'XTickLabel', channelNamesList);
legend({'Encoding', 'Decoding'});
ylabel('Average Theta Power');
title('Avg Theta Power for Encoding v. Decoding');
hold off;















%% helper functions 


function encodeData = getEncodeData(tSel, dataOneChannelSel, encodeStart, encodeEnd)
    encodeData = cell(1, length(encodeStart)); % Store each segment separately
    for i = 1:length(encodeStart)
        idx = (tSel >= encodeStart(i)) & (tSel <= encodeEnd(i)); 
        segment = dataOneChannelSel(idx);
        if ~isempty(segment)
            encodeData{i} = segment(:); % Store each segment as a separate cell
        end
    end
end

function decodeData = getDecodeData(tSel, dataOneChannelSel, decodeStart, decodeEnd)
    decodeData = cell(1, length(decodeStart));
    for i = 1:length(decodeStart)
        idx = (tSel >= decodeStart(i)) & (tSel <= decodeEnd(i)); 
        segment = dataOneChannelSel(idx);
        if ~isempty(segment)
            decodeData{i} = segment(:);
        end
    end
end


function [avgPower, powerSEM] = computeThetaPower(dataSegments, ~)
    powerVals = []; % Store power values across all data points
    
    % Collect power values from all segments
    for i = 1:length(dataSegments)
        if ~isempty(dataSegments{i})
            powerVals = [powerVals; dataSegments{i}.^2]; % Square amplitude for power
        end
    end
    
    % Compute mean power
    avgPower = mean(powerVals, 'omitnan'); 
    
    % Compute SEM using total number of data points
    N = length(powerVals); % Total number of data points across segments
    if N > 1
        powerSEM = std(powerVals, 'omitnan') / sqrt(N); % Standard Error of the Mean (SEM)
    else
        powerSEM = 0; % If only one data point, SEM is zero
    end
end



function [stateInd, stateStartTime, stateEndTime] = ...
    findExpState(expState, expStates, SrlTimes, tRel)

% find timing of state start and end 
stateOn = strcmp(expStates, expState); % true when on, false when off 
stateChange = [0; diff(stateOn)]; 
stateEnd = stateChange < 0; % true -> false (-1) means end
stateStart = stateChange > 0; % false -> true (+1) means start 
stateStartTime = SrlTimes(stateStart); stateEndTime = SrlTimes(stateEnd);

% match start and end times 
stateStartTime = sort(stateStartTime); stateEndTime = sort(stateEndTime);
timeEnd = nan(size(stateStartTime)); 
for ti = 1:length(stateStartTime)
    ts = stateStartTime(ti); 
    te = stateEndTime >= ts; 
    te = stateEndTime(te);
    if ~isempty(te)
        te = min(te); 
        if ti < length(stateStartTime)
            if te <= stateStartTime(ti+1)
                timeEnd(ti) = te;
            end
        else
            timeEnd(ti) = te;
        end
    end
end
stateStartTime = stateStartTime(~isnan(timeEnd)); 
stateEndTime = timeEnd(~isnan(timeEnd)); 

% get timing of this state 
stateInd = false(size(tRel)); 
for ti = 1:length(stateStartTime)
    ts = stateStartTime(ti); te = stateEndTime(ti);
    stateInd((tRel > ts) & (tRel < te)) = true;
end

end
