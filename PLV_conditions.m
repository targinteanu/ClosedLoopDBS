% channelIndices = 1:63;
% numChannels = length(channelIndices);
% phaseData = cell(1, numChannels);  % Store phase data for each channel
% betaData = cell(1, numChannels);
% 
% 
% [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
%     channelName, channelIndex, channelIndexStim, channelNames]...
%     = getRecordedData_NS('datafile001.ns2',1);
% 
% for idx = 1:numChannels
%     ch_raw = dataAllChannels(idx, :); % Raw signal
%     ch_filtered = Myeegfilt(ch_raw, SamplingFreq, 13, 30, 0, 1024); % Filter to beta band (13-30 Hz)
%     betaData{idx} = ch_filtered;
%     [ch_phase, ~] = instPhaseFreq(ch_filtered, SamplingFreq); % Compute phase
%     phaseData{idx} = ch_phase; % Store phase data
%     disp(idx)
% end
% disp('Finished storing phases')



%BETA POWER
numChannels = 63;
avgPower_1 = zeros(numChannels, 1);
avgPower_2 = zeros(numChannels, 1);
avgPower_3 = zeros(numChannels, 1);
avgPower_4 = zeros(numChannels, 1);


SamplingFreq = 1000;  % Hz


for ch = 1:63
    powerVals = [];
    for i = 1:size(intervals_1_sec, 1)
        s = round(intervals_1_sec(i, 1) * SamplingFreq);
        e = round(intervals_1_sec(i, 2) * SamplingFreq);
        s = max(s, 1);
        e = min(e, length(betaData{ch}));
        if e <= s
            continue
        end
        segment = betaData{ch}(s:e); 
        powerVals(end+1) = mean(segment .^ 2);
    end
    avgPower_1(ch) = mean(powerVals);

    powerVals = [];
    for i = 1:size(intervals_2_sec, 1)
        s = round(intervals_2_sec(i, 1) * SamplingFreq);
        e = round(intervals_2_sec(i, 2) * SamplingFreq);
        s = max(s, 1);
        e = min(e, length(betaData{ch}));
        if e <= s
            continue
        end
        segment = betaData{ch}(s:e);  
        powerVals(end+1) = mean(segment .^ 2);
    end
    avgPower_2(ch) = mean(powerVals);

    powerVals = [];
    for i = 1:size(intervals_3_sec, 1)
        s = round(intervals_3_sec(i, 1) * SamplingFreq);
        e = round(intervals_3_sec(i, 2) * SamplingFreq);
        s = max(s, 1);
        e = min(e, length(betaData{ch}));
        if e <= s
            continue
        end
        segment = betaData{ch}(s:e);  
        powerVals(end+1) = mean(segment .^ 2);
    end
    avgPower_3(ch) = mean(powerVals);

    
end






% Stack all three power vectors into a matrix (63 × 3)
allPower = [avgPower_1, avgPower_2, avgPower_3];

% Create channel grid layout (21 rows × 3 columns, top-left = ch 63, bottom-right = ch 1)
channel_order = reshape(63:-1:1, [21, 3]);




% Plot each condition
figure;
for idx = 1:3
    avgPower_per_channel = allPower(:, idx);  % get 63×1 vector for current condition
    gridPower = zeros(21, 3);  % initialize grid
    

    for row = 1:21
        for col = 1:3
        ch = channel_order(row, col);
        gridPower(row, col) = log10(avgPower_per_channel(ch) + eps);  
        end
    end


    % Plot
    subplot(1, 3, idx);
    imagesc(gridPower);
    axis equal tight;
    colormap('hot');
    colorbar;
    title(['Avg Beta Power - Condition ' num2str(idx)]);
    set(gca, 'XTick', [], 'YTick', []);

    % Optional: Add channel numbers
    for row = 1:21
        for col = 1:3
            ch = channel_order(row, col);
            text(col, row, num2str(ch), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'w', 'FontSize', 8);
        end
    end
end














% Load the NEV data
nev = load('PD22N009_nev.mat');

% Get the timestamps and event codes
timestamps_sec = nev.NEV.Data.SerialDigitalIO.TimeStampSec;
events = nev.NEV.Data.SerialDigitalIO.UnparsedData;  

intervals_1_sec = [];
intervals_2_sec = [];
intervals_3_sec = [];

for i = 1:length(events) - 1
    startTime = timestamps_sec(i);
    endTime   = timestamps_sec(i + 1);

    switch events(i)
        case 1
            intervals_1_sec(end+1, :) = [startTime, endTime];
        case 2
            intervals_2_sec(end+1, :) = [startTime, endTime];
        case 3
            intervals_3_sec(end+1, :) = [startTime, endTime];
    end
end

figure;
plot(tRel, Myeegfilt(dataAllChannels(33,:), SamplingFreq, 13, 30, 0, 1024), 'k');
hold on;

h1 = plotIntervalsSeconds(intervals_1_sec, 'r', 0.3);
h2 = plotIntervalsSeconds(intervals_2_sec, 'y', 0.3);
h3 = plotIntervalsSeconds(intervals_3_sec, 'g', 0.3);

xlabel('Time (s)');
ylabel('Amplitude');
title(['Filtered Signal with Task Intervals - Channel ' num2str(33)]);
legend([plot(NaN, 'k'), h1, h2, h3], {'Signal', 'Event 1', 'Event 2', 'Event 3'});


numChannels = 63;
dataLen = length(phaseData{1});
baselineEnd = round(dataLen / 3);  % first third

% Compute PLV over baseline window
PLV_vector_baseline = zeros(numChannels, numChannels);

for ch1 = 1:numChannels
    for ch2 = ch1+1:numChannels
        phase_diff = phaseData{ch1}(1:baselineEnd) - phaseData{ch2}(1:baselineEnd);
        PLV_vector_baseline(ch1, ch2) = mean(exp(1j * phase_diff));
        PLV_vector_baseline(ch2, ch1) = conj(PLV_vector_baseline(ch1,ch2));
    end
end

% Take absolute value to get PLV magnitude
PLV_value_baseline = abs(PLV_vector_baseline);
PLV_AvgValue_baseline = mean(PLV_value_baseline,2);


















intervals_1_samples = round(intervals_1_sec * SamplingFreq);
intervals_2_samples = round(intervals_2_sec * SamplingFreq);
intervals_3_samples = round(intervals_3_sec * SamplingFreq);

PLV_vector_1 = computePLVFromSampleIntervals(phaseData, intervals_1_samples);
PLV_vector_2 = computePLVFromSampleIntervals(phaseData, intervals_2_samples);
PLV_vector_3 = computePLVFromSampleIntervals(phaseData, intervals_3_samples);

PLV_value_1 = abs(PLV_vector_1);
PLV_value_2 = abs(PLV_vector_2);
PLV_value_3 = abs(PLV_vector_3);

avgPLV_1 = mean(PLV_value_1, 2);  % Mean PLV per channel
avgPLV_2 = mean(PLV_value_2, 2);  % Mean PLV per channel
avgPLV_3 = mean(PLV_value_3, 2);  % Mean PLV per channel

avgPLV_1_scaled = (avgPLV_1 - min(avgPLV_1)) / (max(avgPLV_1) - min(avgPLV_1));
avgPLV_2_scaled = (avgPLV_2 - min(avgPLV_2)) / (max(avgPLV_2) - min(avgPLV_2));
avgPLV_3_scaled = (avgPLV_3 - min(avgPLV_3)) / (max(avgPLV_3) - min(avgPLV_3));

avgPLV_baseline_scaled = (PLV_AvgValue_baseline - min(PLV_AvgValue_baseline))/ (max(PLV_AvgValue_baseline) - min(PLV_AvgValue_baseline));

PLV_Conditions_Avg = [avgPLV_1_scaled,avgPLV_2_scaled,avgPLV_3_scaled, avgPLV_baseline_scaled];

figure;
subplot(1,4,1);
imagesc(PLV_value_1); title('PLV - Event 1'); colorbar; caxis([0 1]);

subplot(1,4,2);
imagesc(PLV_value_2); title('PLV - Event 2'); colorbar; caxis([0 1]);

subplot(1,4,3);
imagesc(PLV_value_3); title('PLV - Event 3'); colorbar; caxis([0 1]);

subplot(1,4,4);
imagesc(PLV_value_baseline); title('PLV - Baseline'); colorbar; caxis([0 1]);


figure;
for idx = 1:4
    
    avgPLV_per_channel = PLV_Conditions_Avg(:,idx);
    channel_order = reshape(63:-1:1, [21, 3]);  % 21 rows, 3 columns
    gridPLV = zeros(21, 3);  % Initialize
    for row = 1:21
        for col = 1:3
            ch = channel_order(row, col);
            gridPLV(row, col) = avgPLV_per_channel(ch);
        end
    end
    
    subplot(1,4,idx);
    imagesc(gridPLV);
    axis equal tight;           % Make cells square and remove margins
    colormap('parula');
    colorbar;
    caxis([0 1]);               % PLV range
    title(['Average PLV per Channel - Condition ' num2str(idx)]);
    
    % Remove non-informative axis ticks
    set(gca, 'XTick', [], 'YTick', []);
    
    for row = 1:21
        for col = 1:3
            ch = channel_order(row, col);
            text(col, row, num2str(ch), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'Color', 'w', 'FontSize', 8);
        end
    end

end




% % PARAMETERS
% targetCh = 25;
% SamplingFreq = 1000;
% winLength = SamplingFreq * 1;  % 1 second window
% stepSize = 50;  % Optional downsampling for speed (every 50 ms)
% numChannels = length(phaseData);
% dataLen = length(phaseData{1});
% 
% % Initialize output
% plvTrace = [];
% timeVec = [];
% 
% % Start at 1s (to allow 1s of history)
% for t = winLength:stepSize:dataLen
%     winStart = t - winLength + 1;
%     winEnd = t;
% 
%     plvs = zeros(1, numChannels - 1);
%     c = 1;
%     for ch2 = 1:numChannels
%         if ch2 == targetCh
%             continue
%         end
%         phase_diff = phaseData{targetCh}(winStart:winEnd) - phaseData{ch2}(winStart:winEnd);
%         plvs(c) = abs(mean(exp(1j * phase_diff)));
%         c = c + 1;
%     end
% 
%     plvTrace(end+1) = mean(plvs);
%     timeVec(end+1) = tRel(winEnd);  % Use current time point
% end
% figure; hold on;
% 
% % Plot condition 1 intervals (red)
% for i = 1:size(intervals_1_sec, 1)
%     fill([intervals_1_sec(i,1), intervals_1_sec(i,2), intervals_1_sec(i,2), intervals_1_sec(i,1)], ...
%          [0, 0, 1, 1], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% end
% 
% % Plot condition 2 intervals (green)
% for i = 1:size(intervals_2_sec, 1)
%     fill([intervals_2_sec(i,1), intervals_2_sec(i,2), intervals_2_sec(i,2), intervals_2_sec(i,1)], ...
%          [0, 0, 1, 1], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% end
% 
% % Plot condition 3 intervals (blue)
% for i = 1:size(intervals_3_sec, 1)
%     fill([intervals_3_sec(i,1), intervals_3_sec(i,2), intervals_3_sec(i,2), intervals_3_sec(i,1)], ...
%          [0, 0, 1, 1], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% end
% 
% % Plot PLV trace
% plot(timeVec, plvTrace, 'k-', 'LineWidth', 1.5);
% 
% xlabel('Time (s)');
% ylabel('Avg PLV (last 1 sec)');
% title(['Causal Sliding Avg PLV - Channel ' num2str(targetCh)]);
% legend('Cond 1', 'Cond 2', 'Cond 3', 'PLV Trace');
% ylim([0 1]);
% xlim([tRel(1), tRel(end)]);
% grid on;

 
















function h = plotIntervalsSeconds(intervals_sec, color, alpha)
    yLimits = ylim;
    h = gobjects(size(intervals_sec,1), 1);  % Preallocate handles
    for i = 1:size(intervals_sec, 1)
        x1 = intervals_sec(i,1);
        x2 = intervals_sec(i,2);
        h(i) = fill([x1 x2 x2 x1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
            color, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    end
    h = h(1);  % Return one handle for legend
end




function PLV_vector_matrix = computePLVFromSampleIntervals(phaseData, intervals)
    numChannels = length(phaseData);
    PLV_vector_matrix = zeros(numChannels, numChannels);

    for ch1 = 1:numChannels
        for ch2 = ch1+1:numChannels
            phase_diff_all = [];
            for k = 1:size(intervals, 1)
                idx_start = intervals(k, 1);
                idx_end = intervals(k, 2);
                % Clip to signal bounds
                idx_end = min(idx_end, length(phaseData{ch1}));
                phase_diff = phaseData{ch1}(idx_start:idx_end) - phaseData{ch2}(idx_start:idx_end);
                phase_diff_all = [phase_diff_all, phase_diff];
            end
            PLV_vector_matrix(ch1, ch2) = mean(exp(1j * phase_diff_all));
            PLV_vector_matrix(ch2, ch1) = conj(PLV_vector_matrix(ch1, ch2));
        end
    end
end





