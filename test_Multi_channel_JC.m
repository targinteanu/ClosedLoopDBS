%addpath(genpath('C:\Users\cours\OneDrive\Desktop\Anderson Lab\NPMK-5.5.5.0'));

[dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS('PD22N008.ns2');

%Channels
recordingChannelIndex = 37; % Primary recording channel
neighborChannels = [16, 36, 38, 58]; % Neighbor recording channels - 1, 21, 43, 63
disp(neighborChannels(1))
disp(neighborChannels(2))
disp(neighborChannels(3))
disp(neighborChannels(4))

%Raw signals
ch58_raw = dataAllChannels(recordingChannelIndex, :);
ch_neighbor1_raw = dataAllChannels(neighborChannels(1), :);
ch_neighbor2_raw = dataAllChannels(neighborChannels(2), :);
ch_neighbor3_raw = dataAllChannels(neighborChannels(3), :);
ch_neighbor4_raw = dataAllChannels(neighborChannels(4), :);
ch_neighbors_all_raw = [ch_neighbor1_raw; 
                         ch_neighbor2_raw; 
                         ch_neighbor3_raw; 
                         ch_neighbor4_raw];
ch_neighbors_avg_raw = mean(ch_neighbors_all_raw, 1);

%Filtered signals
ch58_filtered = Myeegfilt(ch58_raw, SamplingFreq, 13, 30, 0, 1024);
ch_neighbor1_filtered = Myeegfilt(ch_neighbor1_raw, SamplingFreq, 13, 30, 0, 1024);
ch_neighbor2_filtered = Myeegfilt(ch_neighbor2_raw, SamplingFreq, 13, 30, 0, 1024);
ch_neighbor3_filtered = Myeegfilt(ch_neighbor3_raw, SamplingFreq, 13, 30, 0, 1024);
ch_neighbor4_filtered = Myeegfilt(ch_neighbor4_raw, SamplingFreq, 13, 30, 0, 1024);
ch_neighbors_all_filtered = [ch_neighbor1_filtered; 
                         ch_neighbor2_filtered; 
                         ch_neighbor3_filtered; 
                         ch_neighbor4_filtered];    
ch_neighbors_avg_filtered = mean(ch_neighbors_all_filtered, 1);
ch_neighbors_avg_raw_filtered = Myeegfilt(ch_neighbors_avg_raw, SamplingFreq, 13, 30, 0, 1024);

%Instantaneous Phase and Frequency
[ch58_Phase, ch58_Freq] = instPhaseFreq(ch58_filtered, SamplingFreq);
[neighborsPhase, neighborsFreq] = instPhaseFreq(ch_neighbors_avg_filtered, SamplingFreq);

[neighbor1Phase, neighbor1Freq] = instPhaseFreq(ch_neighbor1_filtered, SamplingFreq);
[neighbor2Phase, neighbor2Freq] = instPhaseFreq(ch_neighbor2_filtered, SamplingFreq);
[neighbor3Phase, neighbor3Freq] = instPhaseFreq(ch_neighbor3_filtered, SamplingFreq);
[neighbor4Phase, neighbor4Freq] = instPhaseFreq(ch_neighbor4_filtered, SamplingFreq);
neighbor_individual_phases = [neighbor1Phase; 
                         neighbor2Phase; 
                         neighbor3Phase; 
                         neighbor4Phase];
neighbor_individual_freqs = [neighbor1Freq; 
                         neighbor2Freq; 
                         neighbor3Freq; 
                         neighbor4Freq];
neighborsPhase_2 = mean(neighbor_individual_phases, 1);
neighborsFreq_2 = mean(neighbor_individual_freqs, 1);

startTime = t(20000);                  
numSamples = round(SamplingFreq);     
[~, startIdx] = min(abs(t - startTime));
endIdx = min(startIdx + numSamples - 1, length(t));

endingIndex = 250000;
t = t(1:endingIndex);
ch58_raw = ch58_raw(1:endingIndex);
ch58_Phase = ch58_Phase(1:endingIndex);
ch58_Freq = ch58_Freq(1:endingIndex);
ch58_filtered = ch58_filtered(1:endingIndex);
neighborsPhase = neighborsPhase(1:endingIndex);
neighborsPhase_2 = neighborsPhase_2(1:endingIndex);
neighborsFreq = neighborsFreq(1:endingIndex);
neighborsFreq_2 = neighborsFreq_2(1:endingIndex);

%Figure 1
figure;
tiledlayout(6,1);
% ---- Plot 1: Ch58 raw signal vs. Ch58 filtered signal ----
nexttile;
plot(t, ch58_raw, 'k', 'LineWidth', 1); 
hold on;
plot(t, ch58_filtered, 'b', 'LineWidth', 1.2);
title('Ch58: Raw vs. Filtered Signal');
xlabel('Time (s)');
ylabel('uV');
legend('Ch58 - Raw Signal', 'Ch 58 - Filtered Signal');
grid on;

% ---- Plot 2: Ch58 filtered signal vs. its phase ----
nexttile;
hold on;  
plot(t(startIdx:endIdx), ch58_filtered(startIdx:endIdx), 'b', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), ch_neighbor1_filtered(startIdx:endIdx), 'r', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), ch_neighbor2_filtered(startIdx:endIdx), 'g', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), ch_neighbor3_filtered(startIdx:endIdx), 'k', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), ch_neighbor4_filtered(startIdx:endIdx), 'm', 'LineWidth', 1.2); 
ylabel('uV');
xlabel('Time (s)');
title('Ch58 Filtered vs. Neighbor Channels Filtered');
legend('Ch58 Filtered', 'Neighbor 1', 'Neighbor 2', 'Neighbor 3', 'Neighbor 4');
grid on;
hold off;

% ---- Plot 3: Ch58 filtered signal vs. Avg Neighbor Signal Filtered
nexttile;
hold on;
plot(t(startIdx:endIdx), ch58_filtered(startIdx:endIdx), 'b', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), ch_neighbors_avg_filtered(startIdx:endIdx), 'c', 'LineWidth', 1.2); 
ylabel('uV');
xlabel('Time (s)');
title('Ch58 Filtered vs. Avg Neighbor Chs Filtered - Method 1');
legend('Ch58 Filtered', 'Avg Neighbor Chs Filtered'); 
grid on;
hold off;

% ---- Plot 4: Ch58 phase vs. Neighbor phase
nexttile;
yyaxis left;
plot(t(startIdx:endIdx), ch58_Phase(startIdx:endIdx), 'b', 'LineWidth', 1.2);
ylabel('Phase (rad)');
yyaxis right;
plot(t(startIdx:endIdx), neighborsPhase(startIdx:endIdx), 'r', 'LineWidth', 1.2); 
ylabel('Phase (rad)');
title('Ch58 Phase vs. Avg Neighbors Phase - Method 1');
xlabel('Time (s)');
legend('Ch58', 'Avg Neighbor');
grid on;

% ---- Plot 5: Ch58 Freq vs. Neighbor Freq
nexttile;
plot(t(startIdx:endIdx), ch58_Freq(startIdx:endIdx), 'b', 'LineWidth', 1.2); 
ylabel('Freq (Hz)');
hold on;
plot(t(startIdx:endIdx), neighborsFreq(startIdx:endIdx), 'r', 'LineWidth', 1.2);
ylabel('Freq (Hz)')
title('Ch58 Freq vs. Avg Neighbors Freq - Method 1');
xlabel('Time (s)');
legend('Ch58', 'Avg Neighbor');
grid on;

%Figure 2
figure;
tiledlayout(4,1);

nexttile;
hold on;  
plot(t(startIdx:endIdx), ch58_Phase(startIdx:endIdx), 'b', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), neighbor1Phase(startIdx:endIdx), 'r', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), neighbor2Phase(startIdx:endIdx), 'g', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), neighbor3Phase(startIdx:endIdx), 'k', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), neighbor4Phase(startIdx:endIdx), 'm', 'LineWidth', 1.2); 
ylabel('Phase (rad)');
xlabel('Time (s)');
title('Ch58 Phase vs. Neighbor Channels Phases');
legend('Ch58', 'Neighbor 1', 'Neighbor 2', 'Neighbor 3', 'Neighbor 4');
grid on;
hold off;

nexttile;
hold on;  
plot(t(startIdx:endIdx), ch58_Phase(startIdx:endIdx), 'b', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), neighbor1Freq(startIdx:endIdx), 'r', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), neighbor2Freq(startIdx:endIdx), 'g', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), neighbor3Freq(startIdx:endIdx), 'k', 'LineWidth', 1.2); 
plot(t(startIdx:endIdx), neighbor4Freq(startIdx:endIdx), 'm', 'LineWidth', 1.2); 
ylabel('Freq (Hz)');
xlabel('Time (s)');
title('Ch58 Phase vs. Neighbor Channels Freqs');
legend('Ch58', 'Neighbor 1', 'Neighbor 2', 'Neighbor 3', 'Neighbor 4');
grid on;
hold off;

nexttile;
yyaxis left;
plot(t(startIdx:endIdx), ch58_Phase(startIdx:endIdx), 'b', 'LineWidth', 1.2); 
ylabel('Phase (rad)');
yyaxis right;
plot(t(startIdx:endIdx), neighborsPhase_2(startIdx:endIdx), 'r', 'LineWidth', 1.2); 
ylabel('Phase (rad)');
title('Ch58 Phase vs. Avg Neighbors Phase - Method 2');
xlabel('Time (s)');
legend('Ch58', 'Avg Neighbor');
grid on;

nexttile;
yyaxis left;
plot(t(startIdx:endIdx), ch58_Freq(startIdx:endIdx), 'b', 'LineWidth', 1.2); 
ylabel('Freq (Hz)');
yyaxis right;
plot(t(startIdx:endIdx), neighborsFreq_2(startIdx:endIdx), 'r', 'LineWidth', 1.2);
ylabel('Freq (Hz)');
title('Ch58 Freq vs. Avg Neighbors Freq - Method 2');
xlabel('Time (s)');
legend('Ch58', 'Avg Neighbor');
grid on;




%Figure 3
phaseError1 = ch58_Phase - neighborsPhase;
freqError1 = ch58_Freq - neighborsFreq;

phaseError2 = ch58_Phase - neighborsPhase_2;
freqError2 = ch58_Freq - neighborsFreq_2;

MAE_phase1 = mean(abs(phaseError1));
MAE_freq1 = mean(abs(freqError1));
MAE_phase2 = mean(abs(phaseError2));
MAE_freq2 = mean(abs(freqError2));
RMSE_phase1 = sqrt(mean(phaseError1.^2));
RMSE_freq1 = sqrt(mean(freqError1.^2));
RMSE_phase2 = sqrt(mean(phaseError2.^2));
RMSE_freq2 = sqrt(mean(freqError2.^2));

figure;
tiledlayout(4,1);

% ---- Plot 1: Phase Error Distribution ----
nexttile;
histogram(phaseError1, 150, 'FaceColor', 'b', 'EdgeColor', 'k'); % 50 bins
xlabel('Phase Error (radians)');
ylabel('Count');
title(['Phase Error Distribution (Ch58 - Avg Neighbor) - Method 1' ...
    sprintf('\nMAE = %.3f, RMSE = %.3f', MAE_phase1, RMSE_phase1)]);
grid on;

% ---- Plot 2: Frequency Error Distribution ----
nexttile;
histogram(freqError1, 150, 'FaceColor', 'r', 'EdgeColor', 'k'); % 50 bins
xlabel('Frequency Error (Hz)');
ylabel('Count');
title(['Frequency Error Distribution (Ch58 - Avg Neighbor) - Method 1' ...
    sprintf('\nMAE = %.3f, RMSE = %.3f', MAE_freq1, RMSE_freq1)]);
grid on;

% ---- Plot 3: Phase Error Distribution ----
nexttile;
histogram(phaseError2, 150, 'FaceColor', 'b', 'EdgeColor', 'k'); % 50 bins
xlabel('Phase Error (radians)');
ylabel('Count');
title(['Phase Error Distribution (Ch58 - Avg Neighbor) - Method 2' ...
    sprintf('\nMAE = %.3f, RMSE = %.3f', MAE_phase2, RMSE_phase2)]);
grid on;

% ---- Plot 4: Frequency Error Distribution ----
nexttile;
histogram(freqError2, 150, 'FaceColor', 'r', 'EdgeColor', 'k'); % 50 bins
xlabel('Frequency Error (Hz)');
ylabel('Count');
title(['Frequency Error Distribution (Ch58 - Avg Neighbor) - Method 2' ...
    sprintf('\nMAE = %.3f, RMSE = %.3f', MAE_freq2, RMSE_freq2)]);
grid on;






