%PD23N003 - 'pds002.ns2'
%PD22N009 - 'datafile001.ns2'
%PD22N008 - 'PD22N008.ns2'
%PD23N009 - 'PDS001.ns2

patient_data = 'PD22N008.ns2';

%PDN008 - [14, 35, 36, 39, 56, 57, 59, 60, 61, 62]
%PDN009 - [14, 15, 16, 17, 18, 33, 34, 35, 36, 37, 38, 53, 54, 55, 56, 57,
%58, 59]
%PD23N003 - [56, 57, 58]
%PD23N009 - [35, 54, 55, 56]

highlighted_channels =  [14, 35, 36, 56, 57];  % Example list of channels to highlight in green
highlight_region_name = 'Upper Limb Region';

freq_range = 'beta'; %'beta' or 'theta'

[dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS(patient_data,1);

disp(size(dataAllChannels,1))
channelIndices = 1:size(dataAllChannels,1);
numChannels = size(dataAllChannels,1) - 1;
disp(numChannels);

phaseData = cell(1, numChannels);

disp(channelNames');


%%%%%%%% Store phases

for idx = 1:numChannels
    ch_raw = dataAllChannels(idx, :); % Raw signal
    ch_raw_cut = ch_raw(1:round(length(ch_raw)/3)); %Take first half of data (baseline data)

    if strcmp(freq_range, 'theta')
        ch_filtered = Myeegfilt(ch_raw_cut, SamplingFreq, 4, 9, 0, 1024); % Filter to theta band (4-9 Hz)
    elseif strcmp(freq_range, 'beta')
        ch_filtered = Myeegfilt(ch_raw_cut, SamplingFreq, 13, 30, 0, 1024); % Filter to beta band (13-30 Hz)
    else
        error('Unknown frequency range: %s', freq_range);
    end

    ch_one = ch_filtered(1,:);
    [ch_phase, ~] = instPhaseFreq(ch_filtered, SamplingFreq); % Compute phase
    phaseData{idx} = ch_phase; % Store phase data
    disp(idx)
end
disp('Finished storing phases')

%%%%%%%%% Store PLV

% Compute PLV vector matrix (stores mean complex phase difference vectors)
PLV_vector_matrix = zeros(numChannels, numChannels); % Initialize matrix

for ch1 = 1:numChannels
    for ch2 = ch1+1:numChannels  % Only upper triangle since PLV is symmetric
        phase_diff = phaseData{ch1} - phaseData{ch2};  % Compute phase differences
        PLV_vector_matrix(ch1, ch2) = mean(exp(1j * phase_diff)); % Store mean complex vector
        PLV_vector_matrix(ch2, ch1) = conj(PLV_vector_matrix(ch1, ch2)); % Make it symmetric
    end
end

PLV_value_matrix = abs(PLV_vector_matrix);
disp('Done storing PLVs')


best_neighbor_by_plv = zeros(numChannels, 1);
max_plv_values = zeros(numChannels, 1);

for ch = 1:numChannels
    plvs = PLV_value_matrix(ch, :);
    plvs(ch) = -1;  % Exclude self
    [max_plv, best_ch] = max(plvs);
    best_neighbor_by_plv(ch) = best_ch;
    max_plv_values(ch) = max_plv;
end

circ_std = @(x) sqrt(-2 * log(abs(mean(exp(1j * x))))); % circular std
best_neighbor_by_std = zeros(numChannels, 1);
min_std_values = zeros(numChannels, 1);

for ref_ch = 1:numChannels
    std_values = Inf(1, numChannels); % Init with Inf, will be updated

    for ch = 1:numChannels
        if ch == ref_ch
            continue;
        end
        phase_diff = angle(exp(1j * (phaseData{ch} - phaseData{ref_ch})));  % Circular diff
        std_values(ch) = circ_std(phase_diff);
    end

    [min_std, best_std_ch] = min(std_values);
    best_neighbor_by_std(ref_ch) = best_std_ch;
    min_std_values(ref_ch) = min_std;
end


grid_width = 3;
grid_height = 21;
num_channels = grid_width * grid_height;
electrode_positions = cell(num_channels, 1);
channel_num = 1;

for x = 0:grid_width-1
    for y = 0:grid_height-1  % Bottom-up order
        electrode_positions{channel_num} = [x, y]; % Store as a 1x2 array
        channel_num = channel_num + 1;
    end
end



figure; 

% --- Common setup ---
cmap = parula(256);

% Compute centers of boxes
X = zeros(num_channels, 1);
Y = zeros(num_channels, 1);
for i = 1:num_channels
    X(i) = electrode_positions{i}(1) + 0.5;
    Y(i) = electrode_positions{i}(2) + 0.5;
end


% --- Subplot 1: Best PLV Neighbor ---
subplot(1,2,1); hold on;
title(sprintf('Best Neighbor by PLV — %s (%s)', highlight_region_name, patient_data));

% Draw rectangles
for i = 1:num_channels
    pos = electrode_positions{i};
    if ismember(i, highlighted_channels)
        face_col = [0.6 1 0.6];  % light green
    else
        face_col = [0.9 0.9 0.9];  % neutral gray
    end
    
    rectangle('Position', [pos(1), pos(2), 1, 1], ...
              'FaceColor', face_col, 'EdgeColor', 'k');
end

% Draw PLV arrows
for ch = 1:num_channels
    from_pos = electrode_positions{ch};
    target_ch = best_neighbor_by_plv(ch);

    if target_ch < 1 || target_ch > num_channels, continue; end

    to_pos = electrode_positions{target_ch};
    dx = to_pos(1) - from_pos(1);
    dy = to_pos(2) - from_pos(2);
    
    if ch == 3
        continue; 
    end

    quiver(from_pos(1) + 0.5, from_pos(2) + 0.5, ...
           dx, dy, 0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 2);
end

% Label channels
for i = 1:num_channels
    text(X(i), Y(i), num2str(i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k');
end

axis equal; xlim([0, grid_width]); ylim([0, grid_height]);
set(gca, 'XTick', [], 'YTick', []); set(gca, 'YDir', 'normal');


% --- Subplot 2: Best STD Neighbor ---
subplot(1,2,2); hold on;
title(sprintf('Best Neighbor by Circular STD — %s (%s)', highlight_region_name, patient_data));


% Draw rectangles again
for i = 1:num_channels
    pos = electrode_positions{i};
    if ismember(i, highlighted_channels)
        face_col = [0.6 1 0.6];  % light green
    else
        face_col = [0.9 0.9 0.9];  % neutral gray
    end

    
    rectangle('Position', [pos(1), pos(2), 1, 1], ...
              'FaceColor', face_col, 'EdgeColor', 'k');

end

% Draw STD arrows
for ch = 1:num_channels
    from_pos = electrode_positions{ch};
    target_ch = best_neighbor_by_std(ch);

    if target_ch < 1 || target_ch > num_channels, continue; end

    to_pos = electrode_positions{target_ch};
    dx = to_pos(1) - from_pos(1);
    dy = to_pos(2) - from_pos(2);

    if ch == 3
        continue; 
    end

    quiver(from_pos(1) + 0.5, from_pos(2) + 0.5, ...
           dx, dy, 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2);
end

% Label channels
for i = 1:num_channels
    text(X(i), Y(i), num2str(i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'FontWeight', 'bold', 'Color', 'k');
end

axis equal; xlim([0, grid_width]); ylim([0, grid_height]);
set(gca, 'XTick', [], 'YTick', []); set(gca, 'YDir', 'normal');



%ROSE PLOTS
% circ_stats = @(x) deal( ...
%     mean(x), ...  % linear mean (not angle of complex vector)
%     std(x));      % linear std
% 
% figure; 
% num_highlights = numel(highlighted_channels);
% 
% for idx = 1:num_highlights
%     ref_ch = highlighted_channels(idx);
%     if ref_ch == 3
%         continue;  % Skip buggy channel
%     end
% 
%     neighbor_ch = best_neighbor_by_plv(ref_ch);
%     if neighbor_ch < 1 || neighbor_ch > num_channels
%         continue;  % Skip invalid
%     end
% 
%     % Compute linear phase difference (wrapped between [-pi, pi])
%     phase_diff = angle(exp(1j * (phaseData{neighbor_ch} - phaseData{ref_ch})));
% 
%     % Compute linear stats
%     [mu, std_dev] = circ_stats(phase_diff);
%     n = numel(phase_diff);
%     ci95 = 1.96 * std_dev / sqrt(n);  % Linear 95% CI on the mean
% 
%     subplot(1, num_highlights, idx);
% 
%     % --- Create histogram and capture data ---
%     h = polarhistogram(phase_diff, 75, 'Normalization', 'probability');
%     edges = h.BinEdges;
%     counts = h.Values;
%     max_r = max(counts);  % Use max for scaling arrows/lines
% 
%     % Set radial limit (optional: tweak for padding)
%     rlim([0 max_r * 1.1]);
% 
%     hold on;
% 
%     % --- Plot mean direction line ---
%     polarplot([mu mu], [0 max_r], 'r', 'LineWidth', 2);
% 
%     % --- Plot 95% CI lines ---
%     polarplot([mu - ci95, mu - ci95], [0 max_r * 0.8], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
%     polarplot([mu + ci95, mu + ci95], [0 max_r * 0.8], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
% 
%     % --- Title with bounds ---
%     ci_lower = mu - ci95;
%     ci_upper = mu + ci95;
%     title({
%         sprintf('Ch %d vs %d', ref_ch, neighbor_ch), ...
%         sprintf('\\mu = %.2f rad, \\sigma = %.2f', mu, std_dev), ...
%         sprintf('95%% CI = [%.2f, %.2f] rad', ci_lower, ci_upper)
% });
% 
% end
circ_stats = @(x) deal(mean(x), std(x));  % Linear stats (since we’re not using circular CI)
figure;
num_highlights = numel(highlighted_channels);

for idx = 1:num_highlights
    ref_ch = highlighted_channels(idx);
    if ref_ch == 3
        continue;  % Skip buggy
    end

    % --- First subplot (PLV best neighbor) ---
    neighbor_ch_plv = best_neighbor_by_plv(ref_ch);
    if neighbor_ch_plv < 1 || neighbor_ch_plv > num_channels
        continue;
    end

    phase_diff_plv = angle(exp(1j * (phaseData{neighbor_ch_plv} - phaseData{ref_ch})));
    [mu_plv, std_plv] = circ_stats(phase_diff_plv);
    n = numel(phase_diff_plv);
    ci95_plv = 1.96 * std_plv / sqrt(n);
    ci_low_plv = mu_plv - ci95_plv;
    ci_high_plv = mu_plv + ci95_plv;

    subplot(2, num_highlights, idx);
    h = polarhistogram(phase_diff_plv, 75, 'Normalization', 'probability');
    max_r = max(h.Values); rlim([0 max_r * 1.1]);
    hold on;
    polarplot([mu_plv mu_plv], [0 max_r], 'r', 'LineWidth', 2);
    polarplot([ci_low_plv ci_low_plv], [0 max_r * 0.8], '--', 'Color', [0.5 0.5 0.5]);
    polarplot([ci_high_plv ci_high_plv], [0 max_r * 0.8], '--', 'Color', [0.5 0.5 0.5]);
    title({
        sprintf('Ch %d vs PLV Best (Ch %d)', ref_ch, neighbor_ch_plv), ...
        sprintf('\\mu = %.2f, \\sigma = %.2f', mu_plv, std_plv), ...
        sprintf('95%% CI = [%.2f, %.2f]', ci_low_plv, ci_high_plv)
    });


    % --- Second subplot (naive "above" neighbor) ---
    pos = electrode_positions{ref_ch};
    col = pos(1); row = pos(2);
    above_row = row + 1;

    % Find channel at [col, row+1]
    above_ch = [];
    for ch_test = 1:num_channels
        test_pos = electrode_positions{ch_test};
        if test_pos(1) == col && test_pos(2) == above_row
            above_ch = ch_test;
            break;
        end
    end

    if ~isempty(above_ch)
        phase_diff_above = angle(exp(1j * (phaseData{above_ch} - phaseData{ref_ch})));
        [mu_above, std_above] = circ_stats(phase_diff_above);
        ci95_above = 1.96 * std_above / sqrt(n);
        ci_low_above = mu_above - ci95_above;
        ci_high_above = mu_above + ci95_above;

        subplot(2, num_highlights, idx + num_highlights);
        h2 = polarhistogram(phase_diff_above, 75, 'Normalization', 'probability');
        max_r2 = max(h2.Values); rlim([0 max_r2 * 1.1]);
        hold on;
        polarplot([mu_above mu_above], [0 max_r2], 'r', 'LineWidth', 2);
        polarplot([ci_low_above ci_low_above], [0 max_r2 * 0.8], '--', 'Color', [0.5 0.5 0.5]);
        polarplot([ci_high_above ci_high_above], [0 max_r2 * 0.8], '--', 'Color', [0.5 0.5 0.5]);
        title({
            sprintf('Ch %d vs Above (Ch %d)', ref_ch, above_ch), ...
            sprintf('\\mu = %.2f, \\sigma = %.2f', mu_above, std_above), ...
            sprintf('95%% CI = [%.2f, %.2f]', ci_low_above, ci_high_above)
        });

    else
        subplot(2, num_highlights, idx + num_highlights);
        title(sprintf('No above neighbor for Ch %d', ref_ch));
    end
end
























% % Step 1: Get the filtered signals used for phase extraction
% channelIndices = 1:63;
% group_channels = length(channelIndices);
% group_signals = zeros(length(group_channels), length(phaseData{1}));
% 
% for i = 1:length(group_channels)
%     ch = group_channels(i);
%     group_signals(i, :) = dataAllChannels(ch, 1:round(length(dataAllChannels(ch, :))/3)); % same cut
% end
% 
% % Step 2: Filter each channel to beta band
% for i = 1:size(group_signals, 1)
%     group_signals(i, :) = Myeegfilt(group_signals(i, :), SamplingFreq, 13, 30, 0, 1024);
% end
% 
% 
% 
% 
% 
% 
% 
% 
% ref_channel = 24;%24, best is 25, -0.0027 offset, [3,25,45,23] group
% best_neighbor_ch = 25;
% best_neighbor_offset = 0;
% %best_neighbor_offset = 0.0192;  %36 ref, 15 best, 0.0192 offset, [15,37,35,57] group
% 
% % You already have:
% ref_phase = angle(exp(1j * phaseData{ref_channel}));
% best_neighbor_phase = angle(exp(1j * phaseData{best_neighbor_ch}));
% other_ch = 18;
% other_phase = angle(exp(1j * phaseData{other_ch}));
% other_phase_diff = mean(angle(exp(1j * (phaseData{other_ch} - phaseData{ref_channel})))); % circular mean
% 
% % Identify peak times in reference channel
% peak_indices = find(abs(ref_phase) < pi/50);  % Adjust threshold if needed
% 
% group_channels = [3,25,45,23]; %[3,25,45,23] 24 ref 25 best -0.0027 offset
% 
% % Stack the complex phase vectors
% group_complex_phases = cellfun(@(ph) exp(1j * ph), phaseData(group_channels), 'UniformOutput', false);
% group_complex_phases = cat(1, group_complex_phases{:});  % Shape: [nChannels x nTime]
% 
% % Average the complex vectors across channels
% group_avg_vector = mean(group_complex_phases, 1);  % Shape: [1 x nTime] 0.008, .41
% 
% % Convert back to phase
% group_phase = angle(group_avg_vector);  % This is your group mean phase
% 
% 
% % Step 3: Average the raw signals across channels
% signal_avg = mean(group_signals, 1);
% 
% % Step 4: Get the analytic phase of the combined signal
% [combined_phase, ~] = instPhaseFreq(signal_avg, SamplingFreq);
% combined_phase = angle(exp(1j * combined_phase)); % wrapped to [-pi, pi]
% 
% % Step 5: Sample at ref peaks
% signal_avg_group_phase = combined_phase(peak_indices);
% 
% % Step 6: Compute circular stats
% [mu_signalavg, std_signalavg] = circ_stats(signal_avg_group_phase);
% 
% 
% % -----------------------------
% % 1. Best neighbor - with phase offset correction
% corrected_best = angle(exp(1j * (best_neighbor_phase(peak_indices) - best_neighbor_offset)));
% 
% % 2. Best neighbor - without phase correction
% uncorrected_best = angle(exp(1j * best_neighbor_phase(peak_indices)));
% 
% % 3. Other channel (e.g., 18) - with phase offset correction
% corrected_other = angle(exp(1j * (other_phase(peak_indices))));
% 
% % Align group to ref by removing average phase difference
% corrected_group_phase = angle(exp(1j * (group_phase(peak_indices))));
% 
% 
% 
% % Helper function for circular stats
% circ_stats = @(x) deal( ...
%     angle(mean(exp(1j * x))), ...                          % circular mean
%     sqrt(-2 * log(abs(mean(exp(1j * x))))));               % circular std
% 
% % Compute stats for each distribution
% [mu_best, std_best] = circ_stats(corrected_best);
% [mu_best_uncorrected, std_best_uncorrected] = circ_stats(uncorrected_best);
% [mu_other, std_other] = circ_stats(corrected_other);
% [mu_group, std_group] = circ_stats(corrected_group_phase);
% 
% 
% 
% 
% 
% % -----------------------------
% % Plot all 3 rose plots
% subplot(1,5,1);
% polarhistogram(corrected_best, 75, 'Normalization', 'probability');
% title({
%     sprintf('Corrected - Best (Ch %d)', best_neighbor_ch), ...
%     sprintf('\\mu = %.3f rad, \\sigma = %.2f', mu_best, std_best)
% });
% 
% subplot(1,5,2);
% polarhistogram(uncorrected_best, 75, 'Normalization', 'probability');
% title({
%     sprintf('Uncorrected - Best (Ch %d)', best_neighbor_ch), ...
%     sprintf('\\mu = %.3f rad, \\sigma = %.2f', mu_best_uncorrected, std_best_uncorrected)
% });
% 
% subplot(1,5,3);
% polarhistogram(corrected_other, 75, 'Normalization', 'probability');
% title({
%     sprintf('Other (Ch %d)', other_ch), ...
%     sprintf('\\mu = %.3f rad, \\sigma = %.2f', mu_other, std_other)
% });
% 
% subplot(1,5,4);
% polarhistogram(corrected_group_phase, 75, 'Normalization', 'probability');
% title({
%     sprintf('Corrected - Group Avg (Ch %s)', join(string(group_channels), ',')), ...
%     sprintf('\\mu = %.3f rad, \\sigma = %.2f', mu_group, std_group)
% });
% 
% subplot(1,5,5);
% polarhistogram(signal_avg_group_phase, 75, 'Normalization', 'probability');
% title({
%     sprintf('Signal-Averaged Group (Ch %s)', join(string(group_channels), ',')), ...
%     sprintf('\\mu = %.2f rad, \\sigma = %.2f', mu_signalavg, std_signalavg)
% });
% 
