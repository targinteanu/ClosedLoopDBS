% freq_range = 'beta'; %'beta' or 'theta'
% 
% disp(label_str);
% 
% 
% [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
%     channelName, channelIndex, channelIndexStim, channelNames]...
%     = getRecordedData_NS('memory.ns2',1);
% 
% disp(size(dataAllChannels,1))
% channelIndices = 1:size(dataAllChannels,1);
% numChannels = size(dataAllChannels,1);
% phaseData = cell(1, numChannels);
% 
% disp(channelNames');
% 
% reference_ch = find(contains(channelNames, 'RHA3'));
% disp(reference_ch);
% 
% %%%%%%%% Store phases
% 
% for idx = 1:numChannels
%     ch_raw = dataAllChannels(idx, :); % Raw signal
%     ch_raw_cut = PD22N009_BL(idx, :)
%     ch_raw_cut = ch_raw(1:round(length(ch_raw)/3)); %Take first half of data (baseline data)
% 
% 
%     if strcmp(freq_range, 'theta')
%         ch_filtered = Myeegfilt(ch_raw_cut, SamplingFreq, 4, 9, 0, 1024); % Filter to theta band (4-9 Hz)
%     elseif strcmp(freq_range, 'beta')
%         ch_filtered = Myeegfilt(ch_raw_cut, SamplingFreq, 13, 30, 0, 1024); % Filter to beta band (13-30 Hz)
%     else
%         error('Unknown frequency range: %s', freq_range);
%     end
% 
%     ch_one = ch_filtered(1,:);
%     [ch_phase, ~] = instPhaseFreq(ch_filtered, SamplingFreq); % Compute phase
%     phaseData{idx} = ch_phase; % Store phase data
%     disp(idx)
% end
% disp('Finished storing phases')

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

%Plot PLV Value Matrix
figure;
imagesc(PLV_value_matrix); % Display the PLV matrix
colorbar; % Add a color scale
colormap('parula'); % Replace 'viridis' with 'parula' or any other built-in colormap
caxis([0 1]); % Set color scale between 0 and 1 for better interpretation
numChannels = size(PLV_value_matrix, 1);
xticks(1:numChannels);
yticks(1:numChannels);
xticklabels(1:numChannels); 
yticklabels(1:numChannels);
xtickangle(90);
xlabel('Channel');
ylabel('Channel');
title('Phase Locking Value (PLV) Matrix');



%%%%%%%%% Compute PLV Variance Matrix
PLV_variance_matrix = zeros(numChannels, numChannels);

for ch1 = 1:numChannels
    for ch2 = ch1+1:numChannels  % Only compute upper triangle (symmetric)
        phase_diff = phaseData{ch1} - phaseData{ch2};  % Compute phase differences
        PLV_values = abs(exp(1j * phase_diff));  % Compute PLV values (no time windowing)
        PLV_variance_matrix(ch1, ch2) = var(PLV_values);  % Compute variance
        PLV_variance_matrix(ch2, ch1) = PLV_variance_matrix(ch1, ch2);  % Make it symmetric
    end
end

disp('Finished PLV variance matrix');

% Plot PLV Variance Matrix
figure;
imagesc(PLV_variance_matrix); % Display the PLV variance matrix
colorbar; % Add a color scale
colormap('parula'); % Use 'parula' colormap
caxis([0 max(PLV_variance_matrix(:))]); % Set color scale based on max variance for better interpretation
numChannels = size(PLV_variance_matrix, 1);
xticks(1:numChannels);
yticks(1:numChannels);
xticklabels(1:numChannels);
yticklabels(1:numChannels);
xtickangle(90); % Rotate x-axis labels for readability
xlabel('Channel');
ylabel('Channel');
title('PLV Variance Matrix');









%%%%%%%%% Plot Spatial Synchrony Map

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

vector_field = cell(num_channels, 1);
for i = 1:num_channels
    vector_field{i} = [0, 0]; % Initialize each as a 1x2 vector
end

for i = 1:num_channels
    x_i = electrode_positions{i}(1);
    y_i = electrode_positions{i}(2);

    for j = 1:num_channels
        if i == j
            continue; % Skip self-comparison
        end
        
        x_j = electrode_positions{j}(1);
        y_j = electrode_positions{j}(2);
        
        % Compute direction vector
        direction_vector = [x_j - x_i, y_j - y_i];
        norm_factor = norm(direction_vector);
        if norm_factor == 0
            continue; % Avoid division by zero
        end
        direction_vector = direction_vector / norm_factor; % Normalize
        
        % Get phase difference and PLV
        plv_direction = sign(angle(PLV_vector_matrix(i,j)));
        plv_weight = PLV_value_matrix(i, j);

        % Compute weighted phase gradient vector
        weighted_vector = plv_weight * plv_direction * direction_vector;

        % Accumulate vector at electrode position
        vector_field{i} = vector_field{i} + weighted_vector;
    end
end

disp('Done with Spatial Synchrony Field');

% Compute average synchrony (mean PLV magnitude for each channel)
average_sync = mean(abs(PLV_vector_matrix), 2);
% Normalize for color scaling (between 0 and 1)
% average_sync = (average_sync - min(average_sync)) / (max(average_sync) - min(average_sync));

% Extract components for quiver plot (center vectors in grid boxes)
X = cellfun(@(p) p(1) + 0.5, electrode_positions);
Y = cellfun(@(p) p(2) + 0.5, electrode_positions);
U = cellfun(@(v) v(1), vector_field);
V = cellfun(@(v) v(2), vector_field);

% Plot the 21x3 Grid Properly with Blue-to-Yellow Heatmap
figure; hold on;
% Define colormap with 256 colors
cmap = parula(256); 

% Draw heatmap-colored boxes
for i = 1:num_channels
    x_i = electrode_positions{i}(1);
    y_i = electrode_positions{i}(2);
    
    % Get color based on average_sync
    color_idx = round(average_sync(i) * 255) + 1; % Ensure index is valid (1 to 256)
    color_idx = min(max(color_idx, 1), 256); % Clamp values within valid range
    
    % Draw a filled rectangle (box) colored by synchrony
    rectangle('Position', [x_i, y_i, 1, 1], ...
              'FaceColor', cmap(color_idx, :), 'EdgeColor', 'k');
end

% Plot vector field with strong contrast
quiver(X, Y, U, V, 0.8, 'r', 'LineWidth', 2, 'MaxHeadSize', 2); % Keep arrows red


% Label each box with electrode number
arrayfun(@(i) text(X(i), Y(i), num2str(i), 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k'), 1:num_channels);

% Formatting
axis equal; xlim([0, grid_width]); ylim([0, grid_height]);
set(gca, 'XTick', [], 'YTick', []); % Remove x-axis and y-axis ticks
xlabel(''); ylabel(''); % No axis labels
title('Spatial Synchrony Map');

% Reverse Y-direction to match bottom-left first channel
set(gca, 'YDir', 'normal');  

% Add Colorbar
colormap(cmap);
cbar = colorbar;
ylabel(cbar, 'Average PLV');

hold off;


ref_channel = reference_ch; % Reference channel PDN008: 36 - PrG, Left Precentral Gyrus A4ul, area 4 (upper limb region)
avg_phase_diff = zeros(num_channels, 1);
plv_with_ref = zeros(num_channels, 1);

for ch = 1:num_channels
    phase_diff = angle(exp(1j * (phaseData{ch} - phaseData{ref_channel}))); % circular diff
    avg_phase_diff(ch) = mean(phase_diff); % Mean circular phase difference
    plv_with_ref(ch) = abs(mean(exp(1j * phase_diff))); % PLV with channel 33
end

%plv_threshold = 0.85;
%best_neighbors = plv_with_ref > plv_threshold;

% Find best neighbor (excluding self)
plv_with_ref(ref_channel) = -1; % Exclude the reference channel from being its own best neighbor
[max_plv, best_neighbor_ch] = max(plv_with_ref); % Find channel with highest PLV

% Get its average phase difference
best_neighbor_phase_diff = avg_phase_diff(best_neighbor_ch);

% Display results
fprintf('Best neighbor of channel %d is channel %d\n', ref_channel, best_neighbor_ch);
fprintf('Highest PLV: %.3f\n', max_plv);
fprintf('Phase difference with channel %d: %.4f radians\n', ref_channel, best_neighbor_phase_diff);


figure; hold on;

% Use parula colormap for PLV
cmap = parula(256);

% Normalize PLV values for color mapping
plv_norm = plv_with_ref;

for i = 1:num_channels
    x_i = electrode_positions{i}(1);
    y_i = electrode_positions{i}(2);

    % Determine color from normalized PLV
    color_idx = round(plv_norm(i) * 255) + 1;
    color_idx = min(max(color_idx, 1), 256);

    % Draw rectangle
    rectangle('Position', [x_i, y_i, 1, 1], ...
              'FaceColor', cmap(color_idx, :), 'EdgeColor', 'k');

    % --- Label string with both channel name and phase value ---
    chan_name = string(channelNames{i});  % safer than char if cellstr
    phase_val = avg_phase_diff(i);
    label_str = sprintf('%s: %.3f', chan_name, phase_val);  % reliable formatting

    text(x_i + 0.5, y_i + 0.65, chan_name, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');

    text(x_i + 0.5, y_i + 0.35, sprintf('%.3f', phase_val), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 6, 'FontWeight', 'normal', 'Color', 'k');


    % Highlight best neighbor
    if i == best_neighbor_ch
        rectangle('Position', [x_i, y_i, 1, 1], ...
                  'EdgeColor', 'r', 'LineWidth', 3);
    end

    if i == reference_ch
        rectangle('Position', [x_i, y_i, 1, 1], ...
                  'EdgeColor', [0.4 0.8 1], ...  % light blue RGB
                  'LineWidth', 2);
    end
end



% Final plot settings
axis equal; xlim([0, grid_width]); ylim([0, grid_height]);
set(gca, 'XTick', [], 'YTick', []);
xlabel(''); ylabel('');
title('Avg Phase Diff (vs Ch 33) + PLV Heatmap');
set(gca, 'YDir', 'normal');

% Colorbar for PLV
colormap(cmap);
cbar = colorbar;
ylabel(cbar, 'PLV with Channel 33');

hold off;


