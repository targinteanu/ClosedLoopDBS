channelIndices = 1:64;
numChannels = length(channelIndices);
phaseData = cell(1, numChannels);  % Store phase data for each channel


[dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS('PD22N008.ns2',1);


for idx = 1:numChannels
    
    % Raw signal
    ch_raw = dataAllChannels(idx, :);
    % 
    % % Filter to beta band (13-30 Hz)
    ch_filtered = Myeegfilt(ch_raw, SamplingFreq, 13, 30, 0, 1024);

    % Compute phase
    [ch_phase, ~] = instPhaseFreq(ch_filtered, SamplingFreq);
    
    % Store phase data
    phaseData{idx} = ch_phase;
    disp(idx)
end

disp('done with phases')


%% Compute synchrony (average phase difference)
synchronyMatrix = zeros(numChannels);

for i = 1:numChannels
    for j = 1:numChannels
        if i ~= j
            phaseDiff = abs(phaseData{i} - phaseData{j});  % Note: Without taking absolute value we incorrectly map phase difference linearly. Can solve with abs
            synchronyMatrix(i, j) = mean(phaseDiff);       % Average phase difference
            disp(mean(phaseDiff))
        else
            synchronyMatrix(i, j) = 0;  
        end
    end
end

disp('done with phase differences')

%% Plot synchrony matrix
figure;
imagesc(synchronyMatrix);
colorbar;
%clim([-pi, pi]); % Normalize color scale (0 to π, since phase difference is periodic)
xticks(1:numChannels);
yticks(1:numChannels);
xticklabels(channelIndices);
yticklabels(channelIndices);
title('Channel Synchrony Matrix');
xlabel('Channel');
ylabel('Channel');


startTime = t(20000);                  
numSamples = round(SamplingFreq);     
[~, startIdx] = min(abs(t - startTime));
endIdx = min(startIdx + numSamples - 1, length(t));

endingIndex = 250000;
t_subset = t(1:endingIndex);

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

disp('done with PLVs')

% Assuming PLV_matrix is already computed
figure;
imagesc(PLV_value_matrix); % Display the PLV matrix
colorbar; % Add a color scale
colormap('parula'); % Replace 'viridis' with 'parula' or any other built-in colormap
caxis([0 1]); % Set color scale between 0 and 1 for better interpretation

% Label the axes with channel numbers
numChannels = size(PLV_value_matrix, 1);
xticks(1:numChannels);
yticks(1:numChannels);
xticklabels(1:numChannels); % Just numbers, no "Ch"
yticklabels(1:numChannels);

% Rotate x-axis labels for better readability
xtickangle(90);

% Add axis labels and title
xlabel('Channel');
ylabel('Channel');
title('Phase Locking Value (PLV) Matrix');



% Step 1: Find the two channels with the highest PLV magnitude
PLV_magnitudes = abs(PLV_vector_matrix); % Extract PLV magnitudes
PLV_magnitudes(logical(eye(numChannels))) = 0; % Set diagonal to 0 (ignore self-comparison)

% Find the indices of the max PLV value
[maxPLV, maxIdx] = max(PLV_magnitudes(:)); % Get maximum PLV value
[ch1, ch2] = ind2sub(size(PLV_magnitudes), maxIdx); % Convert linear index to channel indices

% Step 2: Compute phase differences for these two channels
phase_diff = angle(exp(1j * (phaseData{ch1} - phaseData{ch2}))); % Correctly wrapped phase differences

% Step 3: Plot the Rose Plot (Polar Histogram)
figure;
polarhistogram(phase_diff, 'BinEdges', linspace(-pi, pi, 100)); % 20 bins from -pi to pi
title(sprintf('Phase Difference Distribution: Ch%d vs Ch%d', ch1, ch2));



% Define grid dimensions
grid_width = 3;
grid_height = 21;
num_channels = grid_width * grid_height;

% Generate electrode coordinates using a 63x1 cell array
electrode_positions = cell(num_channels, 1);
channel_num = 1;
for x = 0:grid_width-1
    for y = 0:grid_height-1
        electrode_positions{channel_num} = [x, y]; % Store as a 1x2 array
        channel_num = channel_num + 1;
    end
end

% Initialize vector field as a 63x1 cell array
vector_field = cell(num_channels, 1);
for i = 1:num_channels
    vector_field{i} = [0, 0]; % Initialize each as a 1x2 vector
end

% Compute phase gradient vectors
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

disp('done with traveling wave vector field');

% Compute average synchrony (mean PLV magnitude for each channel)
average_sync = mean(abs(PLV_vector_matrix), 2);

% Normalize for color scaling (between 0 and 1)
average_sync = (average_sync - min(average_sync)) / (max(average_sync) - min(average_sync));

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
    color_idx = round(color_val * 255) + 1; % Ensure index is valid (1 to 256)
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
set(gca, 'YDir', 'reverse'); grid off;
xlabel('Electrode Column'); ylabel('Electrode Row');
title('Phase Lag Vector Field with Blue-to-Yellow Synchrony Heatmap');
hold off;



%%%%%% Traveling Wave Vector Field %%%%%%%%%%%

% % Define grid dimensions
% grid_width = 3;
% grid_height = 21;
% num_channels = grid_width * grid_height;
% 
% % Generate electrode coordinates using a 63x1 cell array
% electrode_positions = cell(num_channels, 1);
% channel_num = 1;
% for x = 0:grid_width-1
%     for y = 0:grid_height-1
%         electrode_positions{channel_num} = [x, y]; % Store as a 1x2 array
%         channel_num = channel_num + 1;
%     end
% end
% 
% % Initialize vector field as a 63x1 cell array
% vector_field = cell(num_channels, 1);
% for i = 1:num_channels
%     vector_field{i} = [0, 0]; % Initialize each as a 1x2 vector
% end
% 
% % Compute phase gradient vectors
% for i = 1:num_channels
%     x_i = electrode_positions{i}(1);
%     y_i = electrode_positions{i}(2);
% 
%     for j = 1:num_channels
%         if i == j
%             continue; % Skip self-comparison
%         end
% 
%         x_j = electrode_positions{j}(1);
%         y_j = electrode_positions{j}(2);
% 
%         % Compute direction vector
%         direction_vector = [x_j - x_i, y_j - y_i];
%         norm_factor = norm(direction_vector);
%         if norm_factor == 0
%             continue; % Avoid division by zero
%         end
%         direction_vector = direction_vector / norm_factor; % Normalize
% 
%         % Get phase difference and PLV
%         plv_direction = sign(angle(PLV_vector_matrix(i,j)));
%         %print(plv_direction);
%         plv_weight = PLV_value_matrix(i, j);
% 
%         % Compute weighted phase gradient vector
%         weighted_vector = plv_weight * plv_direction * direction_vector;
% 
%         % Accumulate vector at electrode position
%         vector_field{i} = vector_field{i} + weighted_vector;
%     end
% end
% 
% disp('done with traveling wave vector field')
% 
% 
% 
% % Extract components for quiver plot (center vectors in grid boxes)
% X = cellfun(@(p) p(1) + 0.5, electrode_positions);
% Y = cellfun(@(p) p(2) + 0.5, electrode_positions);
% U = cellfun(@(v) v(1), vector_field);
% V = cellfun(@(v) v(2), vector_field);
% 
% % Plot the 21x3 Grid Properly
% figure; hold on;
% arrayfun(@(x) plot([x x], [0 grid_height], 'k', 'LineWidth', 1), 0:grid_width); % Vertical lines
% arrayfun(@(y) plot([0 grid_width], [y y], 'k', 'LineWidth', 1), 0:grid_height); % Horizontal lines
% 
% % Plot vector field with improved visibility
% quiver(X, Y, U, V, 0.8, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
% 
% % Label each box with electrode number
% arrayfun(@(i) text(X(i), Y(i), num2str(i), 'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k'), 1:num_channels);
% 
% % Formatting
% axis equal; xlim([0, grid_width]); ylim([0, grid_height]);
% set(gca, 'YDir', 'reverse'); grid off;
% xlabel('Electrode Column'); ylabel('Electrode Row');
% title('Phase Lag Vector Field');
% hold off;
% 














% t_numeric = second(t) + minute(t) * 60 + hour(t) * 3600;
% 
% % Compute phase difference derivative synchrony matrix
% syncDerivativeMatrix = zeros(numChannels);
% 
% for i = 1:numChannels
%     for j = 1:numChannels
%         if i ~= j
%             % Compute phase difference
%             phaseDiff = phaseData{i} - phaseData{j};
% 
%             % Compute time derivative using gradient
%             phaseDiffDerivative = gradient(phaseDiff, t_numeric); 
% 
%             % Compute average absolute value of derivative
%             syncDerivativeMatrix(i, j) = mean(abs(phaseDiffDerivative));
%         else
%             syncDerivativeMatrix(i, j) = NaN; % NaN for same channel comparison
%         end
%     end
% end
% 
% % Plot synchrony derivative matrix
% figure;
% imagesc(syncDerivativeMatrix);
% colorbar;
% colormap hot;
% caxis([0, max(syncDerivativeMatrix(:))]); % Normalize color scale
% xticks(1:numChannels);
% yticks(1:numChannels);
% xticklabels(channelIndices);
% yticklabels(channelIndices);
% title('Phase Difference Derivative Synchrony');
% xlabel('Channel');
% ylabel('Channel');
