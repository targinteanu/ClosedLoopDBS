%PD23N003 - 'pds002.ns2'
%PD22N009 - 'datafile001.ns2'
%PD22N008 - 'PD22N008.ns2'
%PD23N009 - 'PDS001.ns2

% Inputs for both conditions
patient_data_cond1 = 'datafile001.ns2';  % Condition 1
patient_data_cond2 = 'PDS001.ns2';       % Condition 2
freq_range = 'beta';                     % 'beta' or 'theta'
recollect_data = true;


if recollect_data
    average_sync_1 = compute_average_sync_map(patient_data_cond1, freq_range);
    average_sync_2 = compute_average_sync_map(patient_data_cond2, freq_range);
end

% Define grid layout and positions
grid_width = 3;
grid_height = 21;
num_channels = grid_width * grid_height;
electrode_positions = cell(num_channels, 1);
channel_num = 1;
for x = 0:grid_width-1
    for y = 0:grid_height-1
        electrode_positions{channel_num} = [x, y];
        channel_num = channel_num + 1;
    end
end

% Plotting both synchrony maps side by side
figure;
cmap = parula(256); 

for condition = 1:2
    subplot(1,2,condition); hold on;
    average_sync = eval(sprintf('average_sync_%d', condition));

    for i = 1:num_channels
        pos = electrode_positions{i};
        color_idx = round(average_sync(i) * 255) + 1;
        color_idx = min(max(color_idx, 1), 256);
        rectangle('Position', [pos(1), pos(2), 1, 1], ...
                  'FaceColor', cmap(color_idx,:), 'EdgeColor', 'k');
    end

    % Label each box with electrode number
    arrayfun(@(i) text(X(i), Y(i), num2str(i), 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k'), 1:num_channels);

    axis equal; xlim([0, grid_width]); ylim([0, grid_height]);
    set(gca, 'XTick', [], 'YTick', []);
    title(sprintf('Condition %d Synchrony Map', condition));
end

% Add shared colorbar
h = colorbar('Position', [0.93 0.11 0.02 0.8]); % Adjust position for side colorbar
ylabel(h, 'Average PLV');
colormap(cmap);



%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to compute average synchrony map
function average_sync = compute_average_sync_map(patient_data, freq_range)
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
        channelName, channelIndex, channelIndexStim, channelNames] = ...
        getRecordedData_NS(patient_data, 1);

    numChannels = size(dataAllChannels, 1) - 1;
    channelIndices = 1:numChannels;

    % Store phases
    phaseData = cell(1, numChannels);
    for idx = 1:numChannels
        ch_raw = dataAllChannels(idx, :);
        ch_raw_cut = ch_raw(1:round(length(ch_raw)/3));
        if strcmp(freq_range, 'theta')
            ch_filtered = Myeegfilt(ch_raw_cut, SamplingFreq, 4, 9, 0, 1024);
        elseif strcmp(freq_range, 'beta')
            ch_filtered = Myeegfilt(ch_raw_cut, SamplingFreq, 13, 30, 0, 1024);
        else
            error('Unknown frequency range: %s', freq_range);
        end
        [ch_phase, ~] = instPhaseFreq(ch_filtered, SamplingFreq);
        phaseData{idx} = ch_phase;
    end

    % Compute PLV
    PLV_vector_matrix = zeros(numChannels, numChannels);
    for ch1 = 1:numChannels
        for ch2 = ch1+1:numChannels
            phase_diff = phaseData{ch1} - phaseData{ch2};
            PLV_vector_matrix(ch1, ch2) = mean(exp(1j * phase_diff));
            PLV_vector_matrix(ch2, ch1) = conj(PLV_vector_matrix(ch1, ch2));
        end
    end

    % Compute average synchrony
    average_sync = mean(abs(PLV_vector_matrix), 2);
end
