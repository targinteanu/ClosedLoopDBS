patient_data_cond = 'datafile001.ns2';  %PD22N009 - 'datafile001.ns2'
% patient_data_cond2 = 'PDS001.ns2';    
freq_range = 'beta';   % 'beta' or 'theta'

recollect_artifactRemoved = false;
recollect_avg_sync = false;

nev = load('PD22N009_nev.mat');

% Get the timestamps and event codes
timestamps_sec = nev.NEV.Data.SerialDigitalIO.TimeStampSec;
events = nev.NEV.Data.SerialDigitalIO.UnparsedData;  

baseline_end = 1202530/10;
baseline_interval = {[1 baseline_end]};

condition_intervals = {[2050000 2170675],[2215000 2302550],[2350000 2411100]};


threshold = 0.1;
AR_order = 10;


if recollect_artifactRemoved
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
            channelName, channelIndex, channelIndexStim, channelNames] = ...
            getRecordedData_NS(patient_data_cond, 1);
    
    numChannels = size(dataAllChannels, 1) - 1;
    channelIndices = 1:numChannels;
    
    patient_data = dataAllChannels;
    dataBaseline = dataAllChannels(channelIndices, 1:baseline_end);
    
    patient_data = patient_data(channelIndices, :);
    dataBaseline = dataBaseline(channelIndices, :); 
    patient_data_cond1 = dataBaseline;
    
    % downsample_factor = 20;
    % dataBaseline = dataBaseline(channelIndices, 1:downsample_factor:end);
    
    stim_points = detectStimPoints(patient_data, dataBaseline, threshold);
    patient_data_artifactRemoved = removeArtifactsAR(patient_data, stim_points, dataBaseline, AR_order);
    patient_data_cond2 = patient_data_artifactRemoved;
    
    
end

plot_channel = 59;

original = patient_data(plot_channel, :);
cleaned  = patient_data_artifactRemoved(plot_channel, :);
t = (0:length(original)-1) / 1000;  % Assuming fs = 1000 Hz

figure;
plot(t, original, 'r', 'DisplayName', 'Original'); hold on;
plot(t, cleaned,  'b', 'DisplayName', 'Cleaned');
xlabel('Time (s)');
ylabel('Amplitude');
title(['Full Signal - Channel ', num2str(plot_channel)]);
legend;
grid on;

patient_data_cond2 = patient_data_artifactRemoved;
if recollect_avg_sync
    average_sync_1 = compute_average_sync_map(patient_data_cond1, freq_range, baseline_interval);
    average_sync_2 = compute_average_sync_map(patient_data_cond2, freq_range, condition_intervals);
end


disp(mean(average_sync_1));
disp(mean(average_sync_2));


avgPLV_conditions = {
    average_sync_1;
    average_sync_2
};

figure;
plotPLVGrids(avgPLV_conditions);












function plotPLVGrids(avgPLV_conditions)
% plotPLVGrid - Visualize PLV maps for each condition on a 3x21 grid
%
% Inputs:
%   avgPLV_conditions - cell array of [63 x 1] vectors, one per condition

    num_conditions = length(avgPLV_conditions);

    % Define grid layout
    grid_width = 3;
    grid_height = 21;
    channel_order = reshape(63:-1:1, [21, 3]);  % 21 rows, 3 columns

    for idx = 1:num_conditions
        avgPLV_per_channel = avgPLV_conditions{idx};

        % Create 21x3 grid from linear channel vector
        gridPLV = zeros(21, 3);
        for row = 1:21
            for col = 1:3
                ch = channel_order(row, col);
                gridPLV(row, col) = avgPLV_per_channel(ch);
            end
        end

        % Plotting
        subplot(1, num_conditions, idx);
        imagesc(gridPLV);
        axis equal tight;
        colormap('parula');
        colorbar;
        caxis([0 1]);
        title(['Average PLV per Channel - Condition ' num2str(idx)]);
        set(gca, 'XTick', [], 'YTick', []);

        % Overlay channel numbers
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
end











function patient_data_artifactRemoved = removeArtifactsAR(patient_data, stim_points, dataBaseline, AR_order)
% removeArtifactsAR - Remove stimulation artifacts using AR model
%
% Inputs:
%   patient_data   - matrix of neural data (channels x timepoints)
%   stim_points    - vector of time indices (stim onsets) to remove artifact from
%   dataBaseline   - matrix of baseline data (channels x timepoints) for AR training
%   AR_order       - scalar, order of AR model
%
% Output:
%   patient_data_artifactRemoved - artifact-corrected data

    fs = 1000;               % Sampling rate (Hz)
    artdur = 0.013;           % Duration of artifact (seconds)
    artdur_samples = ceil(artdur * fs);
    context_len = max(AR_order, 100);  

    numChannels = size(patient_data, 1);
    numTimepoints = size(patient_data, 2);

    % Initialize output
    patient_data_artifactRemoved = patient_data;

    % === Precompute AR models ===
    ARmdls = cell(1, numChannels);
    for ch = 1:numChannels
        fprintf('Precomputing AR for CH %d / %d\n', ch, size(patient_data, 1));
        ARmdls{ch} = ar(dataBaseline(ch, :)', AR_order, 'yw');
    end


    for ch = 1:numChannels
        fprintf('Removing artifact from CH %d / %d\n', ch, size(patient_data, 1));

        ARmdl = ARmdls{ch};

        for i = 1:length(stim_points)
            t1 = stim_points(i);
            t2 = min(t1 + ceil(0.95 * artdur_samples), numTimepoints);
            t0 = max(1, t2 - artdur_samples);
            K = t2 - t0;

            if K > 0
                % Use only recent samples for AR forecasting
                start_idx = max(1, t0 - context_len + 1);
                AR_input = patient_data_artifactRemoved(ch, start_idx:t0)';
                ARpred = myFastForecastAR(ARmdl, AR_input, K);

                % Replace artifact with forecast
                patient_data_artifactRemoved(ch, (t0+1):t2) = ARpred';
            end
        end
    end
end






function stim_points = detectStimPoints(patient_data, baseline_data, threshold)
% detectStimPoints - Identify timepoints of stimulation artifacts
%
% Inputs:
%   patient_data   - neural data matrix (channels x timepoints)
%   baseline_data  - baseline segment of data (channels x timepoints)
%   threshold      - fraction of channels that must exceed threshold (e.g., 0.5)
%
% Output:
%   stim_points - vector of time indices where artifact is detected

    % === Derivative across time ===
    d_data = diff(patient_data, 1, 2);            % channels x (timepoints - 1)
    d_data_abs = abs(d_data);

    % === Derivative of baseline ===
    d_baseline = diff(baseline_data, 1, 2);
    d_baseline_abs = abs(d_baseline);

    % === Baseline stats per channel ===
    baseline_mean = mean(d_baseline_abs, 2);      % [channels x 1]
    baseline_std  = std(d_baseline_abs, 0, 2);    % [channels x 1]

    stim_mask = d_data_abs > (baseline_mean + 10 * baseline_std);  % logical mask

    % === Proportion of channels above threshold per timepoint ===
    stim_proportion = mean(stim_mask, 1);         % 1 x (timepoints - 1)

    % === Timepoints where > threshold fraction of channels fire ===
    stim_points = find(stim_proportion > threshold)-1;  % shift for diff()
end




function average_sync = compute_average_sync_map(data, freq_range, condition_intervals)
% compute_average_sync_map - Compute average phase synchrony within given intervals
%
% Inputs:
%   dataAllChannels      - neural data [channels x timepoints]
%   freq_range           - 'theta' or 'beta'
%   condition_intervals  - cell array of [start end] index pairs (e.g., {[s1 e1], [s2 e2], ...})
%
% Output:
%   average_sync         - [numChannels x 1] vector of mean PLV per channel

    numChannels = size(data, 1);
    channelIndices = 1:numChannels;
    SamplingFreq = 1000;

    % === Step 1: Extract phase data per channel ===
    phaseData = cell(1, numChannels);
    for idx = 1:numChannels
        fprintf('Grabbing phase for CH %d / %d\n', idx, numChannels);

        ch_raw = data(idx, :);
        if strcmp(freq_range, 'theta')
            ch_filtered = Myeegfilt(ch_raw, SamplingFreq, 4, 9, 0, 1024);
        elseif strcmp(freq_range, 'beta')
            ch_filtered = Myeegfilt(ch_raw, SamplingFreq, 13, 30, 0, 1024);
        else
            error('Unknown frequency range: %s', freq_range);
        end
        [ch_phase, ~] = instPhaseFreq(ch_filtered, SamplingFreq);
        phaseData{idx} = ch_phase;
    end

    % === Step 2: Compute PLVs within condition intervals ===
    PLV_vector_matrix = zeros(numChannels, numChannels);
    for ch1 = 1:numChannels
        for ch2 = ch1+1:numChannels
            plv_vals = [];
            for k = 1:length(condition_intervals)
                interval = condition_intervals{k};

                i_start = interval(1);
                i_end = interval(2);
            
                phase_diff = phaseData{ch1}(i_start:i_end) - phaseData{ch2}(i_start:i_end);
                plv_vals(end+1) = mean(exp(1j * phase_diff));
            end


            % Average across all intervals
            PLV_vector_matrix(ch1, ch2) = mean(plv_vals);
            PLV_vector_matrix(ch2, ch1) = conj(PLV_vector_matrix(ch1, ch2));
        end
    end

    % === Step 3: Compute per-channel average sync ===
    average_sync = mean(abs(PLV_vector_matrix), 2);
end




