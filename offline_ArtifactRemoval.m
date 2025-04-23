%% offline_ArtifactRemoval 
% Removes both cortical stimulation and DBS artifacts from all channels of
% data recorded as a BlackRock ns2 file

%% parameters to set 

num_channel_threshold = 0.1;
% An outlier value in at least this proportion of channels simultaneously
% will be considered the onset of a DBS pulse. 

AR_order = 10;
% order of AR model 

baseline_start = 1; % sample
baseline_end = 120000; % sample 
% Time index of the start/end of baseline (no stim) period to use for AR 
% model fitting. First 120 seconds are used for subject PD22N009

%% load data 

filename = '';
% This can be empty or point to an ns2 or mat file. If left empty, a dialog
% will prompt the user to select a file. It is assumed that cortical
% stimulation pulses are logged on channel 'ainp1'.

[~, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
        ~, ~, channelIndStim, channelNames] = ...
        getRecordedData_NS(filename, 1);

% remove the stimulation pulse log channel 
dataAllChannels = dataAllChannels(...
    [1:(channelIndStim-1),(channelIndStim+1):(height(dataAllChannels))], :);

%% processing and artifact removal 

numChannels = size(dataAllChannels, 1) - 1;
channelIndices = 1:numChannels; % can adjust this to exclude any channels

patient_data = dataAllChannels;
dataBaseline = dataAllChannels(channelIndices, baseline_start:baseline_end);

patient_data = patient_data(channelIndices, :);
dataBaseline = dataBaseline(channelIndices, :);

% get stim points of DBS 
stim_points = detectStimPoints(patient_data, dataBaseline, num_channel_threshold);

% add in stim points on cortex 
stim_points = [find(StimTrainRec), stim_points]; 
stim_points = sort(unique(stim_points));

% remove artifact
patient_data_artifactRemoved = removeArtifactsAR(patient_data, stim_points, ...
    dataBaseline, AR_order, SamplingFreq);

%% display results 

plot_channel = 59;

original = patient_data(plot_channel, :);
cleaned  = patient_data_artifactRemoved(plot_channel, :);

figure;
plot(t, original, 'r', 'DisplayName', 'Original'); hold on;
plot(t, cleaned,  'b', 'DisplayName', 'Cleaned');
xlabel('Time (s)');
ylabel(channelNames{plot_channel});
title(['Full Signal - Channel ', num2str(plot_channel)]);
legend;
grid on;

%% internal helper functions 
% TO DO: consider making some of these into external functions 

function stim_points = detectStimPoints(patient_data, baseline_data, threshold)
% detectStimPoints - Identify timepoints of stimulation artifacts
% original by Jack in testing4_15_JC.m / edited by Toren 
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
    stim_points = find(stim_proportion > threshold) + 1;  % shift for diff()
end

function patient_data_artifactRemoved = removeArtifactsAR(patient_data, stim_points, dataBaseline, AR_order, fs)
% removeArtifactsAR - Remove stimulation artifacts using AR model
% original by Jack in testing4_15_JC.m / edited by Toren 
%
% Inputs:
%   patient_data   - matrix of neural data (channels x timepoints)
%   stim_points    - vector of time indices (stim onsets) to remove artifact from
%   dataBaseline   - matrix of baseline data (channels x timepoints) for AR training
%   AR_order       - scalar, order of AR model
%   fs             - % Sampling rate (Hz)
%
% Output:
%   patient_data_artifactRemoved - artifact-corrected data

    artdur = 0.02;           % Duration of artifact (seconds)
    artdur_samples = ceil(artdur * fs);
    context_len = max(AR_order, 100);  

    numChannels = size(patient_data, 1);
    numTimepoints = size(patient_data, 2);

    % Initialize output
    patient_data_artifactRemoved = patient_data;

    % === Precompute AR models ===
    ARmdls = cell(1, numChannels);
    for ch = 1:numChannels
        fprintf('Computing AR for CH %d / %d\n', ch, size(patient_data, 1));
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