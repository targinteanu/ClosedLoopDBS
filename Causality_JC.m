clear; close all; clc;

% Step 1: Generate Simulated Multichannel Data
fs = 100;   % Sampling frequency in Hz
t = 0:1/fs:10;  % 10 seconds of data
n = length(t);

% Generate three synthetic causal signals
channel1 = sin(2*pi*10*t) + 0.2*randn(1,n);   % Base signal
channel2 = 0.8 * [0, channel1(1:end-1)] + 0.2*randn(1,n);  % Ch1 → Ch2
channel3 = 0.5 * [0, channel2(1:end-1)] + 0.3*randn(1,n);  % Ch2 → Ch3

% Stack signals into matrix (Rows = Channels, Columns = Time Points)
data = [channel1; channel2; channel3];

% Step 2: Fit a Multivariate Autoregressive (MVAR) Model
p = 10; % Model order (choose via BIC/AIC if needed)
[A_flat, SIG] = arfit(data', p, p);  % AR coefficients and residuals

% Step 3: Reshape AR Coefficients for AsympPDC
nChannels = size(data,1);  % Number of channels

% ✅ Correctly reshape the AR coefficient matrix
A = reshape(A_flat, nChannels, nChannels, p); % Step 1: Reshape to (nChannels, nChannels, p)
A = permute(A, [1, 2, 3]); % Step 2: Keep 3D format for `asymp_pdc`

% Step 4: Compute Partial Directed Coherence (PDC)
nFreqs = 128;    % Number of frequency bins
metric = 'euc';  % Use 'euc' for standard PDC, 'diag' for gPDC, 'info' for iPDC
alpha = 0.05;    % Significance level

% ✅ Correct function call to match asymp_pdc's required format
c = asymp_pdc(data, A, SIG, nFreqs, metric, alpha);

% Step 5: Extract and Visualize PDC
freqs = linspace(0, fs/2, nFreqs); % Generate frequency axis

% Extract PDC at 10 Hz
freq_idx = find(abs(freqs - 10) == min(abs(freqs - 10))); % Closest index to 10 Hz
pdc_10Hz = c.pdc2(:,:,freq_idx); % Extract squared PDC

% Step 6: Display and Visualize PDC at 10 Hz
disp('PDC Matrix at 10 Hz:');
disp(pdc_10Hz);

figure;
imagesc(pdc_10Hz);
colorbar;
xticklabels({'Ch1', 'Ch2', 'Ch3'});
yticklabels({'Ch1', 'Ch2', 'Ch3'});
title('PDC at 10 Hz');
xlabel('Influenced Channel');
ylabel('Source Channel');
