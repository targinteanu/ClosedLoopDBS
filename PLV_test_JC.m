% Define parameters
fs = 1000;  % Sampling frequency (Hz)
T = 1;      % Duration (seconds)
t = linspace(0, T, fs);  % Time vector

% Generate two sinusoids with a 45-degree (π/4 radians) phase difference
freq = 10;  % Frequency of sinusoids in Hz
phase_diff = pi / 4;  % 45 degrees in radians

sin1 = sin(2 * pi * freq * t);  % First sinusoid
sin2 = sin(2 * pi * freq * t + phase_diff);  % Second sinusoid

[phase1, ~] = instPhaseFreq(sin1, fs);
[phase2, ~] = instPhaseFreq(sin2, fs);

% Compute phase difference
phase_diff_signal = phase1 - phase2;

% Compute PLV
PLV = abs(mean(exp(1j * phase_diff_signal))); %The abs is equivalent to an L2 norm in this case apparently

% Display the computed PLV
fprintf('Computed PLV: %.4f\n', PLV);

% Plot the two sinusoids
figure;
plot(t, sin1, 'b', 'LineWidth', 1); hold on;
plot(t, sin2, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Two Sinusoids with a 45° Phase Difference');
legend('Sinusoid 1', 'Sinusoid 2 (45° phase shift)');
grid on;

% Define parameters
fs = 1000;  % Sampling frequency (Hz)
T = 1;      % Duration (seconds)
t = linspace(0, T, fs);  % Time vector

% Generate a base sinusoid
freq = 10;  % Frequency in Hz
base_phase_diff = pi / 4;  % 45 degrees in radians (initial phase difference)

% Add random phase noise that varies over time
random_phase_noise = pi/6 * randn(size(t)); % Phase noise (Gaussian distributed)

% Generate two signals
sin1 = sin(2 * pi * freq * t);  % First sinusoid (clean)
sin2 = sin(2 * pi * freq * t + base_phase_diff + random_phase_noise); % Second sinusoid with random phase variations

% Compute instantaneous phase using Hilbert Transform
analytic_signal1 = hilbert(sin1);
analytic_signal2 = hilbert(sin2);

phase1 = angle(analytic_signal1);
phase2 = angle(analytic_signal2);

% Compute phase difference
phase_diff_signal = phase1 - phase2;

% Compute PLV (should be < 1)
PLV = abs(mean(exp(1j * phase_diff_signal)));

% Display the computed PLV
fprintf('Computed PLV (with phase noise): %.4f\n', PLV);

% Plot the two signals
figure;
plot(t, sin1, 'b', 'LineWidth', 1); hold on;
plot(t, sin2, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Two Sinusoids with Time-Varying Phase Difference');
legend('Clean Sinusoid', 'Sinusoid with Phase Noise');
grid on;






