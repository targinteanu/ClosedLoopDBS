
Patient_ID = 'PD25N006'; %PD25N009
Patient_Folder_Path = 'C:\Users\cours\Downloads\AndersonLabDataAndPackages\Patient Data';

NO_timetable_name = 'SavedTable_Distances_1375Hz1.mat';  
NO_to_BR_ch = 'CMacro_LFP_01'; %CMacro_LFP_AP_T
BR_to_NO_ch = 'AINP3';

ns_files = { ... %PD25N009 - 'baseline001.ns2', ...'DBS004.ns2','memory003.ns2'
    'datafile001.ns2', ... %%PD25N008 - 'baseline001.ns2', 'DepthStim002.ns2','nback004.ns2'
        
    }; 

target_window = [-3 3];  
target_channel = 'CMacro_LFP_01';
minDurSec = 20;    

% ----------- Load NO Data -----------
NO_timetable_path = fullfile(Patient_Folder_Path, Patient_ID, NO_timetable_name);
NO_timetable_OG = load(NO_timetable_path);

NO_timetable_unmodified = sortrows(NO_timetable_OG.Tbl, 'Time');
NO_SamplingFreq = 1375; 
NO_timetable = myRetime(NO_timetable_unmodified, NO_SamplingFreq, 0, 2);
    
NO_targetDist_unmodified = NO_timetable.("targetDistance");
NO_signal_unmodified = NO_timetable.(NO_to_BR_ch); % or other channel

% ----------- Load BR Data -----------
ns_timetables = cell(size(ns_files));
for k = 1:numel(ns_files)
    ns_path = fullfile(Patient_Folder_Path, Patient_ID, ns_files{k});
    disp(['Loading: ', ns_path])
    ns = openNSx(ns_path);
    if k == 1
        BR_SamplingFreq = ns.MetaTags.SamplingFreq;
    end
    
    tt = ns2timetable(ns_path);
    
    % add source file label
    tt.SourceFile = repmat(string(ns_files{k}), height(tt), 1);
    ns_timetables{k} = tt;
end

combined_tt = vertcat(ns_timetables{:});
combined_tt = sortrows(combined_tt, 'Time');
BR_timetable = myRetime(combined_tt, BR_SamplingFreq, 0, 2);


% ----------- Resample NO data to BR sampling rate -----------
t_NO = (0:height(NO_targetDist_unmodified)-1)' / NO_SamplingFreq;
duration = (height(NO_targetDist_unmodified)-1) / NO_SamplingFreq;
t_BR = (0:round(duration*BR_SamplingFreq))' / BR_SamplingFreq;

%%

% Interpolate WITHOUT filtering
NO_targetDist = interp1(t_NO, NO_targetDist_unmodified, t_BR, 'linear'); 

NO_signal     = resample(NO_signal_unmodified, BR_SamplingFreq, NO_SamplingFreq);

% Extract BR signal
BR_signal   = BR_timetable.(BR_to_NO_ch);

%Extract other NO channels
allVars   = NO_timetable.Properties.VariableNames;
sigVars   = setdiff(allVars, "targetDistance");
resampled_data = cellfun(@(v) resample(double(NO_timetable.(v)), ...
                                       BR_SamplingFreq, NO_SamplingFreq), ...
                         sigVars, 'UniformOutput', false);
NO_resampled_matrix = [resampled_data{:}];


%% Visualize All Data Side-by-Side

% --- Step 1: Ensure column vectors ---
BR_signal = double(BR_signal(:));
NO_signal = double(NO_signal(:));

% --- Step 2: Clip signals to [-1000, 2000] ---
BR_signal(BR_signal > 2000)  = 2000;
BR_signal(BR_signal < -2000) = -2000;
NO_signal(NO_signal > 2000)  = 2000;
NO_signal(NO_signal < -2000) = -2000;

% % --- Step 3: Demean (after clipping) ---
BR_signal = BR_signal - mean(BR_signal);
NO_signal = NO_signal - mean(NO_signal);

% --- Step 4: Time vectors ---
t_BR = (0:length(BR_signal)-1)' / BR_SamplingFreq;
t_NO = (0:length(NO_signal)-1)' / BR_SamplingFreq;
t_NO_dist = (0:length(NO_targetDist)-1)' / BR_SamplingFreq;

figure;
subplot(3,1,1);
plot(t_BR, BR_signal, 'b'); hold on;
xlabel('Time (s)'); ylabel('Amplitude');
title('BR\_signal (clipped + demeaned, pre-alignment)');
ylim([-2000 2000]);

% --- Mark file segments ---
segment_lengths = cellfun(@(tt) height(tt), ns_timetables);
segment_edges = [0 cumsum(segment_lengths)] / BR_SamplingFreq;

for k = 1:numel(ns_files)
    % draw vertical line at segment start (except very first point)
    if k > 1
        xline(segment_edges(k), '--r');
    end
    
    % midpoint for label
    mid_time = (segment_edges(k) + segment_edges(k+1)) / 2;
    
    % add filename label
    text(mid_time, 1800, ns_files{k}, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize',8, 'Color','k', 'BackgroundColor','w');
end


subplot(3,1,2);
plot(t_NO, NO_signal, 'r');
xlabel('Time (s)'); ylabel('Amplitude');
title('NO\_signal (clipped + demeaned, pre-alignment)');
ylim([-2000 2000]);

subplot(3,1,3);
plot(t_NO_dist, NO_targetDist, 'r');
xlabel('Time (s)'); ylabel('Amplitude');
title('NO\_targetDist (clipped + demeaned, pre-alignment)');
ylim([-5 25]);

%% First Cross Correlation and Alignment

fs = 1000;
br_win = [1800 2500];
no_win = [1700 2400];

[lag_s, NO_on_BR, TT_overlap] = align_NO_to_BR( ...
    BR_signal, NO_signal, fs, br_win, no_win, ...
    BR_timetable, NO_timetable, NO_resampled_matrix, ...
    sigVars, NO_targetDist, BR_to_NO_ch, NO_to_BR_ch, true);

fprintf('Lag = %.3f s\n', lag_s);
disp(TT_overlap)


%%

% --- full-length plot after alignment ---
figure;
hold on; grid on

% plot BR and NO signals
plot(t_BR_full, BR_timetable.(BR_to_NO_ch), 'b');
plot(NO_on_BR.Time, NO_on_BR.(NO_to_BR_ch), 'r');

% highlight BR segment used for xcorr
patch([t_BR_full(iBR_start) t_BR_full(iBR_end) ...
       t_BR_full(iBR_end) t_BR_full(iBR_start)], ...
      [-1000 -1000 1000 1000], ...
      [0 0 1], 'FaceAlpha',0.1, 'EdgeColor','none');

% highlight NO segment (after alignment) used for xcorr
patch([NO_on_BR.Time(iNO_start) NO_on_BR.Time(iNO_end) ...
       NO_on_BR.Time(iNO_end) NO_on_BR.Time(iNO_start)], ...
      [-1000 -1000 1000 1000], ...
      [1 0 0], 'FaceAlpha',0.1, 'EdgeColor','none');

xlabel('Time (s, BR clock)');
ylabel('Amplitude');
ylim([-1000 1000]);   % <-- force consistent y-limits
legend('BR','NO (aligned)','BR seg (xcorr)','NO seg (xcorr)');
title(sprintf('Full alignment via cross-corr (lag = %.3f s)', lag_s));




%%

fs = 1000;
t_common = (max(BR_timetable.Time(1), NO_on_BR.Time(1)) : ...
            seconds(1/fs) : ...
            min(BR_timetable.Time(end), NO_on_BR.Time(end)))';



% Synchronize only numeric variables
TT_master_num = synchronize(removevars(BR_timetable,"SourceFile"), NO_on_BR, t_common, 'linear');

% Resample the SourceFile column using nearest neighbor
TT_source = retime(BR_timetable(:, "SourceFile"), t_common, 'nearest');

% Merge back
TT_master = [TT_master_num TT_source];


t = TT_master.Properties.RowTimes;
figure;
plot(t, TT_master.(BR_to_NO_ch), 'b'); hold on
plot(t, TT_master.(NO_to_BR_ch), 'r');
legend('BR','NO'); grid on
xlabel('Time'); ylabel('Signal');
title('NO full span with BR zero-padded');
ylim([-1000 1000]);



%% Longest in-threshold NO segment (from TT_master)

fs = 1000;
sigName = 'targetDistance';  % column in TT_master
win = target_window;         % e.g. [-1 4]
win = [-5 3];

frontCut = 900;
backCut = 2525;

TT_longest = extractNthLongestWindow(TT_master, sigName, fs, win, frontCut, backCut, 1);

unique(TT_longest.SourceFile)

TT_longest_clean = removevars(TT_longest, {'SourceFile'});


%%
% Refine alignment inside TT_longest and re-sync (correct lag direction)

fs = 1000;  % TT_master/TT_longest are on BR clock

% --- Cross-correlate inside the current window ---
BR_seg = TT_longest.(BR_to_NO_ch);
NO_seg = TT_longest.(NO_to_BR_ch);

[c,lags] = xcorr(NO_seg - mean(NO_seg), BR_seg - mean(BR_seg), 'coeff');
[~,imax] = max(c);
refined_lag = lags(imax);          % samples
refined_lag_s = refined_lag/fs;    % seconds

fprintf('Refined lag in TT_longest: %d samples (%.6f s)\n', refined_lag, refined_lag_s);

% --- Shift ALL NO channels by the refined lag ---
% IMPORTANT: xcorr(NO,BR): positive lag => NO lags BR => move NO earlier => subtract
NO_on_BR_refined = NO_on_BR;
NO_on_BR_refined.Time = NO_on_BR.Time - seconds(refined_lag_s);

% --- Rebuild a synchronized master table on the refined alignment ---
t_common = (max(BR_timetable.Time(1), NO_on_BR_refined.Time(1)) : ...
            seconds(1/fs) : ...
            min(BR_timetable.Time(end), NO_on_BR_refined.Time(end)))';

TT_master_num = synchronize(removevars(BR_timetable,"SourceFile"), ...
                            NO_on_BR_refined, t_common, 'linear');
TT_source = retime(BR_timetable(:, "SourceFile"), t_common, 'nearest');
TT_master_refined = [TT_master_num TT_source];

% --- (Optional) re-extract longest usable window on refined alignment ---
TT_longest_refined = extractNthLongestWindow(TT_master_refined, ...
                                          'targetDistance', fs, target_window, frontCut, backCut,1);

% --- Quick visual check around the old window ---
figure; 
subplot(2,1,1); plot(lags/fs, c); grid on
xlabel('Lag (s)'); ylabel('r'); title('Refined XCorr in TT\_longest');

subplot(2,1,2); 
plot(TT_longest.Time, TT_longest.(BR_to_NO_ch), 'b'); hold on
plot(TT_longest.Time, TT_longest.(NO_to_BR_ch), 'r'); 
title('Before refined shift'); legend('BR','NO'); grid on; ylim([-1000 1000]);

figure; 
plot(TT_master_refined.Time, TT_master_refined.(BR_to_NO_ch), 'b'); hold on
plot(TT_master_refined.Time, TT_master_refined.(NO_to_BR_ch), 'r');
title('After refined shift (full span)'); legend('BR','NO (refined)'); grid on; ylim([-1000 1000]);



%%
signals = {"targetDistance","CMacro_LFP_01","AINP3","LS31"};
titles  = {"Target Distance (thresholds shown)", ...
           "NO signal: CMacro_LFP_AP_T", ...
           "BR signal: AINP2", ...
           "LS31"};
colors  = {'k','b','r','m'};   % added magenta for the NO_to_BR_ch

figure;
for i = 1:numel(signals)
    subplot(numel(signals),1,i)
    plot(TT_longest_refined.Time, TT_longest_refined.(signals{i}), colors{i});
    title(titles{i}); grid on;
    ylabel(signals{i}, 'Interpreter','none');
    if i == 1
        yline(win(1),'--r'); yline(win(2),'--r');
    end
    %%ylim([-1000 1000]); % optional if you want consistent scaling
end
xlabel('Time');

TT_longest_refined = removevars(TT_longest_refined, {'SourceFile'});


%%
[Phase_Data, Channel_Names] = obtainPhase_tt(TT_longest_refined, 1000, 'beta');
disp(Channel_Names)
disp('Obtained Phase Data')
PLV_VectorMatrix = compute_PLV_VectorMatrix(Phase_Data); %Get the PLV vector matrix for every channel (64x64)
PLV_Average = getAveragePLV(PLV_VectorMatrix); %Get the Average PLV for every channel (64x1)
disp('Calculated PLV Vector Matrix')
%%
visualize_PLV_Matrix(PLV_VectorMatrix, Channel_Names)


%%
fs = 1000;            % Hz (adjust if needed)
duration = 5;        % seconds
nSamp = duration * fs;

chA = 31; chB = 70;

t = (0:nSamp-1)/fs;

figure;
plot(t, Phase_Data{chA}(1:nSamp), 'b'); hold on;
plot(t, Phase_Data{chB}(1:nSamp), 'r');
xlabel('Time (s)');
ylabel('Phase (rad)');
title(sprintf('Phase of Ch %d and Ch %d (first %d s)', chA, chB, duration));
legend(sprintf('Ch %d', chA), sprintf('Ch %d', chB));
grid on;



%%
NO_to_BR_ch = 'CMacro_LFP_01'; %CMacro_LFP_AP_T AINP3
ref_idx = find(strcmp(Channel_Names, NO_to_BR_ch));
fprintf('"%s" is at channel index %d\n', Channel_Names{ref_idx}, ref_idx);
disp(ref_idx)
ref_channel = ref_idx;
PLV_Relative_ref = getPLV_wrt_Ref(PLV_VectorMatrix, ref_channel); %Get PLV relative to ref channel (64x1)
disp('Calculated PLV wrt Reference')

%%
electrode_data_name = 'electrode_data(Sheet1).csv';
electrode_data_filepath = fullfile(Patient_Folder_Path, Patient_ID, electrode_data_name);
disp(electrode_data_filepath)
meta = readtable(electrode_data_filepath, 'VariableNamingRule','preserve');


plot_PLV_Relative_Brainnetome(Phase_Data(1:63), Phase_Data{ref_channel}, PLV_Relative_ref, ref_channel, Channel_Names, meta);


%%

plot_PLV_Relative(Phase_Data(1:63),Phase_Data{ref_channel},PLV_Relative_ref, ref_channel, Channel_Names);  %Plot a 21x3 grid - color corresponds to PLV wrt ref channel, # corresponds to avg. phase difference wrt ref channel


%%

% Phase_Data = obtainPhase(ns_path, 'beta', 1, 100000); % (Baseline data 1-10,000) Get the phase data from every channel (
% Channel_Names = getChannelNames(ns_path);
% 
% disp('Obtained Phase Data')
% PLV_VectorMatrix = compute_PLV_VectorMatrix(Phase_Data); %Get the PLV vector matrix for every channel (64x64)
% PLV_Average = getAveragePLV(PLV_VectorMatrix); %Get the Average PLV for every channel (64x1)
% disp('Calculated PLV Vector Matrix')
% 
% ref_channel = 34;
% PLV_Relative_ref = getPLV_wrt_Ref(PLV_VectorMatrix, ref_channel); %Get PLV relative to ref channel (64x1)
% disp('Calculated PLV wrt Reference')
% 
% plot_PLV_Relative(Phase_Data,PLV_Relative_ref, ref_channel, Channel_Names);


%%

lh_pial_path = [Patient_Folder_Path Patient_ID '\freesurfer\surf\lh.pial'];
ACPC_Path = [Patient_Folder_Path Patient_ID '\Subject' Patient_ID '_elec_acpc_fr.mat'];
FSAVG_Path = [Patient_Folder_Path Patient_ID '\Subject' Patient_ID '_elec_fsavg_frs.mat'];

disp(lh_pial_path)
disp(ACPC_Path)

load_ACPC_FR_mesh(ACPC_Path, lh_pial_path, PLV_Relative_ref);
load_FSAVG_FRS_mesh(FSAVG_Path, lh_pial_path, PLV_Relative_ref);


min_val = min(PLV_Average);
max_val = max(PLV_Average);
PD22N009_PLV_Average_Normalized = (PLV_Average - min_val) / (max_val - min_val);
PD22N009_PLV_Average_Normalized(PD22N009_PLV_Average_Normalized >= 1) = 0.999;

load_ACPC_FR_mesh('SubjectPD22N009_elec_acpc_fr.mat', 'PD22N009\PD22N009\freesurfer\surf\lh.pial', PD22N009_PLV_Average_Normalized);

central_channel = 34;

best_neighbor = getBestNeighborByPLV(PD22N009_PLV_VectorMatrix, central_channel);

comparison_channels = [best_neighbor, 12,14,33,35,54,55,56];

compareRosePlots(PD22N009_Phase_Data, central_channel, comparison_channels);



%PD22N008%%%%%%%%%%%%%%%%%%%%%%

% PD22N008_Phase_Data = obtainPhase(PD22N008_data, 'beta', 1, 10000); %Get the phase data from every channel
% PD22N008_Channel_Names = getChannelNames(PD22N008_data);
% 
% PD22N008_PLV_VectorMatrix = compute_PLV_VectorMatrix(PD22N008_Phase_Data); %Get the PLV vector matrix for every channel (64x64)
% PD22N008_PLV_Average = getAveragePLV(PD22N008_PLV_VectorMatrix); %Get the Average PLV for every channel (64x1)
% ref_channel = 34;
% PD22N008_PLV_Relative_ref = getPLV_wrt_Ref(PD22N008_PLV_VectorMatrix, ref_channel); %Get PLV relative to ref channel (64x1)
% 
% plot_PLV_Relative(PD22N008_Phase_Data,PD22N008_PLV_Relative_ref, ref_channel, PD22N008_Channel_Names);  %Plot a 21x3 grid - color corresponds to PLV wrt ref channel, # corresponds to avg. phase difference wrt ref channel
% 
% load_ACPC_FR_mesh('SubjectPD22N009_elec_acpc_fr.mat', 'PD22N009\PD22N008\freesurfer\surf\lh.pial', PD22N009_PLV_Relative_ref);
% 
% load_FSAVG_FRS_mesh('SubjectPD22N009_elec_fsavg_frs.mat', 'PD22N009\PD22N008\freesurfer\surf\lh.pial', PD22N009_PLV_Relative_ref);
% 
% min_val = min(PD22N008_PLV_Average);
% max_val = max(PD22N008_PLV_Average);
% PD22N008_PLV_Average_Normalized = (PD22N008_PLV_Average - min_val) / (max_val - min_val);
% PD22N008_PLV_Average_Normalized(PD22N008_PLV_Average_Normalized >= 1) = 0.999;
% 
% load_ACPC_FR_mesh('SubjectPD22N008_elec_acpc_fr.mat', 'lh.pial.T1', PD22N008_PLV_Average_Normalized);








function [lag_s, NO_on_BR, TT_overlap] = ...
    align_NO_to_BR(BR_signal, NO_signal, fs, br_win, no_win, ...
                   BR_timetable, NO_timetable, NO_resampled_matrix, ...
                   sigVars, NO_targetDist, BR_to_NO_ch, NO_to_BR_ch, doPlot)
%ALIGN_NO_TO_BR Cross-correlate two independent windows and rebuild NO timeline
%
%   [lag_s, NO_on_BR, TT_overlap] = align_NO_to_BR(...)
%
%   Inputs:
%       BR_signal           - BR signal vector (double)
%       NO_signal           - NO signal vector (double, resampled to BR fs)
%       fs                  - sampling rate (Hz)
%       br_win, no_win      - [start end] windows in seconds
%       BR_timetable        - timetable for BR (with Time + signals)
%       NO_timetable        - original NO timetable (before resample)
%       NO_resampled_matrix - all NO channels resampled to BR fs
%       sigVars             - names of NO channels
%       NO_targetDist       - target distance variable
%       BR_to_NO_ch         - string, BR channel of interest
%       NO_to_BR_ch         - string, NO channel of interest
%       doPlot              - true/false, plot output
%
%   Outputs:
%       lag_s       - estimated lag (s) from xcorr
%       NO_on_BR    - aligned NO timetable (on BR clock)
%       TT_overlap  - struct with overlap indices and times

    % --- Convert windows to sample indices ---
    iBR_start = round(br_win(1)*fs) + 1;
    iBR_end   = round(br_win(2)*fs);
    iNO_start = round(no_win(1)*fs) + 1;
    iNO_end   = round(no_win(2)*fs);

    % --- Extract segments ---
    BR_seg = BR_signal(iBR_start:iBR_end);
    NO_seg = NO_signal(iNO_start:iNO_end);

    % --- Cross correlation ---
    [c, lags] = xcorr(NO_seg - mean(NO_seg), BR_seg - mean(BR_seg), 'coeff');
    [~, iMax] = max(c);
    lag = lags(iMax);
    lag_s = lag / fs;

    % --- Overlap calculation (index alignment) ---
    if lag >= 0
        i0_BR = 1;          i0_NO = lag + 1;
    else
        i0_BR = -lag + 1;   i0_NO = 1;
    end
    L = min(numel(BR_seg)-i0_BR+1, numel(NO_seg)-i0_NO+1);

    % --- Build NO timeline on BR clock ---
    nNO = size(NO_resampled_matrix,1);
    t0_NO = NO_timetable.Time(1);  % absolute first NO timestamp
    t_NO_resamp = t0_NO + seconds((0:nNO-1)'/fs);

    % Anchor point alignment
    iNO_anchor = iNO_start + lag;
    iNO_anchor = max(1, min(nNO, iNO_anchor));
    t_BR_full = BR_timetable.Time;
    t_BR_anchor_abs = t_BR_full(iBR_start);
    t_NO_anchor_abs = t_NO_resamp(iNO_anchor);
    offset = t_BR_anchor_abs - t_NO_anchor_abs;
    t_NO_on_BR = t_NO_resamp + offset;

    % Rebuild aligned NO timetable
    NO_on_BR = [ ...
        array2timetable(NO_resampled_matrix, ...
            'RowTimes', t_NO_on_BR, ...
            'VariableNames', sigVars), ...
        timetable(NO_targetDist, ...
            'RowTimes', t_NO_on_BR, ...
            'VariableNames', {'targetDistance'}) ...
    ];

    % --- Package overlap info ---
    TT_overlap.br_idx = [iBR_start + (i0_BR-1), iBR_start + (i0_BR-1) + L - 1];
    TT_overlap.no_idx = [iNO_start + (i0_NO-1), iNO_start + (i0_NO-1) + L - 1];
    TT_overlap.t_BR   = t_BR_full(TT_overlap.br_idx);
    TT_overlap.t_NO   = t_NO_on_BR(TT_overlap.no_idx);

    % --- Plot if requested ---
    if nargin > 10 && doPlot
        BR_full_sig = BR_timetable.(BR_to_NO_ch);
        NO_full_sig = NO_on_BR.(NO_to_BR_ch);

        figure; hold on; grid on
        plot(t_BR_full, BR_full_sig, 'b');
        plot(NO_on_BR.Time, NO_full_sig, 'r');
        ylim([-1000 1000]);
        xlabel('Time (BR clock)'); ylabel('Amplitude');
        legend('BR','NO (aligned)','Location','best');
        title(sprintf('Full alignment via xcorr (lag = %.3f s)', lag_s));

        % Highlight overlap windows
        patch([TT_overlap.t_BR(1) TT_overlap.t_BR(2) TT_overlap.t_BR(2) TT_overlap.t_BR(1)], ...
              [-1000 -1000 1000 1000], [0 0 1], 'FaceAlpha', 0.1,'EdgeColor','none');
        patch([TT_overlap.t_NO(1) TT_overlap.t_NO(2) TT_overlap.t_NO(2) TT_overlap.t_NO(1)], ...
              [-1000 -1000 1000 1000], [1 0 0], 'FaceAlpha', 0.1,'EdgeColor','none');
    end
end


function [lag_s, c, lags, t_ov, BR_aligned, NO_aligned] = ...
    xcorr_independent_windows(BR_signal, NO_signal, fs, br_win, no_win, doPlot)
%XCORR_INDEPENDENT_WINDOWS Cross-correlate BR and NO signals over independent windows
%
%   [lag_s, c, lags, t_ov, BR_aligned, NO_aligned] = ...
%       xcorr_independent_windows(BR_signal, NO_signal, fs, br_win, no_win, doPlot)
%
%   Inputs:
%       BR_signal - BR signal vector (double)
%       NO_signal - NO signal vector (double)
%       fs        - sampling rate (Hz)
%       br_win    - [start_sec end_sec] window for BR
%       no_win    - [start_sec end_sec] window for NO
%       doPlot    - true/false, whether to plot results
%
%   Outputs:
%       lag_s       - estimated lag (seconds) to align NO to BR
%       c, lags     - xcorr coefficients and lag indices
%       t_ov        - time vector (s) for overlapped aligned segments
%       BR_aligned  - aligned BR segment
%       NO_aligned  - aligned NO segment
%

    % --- Convert windows to sample indices ---
    iBR_start = round(br_win(1) * fs) + 1;
    iBR_end   = round(br_win(2) * fs);
    iNO_start = round(no_win(1) * fs) + 1;
    iNO_end   = round(no_win(2) * fs);

    % --- Extract segments ---
    BR_seg = BR_signal(iBR_start:iBR_end);
    NO_seg = NO_signal(iNO_start:iNO_end);

    % --- Cross correlation ---
    [c, lags] = xcorr(NO_seg - mean(NO_seg), BR_seg - mean(BR_seg), 'coeff');
    [~, iMax] = max(c);
    lag = lags(iMax);
    lag_s = lag / fs;

    % --- Compute overlap after alignment ---
    if lag >= 0
        i0_BR = 1;          i0_NO = lag + 1;
    else
        i0_BR = -lag + 1;   i0_NO = 1;
    end
    L = min(numel(BR_seg) - i0_BR + 1, numel(NO_seg) - i0_NO + 1);

    t_ov = (0:L-1) / fs;
    BR_aligned = BR_seg(i0_BR:i0_BR+L-1);
    NO_aligned = NO_seg(i0_NO:i0_NO+L-1);

    % --- Plot if requested ---
    if nargin > 5 && doPlot
        figure;
        subplot(2,1,1);
        plot(lags/fs, c), grid on
        xlabel('Lag (s)'); ylabel('r'); title('Cross-correlation');

        subplot(2,1,2);
        plot(t_ov, BR_aligned, 'b'); hold on
        plot(t_ov, NO_aligned, 'r');
        grid on
        title(sprintf('Aligned overlap (lag = %.3f s)', lag_s));
        legend('BR','NO');
        xlabel('Time (s)');
    end
end




function TT_window = extractNthLongestWindow(TT_master, sigName, fs, win, trimFrontSec, trimBackSec, nSelect)
%EXTRACTNTHLONGESTWINDOW Keep the n-th longest contiguous window within threshold,
%   with optional separate trimming of front and back.
%
%   TT_window = extractNthLongestWindow(TT_master, sigName, fs, win, trimFrontSec, trimBackSec, nSelect)
%
%   Inputs:
%     TT_master    - timetable (full data)
%     sigName      - name of signal variable to threshold (string/char)
%     fs           - sampling rate (Hz)
%     win          - [lower upper] threshold window
%     trimFrontSec - seconds to cut from the beginning of the window (default = 0)
%     trimBackSec  - seconds to cut from the end of the window (default = 0)
%     nSelect      - which longest window to select (1 = longest, 2 = 2nd longest, etc.)
%
%   Output:
%     TT_window - timetable (subset of TT_master)

    if nargin < 5, trimFrontSec = 0; end
    if nargin < 6, trimBackSec  = 0; end
    if nargin < 7, nSelect = 1; end

    trimFrontSamp = round(trimFrontSec * fs);
    trimBackSamp  = round(trimBackSec  * fs);

    % --- Step 1: Mask valid samples ---
    sig = TT_master.(sigName);
    sig = sig(:);
    valid = (sig >= win(1)) & (sig <= win(2)) & (sig ~= 0) & isfinite(sig);

    % --- Step 2: Find contiguous runs ---
    d = diff([false; valid; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    % --- Step 3: Rank runs by duration ---
    if isempty(starts)
        TT_window = TT_master([],:);  % empty timetable
        return;
    end
    dur = (ends - starts + 1) / fs;
    [durSorted, order] = sort(dur, 'descend');

    if nSelect > numel(order)
        warning('Requested nSelect=%d but only %d valid segments found', nSelect, numel(order));
        TT_window = TT_master([],:);
        return;
    end

    idxChoice = order(nSelect);
    S = starts(idxChoice);
    E = ends(idxChoice);

    % --- Step 4: Apply trimming ---
    S = S + trimFrontSamp;
    E = E - trimBackSamp;

    if S >= E
        warning('Trim too aggressive: window collapsed');
        TT_window = TT_master([],:);
        return;
    end

    % --- Step 5: Slice master timetable ---
    TT_window = TT_master(S:E, :);

    % --- Optional: Plot ---
    figure;
    t = TT_master.Properties.RowTimes;
    plot(t, sig, 'k'); hold on; grid on;
    yline(win(1), '--r','Lower'); yline(win(2), '--r','Upper');
    yl = ylim;
    patch([t(S) t(E) t(E) t(S)], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.85 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');
    uistack(findobj(gca,'Type','line'),'top');
    xlabel('Time'); ylabel(sigName);
    title(sprintf('%d-th longest segment = %.2f s (trim = %.2f s front, %.2f s back)', ...
        nSelect, durSorted(nSelect), trimFrontSec, trimBackSec));
end






function load_ACPC_FR_mesh(acpc_path, pial_lh_path, data)

    % === Step 1: Load Electrode Positions and Mesh ===
    elec_acpc_fr = load(acpc_path).elec_acpc_fr;
    pial_lh = ft_read_headshape(pial_lh_path);
    pial_lh.coordsys = 'acpc';
    

    pial_lh = ft_convert_units(pial_lh, 'mm');
    elec_acpc_fr = ft_convert_units(elec_acpc_fr, 'mm');

    % === Step 2: Identify reference channels ===
    ref_idx = isnan(data);          % Logical array: true if it's reference
    data_clean = data;              % Copy data
    data_clean(ref_idx) = 0;         % Set NaNs to 0 for interpolation

    % === Step 3: Build source structure ===
    source = struct();
    source.pos = elec_acpc_fr.chanpos;
    source.pow = data_clean;
    source.dimord = 'pos';
    source.coordsys = 'acpc';

    % === Step 4: Interpolate electrode data ===
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'sphere_weighteddistance';
    cfg.sphereradius = 3;
    interp_source = ft_sourceinterpolate(cfg, source, pial_lh);

    % === Step 5: Plot the Interpolated Data ===
    cfg = [];
    cfg.funparameter = 'pow';
    cfg.funcolorlim = [0 1];            % normal 0-1 scaling
    cfg.method = 'surface';
    cfg.interpmethod = 'nearest';
    cfg.sphereradius = 3;
    cfg.camlight = 'no';

    ft_sourceplot(cfg, interp_source, pial_lh);

    view([-55 10]);
    material dull;
    lighting gouraud;
    camlight;

    % === Step 6: Normal color map (pure parula) ===
    colormap(parula);   % NO special red added

    % === Step 7: Plot electrodes ===
    % Plot normal electrodes
    elec_plot_colors = data_clean;   % Use interpolated data (0–1)
    ft_plot_sens(elec_acpc_fr, 'elecsize', 35, 'facecolor', elec_plot_colors, 'edgecolor', 'black');

    % === Step 8: Overlay reference electrodes manually ===
    if any(ref_idx)
        ref_positions = elec_acpc_fr.chanpos(ref_idx, :);
        hold on;
        scatter3(ref_positions(:,1), ref_positions(:,2), ref_positions(:,3), ...
            100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold off;
    end

end



function load_FSAVG_FRS_mesh(fsavg_path, pial_lh_path, data)

    % === Step 1: Load Electrode Positions and Mesh ===
    elec_fsavg = load(fsavg_path).elec_fsavg_frs;
    pial_lh = ft_read_headshape(pial_lh_path);
    pial_lh.coordsys = 'fsaverage';

    pial_lh = ft_convert_units(pial_lh, 'mm');
    elec_fsavg = ft_convert_units(elec_fsavg, 'mm');

    % === Step 2: Prepare Electrode Data ===
    ref_idx = isnan(data);           % Find reference channels
    elec_values = data;

    % Replace NaNs with a very large dummy value
    dummy_value = 2;  % well above PLV 1
    elec_values(ref_idx) = dummy_value;  % mark reference channels

    % === Step 3: Build source structure ===
    source = struct();
    source.pos = elec_fsavg.chanpos;
    source.inside = true(size(source.pos,1),1);
    source.pow = elec_values;
    source.dimord = 'pos';
    source.coordsys = 'fsaverage';

    % === Step 4: Interpolate electrode data ===
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'sphere_weighteddistance';
    cfg.sphereradius = 8;
    interp_source = ft_sourceinterpolate(cfg, source, pial_lh);

    % === Step 5: After interpolation, restore reference channel to a clean value ===
    % Any interpolated value >1.5 is from dummy (reference channels)
    ref_mask = interp_source.pow > 1.5;
    interp_source.pow(ref_mask) = 1.001;  % set a standard "red" value

    % === Step 6: Plot brain surface ===
    cfg = [];
    cfg.funparameter = 'pow';
    cfg.funcolorlim = [0 1];  % PLVs go between 0 and 1
    cfg.method = 'surface';
    cfg.interpmethod = 'nearest';
    cfg.sphereradius = 3;
    cfg.camlight = 'no';

    ft_sourceplot(cfg, interp_source, pial_lh);

    view([-55 10]);
    material dull;
    lighting gouraud;
    camlight;

    % === Step 7: Build custom colormap ===
    base_cmap = parula(256);
    red_color = [1 0 0];  % pure red
    cmap_with_red = [base_cmap; red_color];
    colormap(gca, cmap_with_red);

    elec_plot_colors = data;
    elec_plot_colors(ref_idx) = 1.001;  % manually set reference channels to red
    ft_plot_sens(elec_fsavg, 'elecsize', 35, 'facecolor', elec_plot_colors, 'edgecolor', 'black');

end


function chNames = getChannelNames(raw_data)
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
        channelName, channelIndex, channelIndexStim, channelNames] = ...
        getRecordedData_NS(raw_data,1);

    chNames = channelNames;
end

function [phaseData, channelNames] = obtainPhase_tt(TT, SamplingFreq, freq_range)

    % Channel names
    channelNames = TT.Properties.VariableNames;
    numChannels = numel(channelNames);

    % Allocate output
    phaseData = cell(1, numChannels);

    % Loop over channels
    for idx = 1:numChannels
        ch_raw = TT.(channelNames{idx})(:).';

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
end


function phaseData = obtainPhase(raw_data, freq_range, start_idx, stop_idx)
   
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
        channelName, channelIndex, channelIndexStim, channelNames] = ...
        getRecordedData_NS(raw_data,1);

    numChannels = size(dataAllChannels, 1) - 1;

    phaseData = cell(1, numChannels);
    for idx = 1:numChannels
        ch_raw = dataAllChannels(idx, start_idx:stop_idx);
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

end


function PLV_Relative_Vector = getPLV_wrt_Ref(PLV_VectorMatrix, ref_channel)
    % Extracts PLV of all channels with respect to the reference channel
    % Inputs:
    %   - PLV_VectorMatrix: NxN PLV matrix
    %   - ref_channel: reference channel index (scalar)
    % Output:
    %   - PLV_Relative_Vector: Nx1 vector of PLV values relative to ref_channel

    PLV_Relative_Vector = abs(PLV_VectorMatrix(:, ref_channel)); % take absolute value in case of complex numbers
    PLV_Relative_Vector(ref_channel) = NaN; % ignore self-PLV (you could also set to 0 if you want)
end



function PLV_VectorMatrix = compute_PLV_VectorMatrix(phase_data)
    
    numChannels = length(phase_data);
    PLV_vector_matrix = zeros(numChannels, numChannels);
    for ch1 = 1:numChannels
        for ch2 = ch1+1:numChannels
            phase_diff = phase_data{ch1} - phase_data{ch2};
            PLV_vector_matrix(ch1, ch2) = mean(exp(1j * phase_diff));
            PLV_vector_matrix(ch2, ch1) = conj(PLV_vector_matrix(ch1, ch2));
        end
    end

    PLV_VectorMatrix = PLV_vector_matrix;
end

function visualize_PLV_Matrix(PLV_VectorMatrix, channelNames)
%VISUALIZE_PLV_MATRIX Plot PLV matrix as a heatmap
%
%   visualize_PLV_Matrix(PLV_VectorMatrix, channelNames)
%
%   Inputs:
%     PLV_VectorMatrix - NxN complex PLV matrix (from compute_PLV_VectorMatrix)
%     channelNames     - optional cell array of channel labels (1xN)
%
%   Note: Takes abs() of PLV values for visualization.

    % Take absolute value (strength of coupling)
    PLV_abs = abs(PLV_VectorMatrix);

    % Plot
    figure;
    imagesc(PLV_abs);
    axis square;
    cmap = parula(256);   % <-- same colormap as plot_PLV_Relative
    colormap(cmap);
    colorbar;
    caxis([0 1]);  % PLV values are between 0 and 1

    % Label axes if channelNames provided
    numCh = size(PLV_abs,1);
    if nargin > 1 && length(channelNames) == numCh
        xticks(1:numCh); yticks(1:numCh);
        xticklabels(channelNames); yticklabels(channelNames);
        xtickangle(90);
    end

    xlabel('Channel'); ylabel('Channel');
    title('Phase Locking Value (PLV) Matrix');
end



function PLV_Average = getAveragePLV(vector_matrix)
    PLV_Average = mean(abs(vector_matrix), 2);
end


function plot_PLV_Relative(phase_data, ref_phase, PLV_relative_vector, ref_channel, channelNames)
    numChannels = length(phase_data);
    avg_phase_diff = zeros(numChannels, 1);

    ref_name = channelNames{ref_channel};

    for ch = 1:numChannels
        phase_diff = angle(exp(1j * (phase_data{ch} - ref_phase)));
        avg_phase_diff(ch) = mean(phase_diff);
    end
    
    % Best neighbor determination
    PLV_copy = PLV_relative_vector;
    PLV_copy(ref_channel) = -1;
    PLV_copy(63:end) = -1;
    [max_plv, best_neighbor_ch] = max(PLV_copy);
    best_neighbor_phase_diff = avg_phase_diff(best_neighbor_ch);
    
    fprintf('Best neighbor of channel %d is channel %d\n', ref_channel, best_neighbor_ch);
    fprintf('Highest PLV: %.3f\n', max_plv);
    fprintf('Phase difference with channel %d: %.4f radians\n', ...
            ref_channel, best_neighbor_phase_diff);
    
    % Layout
    cmap = parula(256);
    grid_width = 3; grid_height = 21;
    num_channels = grid_width * grid_height;
    electrode_positions = cell(num_channels,1);
    channel_num = 1;
    for x = 0:grid_width-1
        for y = 0:grid_height-1
            electrode_positions{channel_num} = [x,y];
            channel_num = channel_num+1;
        end
    end

    % Normalized values
    plv_norm = PLV_relative_vector;
    plv_norm(isnan(plv_norm)) = 0;
    plv_norm = plv_norm / max_plv;

    % ---- Subplots ----
    figure;
    for subplot_idx = 1:2
        subplot(1,2,subplot_idx); hold on;

        for i = 1:numChannels
            x_i = electrode_positions{i}(1);
            y_i = electrode_positions{i}(2);

            if subplot_idx == 1
                % Raw PLV scaling
                color_val = PLV_relative_vector(i);
                title_str = sprintf('PLV Map wrt %s (Raw)', ref_name);
                label_val = PLV_relative_vector(i);
            else
                % Normalized PLV scaling
                color_val = plv_norm(i);
                title_str = sprintf('PLV Map wrt %s (Normalized)', ref_name);
                label_val = plv_norm(i);
            end

            % Map to colormap index
            color_idx = round(color_val * 255) + 1;
            color_idx = min(max(color_idx,1),256);

            if i == ref_channel
                rectangle('Position',[x_i,y_i,1,1], ...
                          'FaceColor','r','EdgeColor','k');
            else
                rectangle('Position',[x_i,y_i,1,1], ...
                          'FaceColor',cmap(color_idx,:),'EdgeColor','k');
            end

            % Text
            text(x_i+0.5,y_i+0.65,string(channelNames{i}), ...
                 'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
            text(x_i+0.5,y_i+0.35,sprintf('%.3f',label_val), ...
                 'HorizontalAlignment','center','FontSize',6);
        end

        % Highlight best neighbor
        x_best = electrode_positions{best_neighbor_ch}(1);
        y_best = electrode_positions{best_neighbor_ch}(2);
        rectangle('Position',[x_best,y_best,1,1], ...
                  'EdgeColor',[0 0.8 0],'LineWidth',2);

        axis equal; xlim([0,grid_width]); ylim([0,grid_height]);
        set(gca,'XTick',[],'YTick',[],'YDir','normal');
        title(title_str, 'Interpreter','none');
        colormap(cmap); colorbar;
        cbar = colorbar;
        ylabel(cbar, sprintf('PLV with respect to %s', 'Interpreter','none'));
    end
end





function compareRosePlots(phaseData, ref_channel, comparison_channels)
% Compare reference channel to others and plot rose histograms with 95% CI and ±1 STD shading

    figure;
    numComparisons = numel(comparison_channels);

    % --- Decide layout ---
    if numComparisons <= 4
        nRows = 1;
        nCols = numComparisons;
    else
        nRows = 2;
        nCols = ceil(numComparisons / 2);
    end

    for idx = 1:numComparisons
        comp_ch = comparison_channels(idx);

        ref_phase = angle(exp(1j * phaseData{ref_channel}));
        comp_phase = angle(exp(1j * phaseData{comp_ch}));

        phase_diff = angle(exp(1j * (comp_phase - ref_phase)));
    
        [mu, ul, ll] = circ_mean(phase_diff, [], 2);   % Mean, upper, lower 95% CI
        [~, sigma] = circ_std(phase_diff, [], [], 2);  % Circular standard deviation

        % --- Plot ---
        subplot(nRows, nCols, idx);

        % Rose histogram
        h = polarhistogram(phase_diff, 75, 'Normalization', 'probability');
        max_r = max(h.Values); 
        rlim([0 max_r * 1.1]);
        hold on;

        % --- Shaded ±1 STD Cone ---
        theta_fill = linspace(mu - sigma, mu + sigma, 100);
        r_fill = ones(size(theta_fill)) * max_r * 0.8;
        polarplot(theta_fill, r_fill, 'Color', [0.6 0.8 1], 'LineWidth', 6);

        % --- Mean Direction ---
        polarplot([mu mu], [0 max_r], 'r-', 'LineWidth', 2);

        % --- 95% Confidence Interval Lines ---
        polarplot([ul ul], [0 max_r * 0.8], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        polarplot([ll ll], [0 max_r * 0.8], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

        % --- Title ---
        title({
            sprintf('Reference %d vs Neighbor %d', ref_channel, comp_ch), ...
            sprintf('\\mu = %.3f rad', mu), ...
            sprintf('\\sigma = %.3f rad', sigma)
            }, 'FontSize', 10);

    end
end


function plot_PLV_Relative_Brainnetome(phase_data, ref_phase, PLV_relative_vector, ...
                                       ref_channel, channelNames, meta)

    numChannels = length(phase_data);
    avg_phase_diff = zeros(numChannels, 1);

    ref_name = channelNames{ref_channel};

    for ch = 1:numChannels
        phase_diff = angle(exp(1j * (phase_data{ch} - ref_phase)));
        avg_phase_diff(ch) = mean(phase_diff);
    end
    
    % Best neighbor determination (kept same as your reference)
    PLV_copy = PLV_relative_vector;
    PLV_copy(ref_channel) = -1;
    PLV_copy(63:end) = -1;
    [max_plv, best_neighbor_ch] = max(PLV_copy);
    best_neighbor_phase_diff = avg_phase_diff(best_neighbor_ch);
    
    fprintf('Best neighbor of channel %d is channel %d\n', ref_channel, best_neighbor_ch);
    fprintf('Highest PLV: %.3f\n', max_plv);
    fprintf('Phase difference with channel %d: %.4f radians\n', ...
            ref_channel, best_neighbor_phase_diff);
    
    % --- Layout (same as yours) ---
    cmap = parula(256);
    grid_width = 3; grid_height = 21;
    num_channels = grid_width * grid_height;
    electrode_positions = cell(num_channels,1);
    channel_num = 1;
    for x = 0:grid_width-1
        for y = 0:grid_height-1
            electrode_positions{channel_num} = [x,y];
            channel_num = channel_num+1;
        end
    end

    % --- Normalized values ---
    plv_norm = PLV_relative_vector;
    plv_norm(isnan(plv_norm)) = 0;
    plv_norm = plv_norm / max_plv;

    % --- Brainnetome mapping ---
    cleanLabels = regexprep(string(meta.Brainnetome), '\s*\(\d+%?\)\s*$', ''); % strip "(11%)"
    keyElectrodes = lower(strtrim(string(meta.Electrode)));
    brain_map = containers.Map(keyElectrodes, cleanLabels);

    regions = unique(cleanLabels);
    regionColors = lines(numel(regions));
    regionColor = containers.Map(regions, num2cell(regionColors,2));
    legendHandles = gobjects(numel(regions),1);

    % ---- Subplots ----
    figure;
    for subplot_idx = 1:2
        subplot(1,2,subplot_idx); hold on;

        for i = 1:numChannels
            x_i = electrode_positions{i}(1);
            y_i = electrode_positions{i}(2);

            if subplot_idx == 1
                % Raw PLV scaling
                color_val = PLV_relative_vector(i);
                title_str = sprintf('PLV Map wrt %s (Raw)', ref_name);
                label_val = PLV_relative_vector(i);
            else
                % Normalized PLV scaling
                color_val = plv_norm(i);
                title_str = sprintf('PLV Map wrt %s (Normalized)', ref_name);
                label_val = plv_norm(i);
            end

            % Fill color
            color_idx = round(color_val * 255) + 1;
            color_idx = min(max(color_idx,1),256);
            faceCol = cmap(color_idx,:);

            % Border color
            chKey = lower(strtrim(channelNames{i}));
            if isKey(brain_map, chKey)
                region = brain_map(chKey);
                edgeCol = regionColor(region);
            else
                region = "Unlabeled";
                edgeCol = 'k';
            end

            % Draw rectangle
            if i == ref_channel
                rectangle('Position',[x_i,y_i,1,1], ...
                          'FaceColor','r','EdgeColor','k','LineWidth',2);
            else
                rectangle('Position',[x_i,y_i,1,1], ...
                          'FaceColor',faceCol,'EdgeColor',edgeCol,'LineWidth',2);
            end

            % Add proxy for legend
            if subplot_idx==1 && region ~= "Unlabeled"
                idx = find(strcmp(regions, region));
                if isempty(legendHandles(idx)) || ~isgraphics(legendHandles(idx))
                    legendHandles(idx) = plot(nan,nan,'s', ...
                        'MarkerSize',10,'MarkerEdgeColor',edgeCol, ...
                        'MarkerFaceColor','w','LineWidth',2);
                end
            end

            % Text
            text(x_i+0.5,y_i+0.65,string(channelNames{i}), ...
                 'HorizontalAlignment','center','FontSize',7,'FontWeight','bold');
            text(x_i+0.5,y_i+0.35,sprintf('%.3f',label_val), ...
                 'HorizontalAlignment','center','FontSize',6);
        end

        % (keep your neighbor outline, or remove if not needed)
        x_best = electrode_positions{best_neighbor_ch}(1);
        y_best = electrode_positions{best_neighbor_ch}(2);
        rectangle('Position',[x_best,y_best,1,1], ...
                  'EdgeColor',[0 0.8 0],'LineWidth',2);

        axis equal; xlim([0,grid_width]); ylim([0,grid_height]);
        set(gca,'XTick',[],'YTick',[],'YDir','normal');
        title(title_str, 'Interpreter','none');
        colormap(cmap); colorbar;
        ylabel(colorbar, sprintf('PLV with respect to %s', ref_name), 'Interpreter','none');

        % Add legend
        if subplot_idx == 1
            validIdx = isgraphics(legendHandles);
            legend(legendHandles(validIdx), regions(validIdx), ...
                   'Location','eastoutside','Interpreter','none');
        end
    end
end
