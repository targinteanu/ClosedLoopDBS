    %% Patient presets
    % =========================
    presets = struct();
    
    presets.PD25N006 = struct( ...
        "Patient_ID",        "PD25N006", ...
        "NO_timetable_name", "SavedTable_Distances_1375Hz1.mat", ...
        "NO_SamplingFreq",   1375, ...
        "NO_to_BR_ch",       "CMacro_LFP_01", ...
        "BR_to_NO_ch",       "AINP3", ...
        "ns_files",          {{'datafile001.ns2'}}, ...
        "use_ns2_only",      true, ...
        "target_window",     [-5 3], ...
        "target_channel",    "CMacro_LFP_01", ...
        "minDurSec",         20, ...
        "br_win",       [1800 2500], ...
        "no_win",       [1700 2400], ...
        "frontCut",    60, ...
        "backCut",     3800, ...
        "electrode_data_name", "electrode_data(Sheet1).csv", ...
        "badChannel",        [] ...
    );
    
    presets.PD25N008 = struct( ...
        "Patient_ID",        "PD25N008", ...
        "NO_timetable_name", "SavedTable_Distances_1375Hz1.mat", ...
        "NO_SamplingFreq",   1375, ...
        "NO_to_BR_ch",       "CMacro_LFP_AP_T", ...
        "BR_to_NO_ch",       "AINP2", ...
        "ns_files",          {{'baseline001.ns2','DepthStim002.ns2','nback004.ns2'}}, ...
        "use_ns2_only",      true, ...
        "target_window",     [-5 3], ...
        "target_channel",    "CMacro_LFP_AP_T", ...
        "minDurSec",         20, ...
        "br_win",       [800 1500], ...
        "no_win",       [400 1100], ...
        "frontCut",    20, ...
        "backCut",     20, ...
        "electrode_data_name", "PD25N008_electrode_data", ...
        "badChannel",        [48] ...
    );
    
    presets.PD24N009 = struct( ...
        "Patient_ID",        "PD24N009", ...
        "NO_timetable_name", "SavedTable_Distances_1375Hz1.mat", ...
        "NO_SamplingFreq",   1375, ...
        "NO_to_BR_ch",       "CMacro_LFP_01", ...
        "BR_to_NO_ch",       "AINP4_ns5", ...              % you can override later for ns5 case
        "ns_files",          {{'datafile001.ns2','datafile001.ns5'}}, ...
        "use_ns2_only",      false, ...
        "target_window",     [-5 3], ...
        "target_channel",    "CMacro_LFP_01", ...
        "minDurSec",         20, ...
        "br_win",       [600 850], ...
        "no_win",       [2900 3150], ...
        "frontCut",    600, ...
        "backCut",     440, ...
        "electrode_data_name", "PD24N009_electrode_data", ...
        "badChannel",        [] ...
    );
    
    %%
    
    % Choose patient here
    % =========================
    Patient_Folder_Path = 'C:\Users\cours\Downloads\AndersonLabDataAndPackages\Patient Data';
    whichPatient = 'PD25N008';
    
    P = presets.(whichPatient);
    
    
    Patient_ID        = char(P.Patient_ID);
    NO_timetable_name = P.NO_timetable_name;
    NO_SamplingFreq   = P.NO_SamplingFreq;
    ns_files          = P.ns_files;
    use_ns2_only      = P.use_ns2_only;
    target_channel    = P.target_channel;
    minDurSec         = P.minDurSec;
    NO_to_BR_ch       = P.NO_to_BR_ch;
    
    
    BR_SamplingFreq = 1000;
    
    
    % ----------- Load NO Data -----------
    NO_timetable_path = fullfile(Patient_Folder_Path, Patient_ID, NO_timetable_name);
    NO_timetable_OG = load(NO_timetable_path); 
    
    % Sort time chronologically
    NO_timetable_unmodified = sortrows(NO_timetable_OG.Tbl, 'Time'); 
    
    % Make sure NO is at its supposed sampling frequency (and filling in 0
    % otherwise)
    NO_timetable = myRetime(NO_timetable_unmodified, NO_SamplingFreq, 0, 2);
    
    NO_targetDist_unmodified = NO_timetable.("targetDistance");
    
    %The NO signal of the channel sent to BlackRock
    NO_signal_unmodified = NO_timetable.(NO_to_BR_ch); % or other channel
    
    % ----------- Load BR Data -----------
    ns_timetables = cell(size(ns_files));
    for k = 1:numel(ns_files)
        ns_path = fullfile(Patient_Folder_Path, Patient_ID, ns_files{k});
        disp(['Loading: ', ns_path])
        ns = openNSx(ns_path);   
        fs_curr = ns.MetaTags.SamplingFreq;
        tt = ns2timetable(ns_path);
    
        if fs_curr > 1000
            disp(['Downsampling ' ns_files{k} ' from ' num2str(fs_curr) ' Hz to 1 kHz...']);
            data_ds = resample(double(tt{:, :}), 1000, fs_curr);
            t_ds = tt.Time(1) + seconds((0:size(data_ds,1)-1)'/1000);
            tt = array2timetable(data_ds, 'RowTimes', t_ds, ...
                         'VariableNames', tt.Properties.VariableNames);
            fs_curr = 1000;
        else
            tt.SourceFile = repmat(string(ns_files{k}), height(tt), 1);
        end
        ns_timetables{k} = tt;
    end
    
    %Make timetable of all the ns files, and sort time chronologically
    combined_tt = vertcat(ns_timetables{:});
    combined_tt = sortrows(combined_tt, 'Time');
    
    
    %% Just for PD24N009
    
    Patient_Folder_Path = 'C:\Users\cours\Downloads\AndersonLabDataAndPackages\Patient Data';
    Patient_ID = 'PD24N009';
    BR_to_NO_ch = 'AINP4_ns5';
    ns2_file = fullfile(Patient_Folder_Path, Patient_ID, 'datafile001.ns2');
    ns5_file = fullfile(Patient_Folder_Path, Patient_ID, 'datafile001.ns5');
    
    T2 = ns2timetable(ns2_file);
    T5 = ns2timetable(ns5_file);
    T2 = sortrows(T2,'Time'); T5 = sortrows(T5,'Time');
    
    % Trim to common overlap
    t0 = max(T2.Time(1), T5.Time(1));
    t1 = min(T2.Time(end), T5.Time(end));
    T2 = T2(T2.Time>=t0 & T2.Time<=t1,:);
    T5 = T5(T5.Time>=t0 & T5.Time<=t1,:);
    
    % Resample ns5 to ns2 timestamps (1kHz)
    T5r = retime(T5, T2.Time, 'linear');
    
    % Rename ns5 vars to avoid name conflicts
    T5r = renamevars(T5r, T5r.Properties.VariableNames, ...
                     append(T5r.Properties.VariableNames,"_ns5"));
    
    % Merge horizontally
    combined_tt = [T2 T5r];
    
    
    
    %% Set up resampled matrices
    
    BR_to_NO_ch       = P.BR_to_NO_ch;
    % Make sure BR is at its supposed sampling frequency (and filling in 0
    % otherwise)
    BR_timetable = myRetime(combined_tt, BR_SamplingFreq, 0, 2);
    
    % ----------- Resample NO data to BR sampling rate -----------
    t_NO = (0:height(NO_targetDist_unmodified)-1)' / NO_SamplingFreq;
    duration = (height(NO_targetDist_unmodified)-1) / NO_SamplingFreq;
    t_BR = (0:round(duration*BR_SamplingFreq))' / BR_SamplingFreq;
    minLen = min(numel(NO_targetDist_unmodified), numel(t_NO));
    NO_targetDist_unmodified = NO_targetDist_unmodified(1:minLen);
    t_NO = t_NO(1:minLen);
    
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
    
    
    BR_signal = double(BR_signal(:));
    NO_signal = double(NO_signal(:));
    BR_signal(BR_signal > 2000)  = 2000;
    BR_signal(BR_signal < -2000) = -2000;
    NO_signal(NO_signal > 2000)  = 2000;
    NO_signal(NO_signal < -2000) = -2000;
    BR_signal = BR_signal - mean(BR_signal);
    NO_signal = NO_signal - mean(NO_signal);
    
    t_BR = (0:length(BR_signal)-1)' / BR_SamplingFreq;
    t_NO = (0:length(NO_signal)-1)' / BR_SamplingFreq;
    t_NO_dist = (0:length(NO_targetDist)-1)' / BR_SamplingFreq;
    
    figure;
    subplot(3,1,1);
    plot(t_BR, BR_signal, 'b'); hold on;
    xlabel('Time (s)'); ylabel('Amplitude');
    title('BR\_signal (clipped + demeaned, pre-alignment)');
    ylim([-2000 2000]);
    segment_lengths = cellfun(@(tt) height(tt), ns_timetables);
    segment_lengths = segment_lengths(:).';           % force row vector
    segment_edges   = [0, cumsum(segment_lengths)] / BR_SamplingFreq;
    
    for k = 1:numel(ns_files)
        if k > 1
            xline(segment_edges(k), '--r');
        end
    
        mid_time = (segment_edges(k) + segment_edges(k+1)) / 2;
    
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
    
    fs        = 1000;
    br_win    = P.br_win;
    no_win    = P.no_win;
    BR_to_NO_ch       = P.BR_to_NO_ch;
    
    % br_win = [1800 2500];   %%PD25N006 - br_win = [1800 2500], no_win = [1700 2400] %%PD25N008 - [800 1500] [400 1100]
    % no_win = [1700 2400];   %%PD25N009 - br_win = [500 1000], no_win = [2800 3300] %%PD24N009 - [600 850] [2900 3150]
    
    if use_ns2_only == false
        minLen = min(numel(NO_targetDist), size(NO_resampled_matrix,1));
        NO_targetDist = NO_targetDist(1:minLen);
        NO_resampled_matrix = NO_resampled_matrix(1:minLen,:);
    end
    
    [lag_s, NO_on_BR, TT_overlap] = align_NO_to_BR( ...
        BR_signal, NO_signal, fs, br_win, no_win, ...
        BR_timetable, NO_timetable, NO_resampled_matrix, ...
        sigVars, NO_targetDist, BR_to_NO_ch, NO_to_BR_ch, true);
    
    fprintf('Lag = %.3f s\n', lag_s);
    disp(TT_overlap)
    
    
    %% Building Full Synchronized Table
    
    fs = 1000;
    t_common = (max(BR_timetable.Time(1), NO_on_BR.Time(1)) : ...
                seconds(1/fs) : ...
                min(BR_timetable.Time(end), NO_on_BR.Time(end)))';
    
    if any(strcmp("SourceFile", BR_timetable.Properties.VariableNames))
        TT_master_num = synchronize(removevars(BR_timetable,"SourceFile"), NO_on_BR, t_common, 'linear');
        % Resample SourceFile separately (nearest)
        TT_source = retime(BR_timetable(:, "SourceFile"), t_common, 'nearest');
        % Merge back in
        TT_master = [TT_master_num TT_source];
    else
        % No SourceFile column → just synchronize everything directly
        TT_master = synchronize(BR_timetable, NO_on_BR, t_common, 'linear');
    end
    
    t = TT_master.Properties.RowTimes;
    figure;
    plot(t, TT_master.(BR_to_NO_ch), 'b'); hold on
    plot(t, TT_master.(NO_to_BR_ch), 'r');
    legend('BR','NO'); grid on
    xlabel('Time'); ylabel('Signal');
    title('NO full span with BR zero-padded');
    ylim([-1000 1000]);
    
    
    
    %% Longest NO segment within the Target Distance Threshold (from TT_master)
    
    fs = 1000;
    sigName = 'targetDistance'; 
    win     = P.target_window;
    frontCut    = P.frontCut;
    backCut     = P.backCut;
    
    % frontCut = 60; %60  %%PD25N006 - frontCut = 900, backCut = 2525 %%PD25N008 - 20,20 %%PD24N009 - 600, 440
    % backCut = 3800; %3800
    TT_longest = extractNthLongestWindow(TT_master, sigName, fs, win, frontCut, backCut, 1, BR_to_NO_ch);
    
    if any(strcmp('SourceFile', TT_longest.Properties.VariableNames))
        TT_longest_clean = removevars(TT_longest, {'SourceFile'});
    else
        TT_longest_clean = TT_longest;
    end
    
    
    %% 2nd Cross Correlation on the Isolated NO Segment
    
    fs = 1000;  
    
    BR_seg = TT_longest.(BR_to_NO_ch);
    NO_seg = TT_longest.(NO_to_BR_ch);
    
    [c,lags] = xcorr(NO_seg - mean(NO_seg), BR_seg - mean(BR_seg), 'coeff');
    [~,imax] = max(c);
    refined_lag = lags(imax);          % samples
    refined_lag_s = refined_lag/fs;    % seconds
    
    fprintf('Refined lag in TT_longest: %d samples (%.6f s)\n', refined_lag, refined_lag_s);
    
    NO_on_BR_refined = NO_on_BR;
    NO_on_BR_refined.Time = NO_on_BR.Time - seconds(refined_lag_s);
    
    t_common = (max(BR_timetable.Time(1), NO_on_BR_refined.Time(1)) : ...
                seconds(1/fs) : ...
                min(BR_timetable.Time(end), NO_on_BR_refined.Time(end)))';
    
    if any(strcmp("SourceFile", BR_timetable.Properties.VariableNames))
        TT_master_num = synchronize(removevars(BR_timetable,"SourceFile"), NO_on_BR_refined, t_common, 'linear');
        TT_source = retime(BR_timetable(:, "SourceFile"), t_common, 'nearest');
        TT_master_refined = [TT_master_num TT_source];
    else
        TT_master_refined = synchronize(BR_timetable, NO_on_BR_refined, t_common, 'linear');
    end
    
    
    TT_longest_refined = extractNthLongestWindow(TT_master_refined, ...
                                              'targetDistance', fs, win, frontCut, backCut,1);
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
    
    %% Plot synced signals against each other
    
    ECoG_ch = 'LS33';  
    
    signals = {"targetDistance", NO_to_BR_ch, BR_to_NO_ch, ECoG_ch};
    
    titles  = { ...
        "Target Distance (thresholds shown)", ...
        sprintf("NO signal: %s", NO_to_BR_ch), ...
        sprintf("Subthalamic Nucleus (STN) Recording: %s", BR_to_NO_ch), ...
        sprintf("ECoG Channel: %s", ECoG_ch) ...
        };
    colors  = {'k','m','r','b'};   
    
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
    
    %% ADF Test (testing for stationarity)
    % Extract signals
    x_ecog = TT_longest_refined.LS33;
    x_lfp  = TT_longest_refined.CMacro_LFP_01;
    
    % Remove NaNs if present
    % x_ecog = x_ecog(~isnan(x_ecog));
    % x_lfp  = x_lfp(~isnan(x_lfp));
    [h_ecog,p_ecog,stat_ecog,cval_ecog] = adftest(x_ecog);
    [h_lfp, p_lfp, stat_lfp, cval_lfp ] = adftest(x_lfp);
    
    fprintf('ECoG LS33: h = %d, p = %.5f\n', h_ecog, p_ecog);
    fprintf('LFP  CMacroLFP01: h = %d, p = %.5f\n', h_lfp, p_lfp);
    fprintf('ECoG LS33: ADF stat = %.4f\n', stat_ecog);
    fprintf('LFP  CMacroLFP01: ADF stat = %.4f\n', stat_lfp);
    
    
    t = TT_longest_refined.Time;
    
    figure;
    subplot(2,1,1)
    plot(t, x_ecog, 'k')
    title('ECoG LS33 (raw)')
    ylabel('Amplitude')
    grid on
    
    subplot(2,1,2)
    plot(t, x_lfp, 'r')
    title('LFP CMacroLFP01 (raw)')
    ylabel('Amplitude')
    xlabel('Time')
    grid on
    
    
    %% Get rid of Source File label
    
    TT_longest_refined = removevars(TT_longest_refined, {'SourceFile'});
    
    
    %% Grab NO time segment on its own vector (not required)
    % --- Infer offset using REFINED alignment (IMPORTANT) ---
    NO_to_BR_offset = NO_on_BR_refined.Time(1) - NO_timetable.Time(1);
    
    % --- Map TT_longest_refined bounds back to NO time ---
    tNO_start = TT_longest_refined.Time(1)  - NO_to_BR_offset;
    tNO_end   = TT_longest_refined.Time(end) - NO_to_BR_offset;
    
    % --- Plot in NO's native time space ---
    figure; hold on;
    plot(NO_timetable.Time, NO_timetable.(NO_to_BR_ch), 'k');
    xline(tNO_start, 'r', 'LineWidth', 2);
    xline(tNO_end,   'r', 'LineWidth', 2);
    
    % --- Overlay NO events ---
    E = NO_timetable.Properties.Events;
    xline(E.Time, '--b');
    
    grid on;
    xlabel('NO time (native)');
    ylabel('NO signal');
    title('TT\_longest\_refined window in NO timetable time');
    legend('NO signal','Window start','Window end','Events');
    
    % --- Print absolute NO times ---
    fprintf('NO absolute time:\n');
    fprintf('  Start: %.3f s\n', seconds(tNO_start));
    fprintf('  End:   %.3f s\n', seconds(tNO_end));
    fprintf('  Dur:   %.3f s\n', seconds(tNO_end - tNO_start));
    
    % --- Get ACTIVE event at window start (not nearest!) ---
    iEvt = find(E.Time <= tNO_start, 1, 'last');
    
    if isempty(iEvt)
        event_label = "<none>";
    else
        event_label = E.EventLabels(iEvt);
    end
    
    fprintf('Event at window start: %s\n', string(event_label));
    
    
    % --- Find all events that occur INSIDE the window ---
    idxEventsInWindow = E.Time > tNO_start & E.Time < tNO_end;
    
    events_crossed = E.EventLabels(idxEventsInWindow);
    
    if isempty(events_crossed)
        fprintf('Window stays within a single event.\n');
    else
        fprintf('Window crosses %d additional event(s):\n', numel(events_crossed));
        disp(string(events_crossed));
    end
    
    
    
    
    
    %% Calculate PLV vector matrices
    
    [Phase_Data, Channel_Names] = obtainPhase_tt(TT_longest_refined, 1000, 'beta');
    disp(Channel_Names)
    disp('Obtained Phase Data')
    PLV_VectorMatrix = compute_PLV_VectorMatrix(Phase_Data); %Get the PLV vector matrix for every channel (64x64)
    PLV_Average = getAveragePLV(PLV_VectorMatrix); %Get the Average PLV for every channel (64x1)
    disp('Calculated PLV Vector Matrix')
    
    
    %% Visualize PLV Vector Matrix
    visualize_PLV_Matrix(PLV_VectorMatrix, Channel_Names)
    
    
    %% Sanity Check Phase
    fs = 1000;          % Hz
    duration = 5;       % seconds
    nSamp = duration * fs;
    
    chA = 33; chB = 70;
    
    start_time = 0.525;                   % seconds
    start_idx  = round(start_time * fs) + 1;   % MATLAB 1-based indexing
    
    t = (0:nSamp-1) / fs;
    
    figure;
    t_shifted = t(start_idx:end) - t(start_idx);
    
    plot(t_shifted, Phase_Data{chA}(start_idx:nSamp), 'b'); hold on;
    plot(t_shifted, Phase_Data{chB}(start_idx:nSamp), 'r');
    
    xlabel('Time (s)');
    ylabel('Phase (rad)');
    
    
    xlabel('Time (s)');
    ylabel('Phase (rad)');
    
    labelA = 'ECoG Recording (Ch 33)';
    labelB = 'STN Recording';
    
    title(sprintf('Phase of %s and %s (250ms window)', labelA, labelB));
    legend(labelA, labelB);
    
    grid on;
    
    
    
    %% Calculate PLV wrt Reference Channel
    
    %NO_to_BR_ch = 'CMacro_LFP_01'; %CMacro_LFP_AP_T AINP3 PD25N006 - CMacro_LFP_01 PD24N009 - CMacro_LFP_01
    NO_to_BR_ch       = P.NO_to_BR_ch;
    
    ref_idx = find(strcmp(Channel_Names, NO_to_BR_ch));
    fprintf('"%s" is at channel index %d\n', Channel_Names{ref_idx}, ref_idx);
    disp(ref_idx)
    ref_channel = ref_idx;
    PLV_Relative_ref = getPLV_wrt_Ref(PLV_VectorMatrix, ref_channel); %Get PLV relative to ref channel (64x1)
    disp('Calculated PLV wrt Reference')
    
    
    %% PLV Summary Statistics
    
    nSur = 100;
    
    channels = [19 20 21 40 41 42 61 62 63];%setdiff(1:63, 48);   % exclude bad channel 48
    
    
    % ---------- Reference channel ----------
    NO_to_BR_ch = 'CMacro_LFP_01';
    ref_channel = find(strcmp(Channel_Names, NO_to_BR_ch));
    fprintf('"%s" is at channel index %d\n', Channel_Names{ref_channel}, ref_channel);
    
    % ---------- Real PLV ----------
    PLV = getPLV_wrt_Ref(PLV_VectorMatrix, ref_channel);
    PLV_raw = PLV(channels);
    
    % ---------- Normalized PLV (relative to strongest channel) ----------
    PLV_norm = PLV_raw / max(PLV_raw);
    
    % ---------- Brainnetome labels ----------
    % ---------- Brainnetome labels (MATCH plot_PLV_Relative_Brainnetome) ----------
    cleanLabels = regexprep(string(meta.Brainnetome), '\s*\(\d+%?\)\s*$', '');
    
    % Extract electrode numbers from meta (e.g. "1", "L01", "LS01")
    metaNums = regexp(string(meta.Electrode), '\d+', 'match', 'once');
    metaNums = string(pad(metaNums, 2, 'left', '0'));   % zero-pad: "1" → "01"
    
    % Build electrode-number → region map
    brain_map = containers.Map(metaNums, cleanLabels);
    
    % Assign region to each channel using Channel_Names
    region_labels = strings(numel(channels),1);
    for i = 1:numel(channels)
        chNum = regexp(string(Channel_Names{channels(i)}), '\d+', 'match', 'once');
        if isKey(brain_map, chNum)
            region_labels(i) = brain_map(chNum);
        else
            region_labels(i) = "Unlabeled";
        end
    end
    
    regions = unique(region_labels);
    
    
    % ---------- Highest-AVG raw PLV region ----------
    meanPLV = zeros(numel(regions),1);
    for r = 1:numel(regions)
        meanPLV(r) = mean(PLV_raw(strcmp(region_labels, regions{r})), 'omitnan');
    end
    
    [~,ridx] = max(meanPLV);
    topRegion = regions{ridx};
    
    idx_top   = strcmp(region_labels, topRegion);
    idx_other = ~idx_top;
    
    % ---------- Surrogate PLV (pipeline-consistent) ----------
    PLV_surrogate = nan(numel(channels), nSur);
    
    ref_phase = Phase_Data{ref_channel};
    T = numel(ref_phase);
    
    for s = 1:nSur
        shift = randi(T);
    
        Phase_Data_surr = Phase_Data;
        Phase_Data_surr{ref_channel} = circshift(ref_phase, shift);
    
        PLV_VM_surr = compute_PLV_VectorMatrix(Phase_Data_surr);
        PLV_tmp = getPLV_wrt_Ref(PLV_VM_surr, ref_channel);
    
        PLV_surrogate(:,s) = PLV_tmp(channels);
    
        % ---- Progress print ----
        if mod(s,10)==0 || s==1 || s==nSur
            fprintf('  Surrogate %3d / %3d completed\n', s, nSur);
        end
    end
    
    disp('Finished surrogate PLV generation.');
    
    % ---------- Surrogate means ----------
    mu_surr = mean(PLV_surrogate,2);
    
    % ---------- Highest PLV channel ----------
    [~,ch_max] = max(PLV_raw);
    
    % ---------- Metrics ----------
    Top_mean_raw   = mean(PLV_raw(idx_top),  'omitnan');
    Top_mean_norm  = mean(PLV_norm(idx_top), 'omitnan');
    
    Other_mean_raw  = mean(PLV_raw(idx_other),  'omitnan');
    Other_mean_norm = mean(PLV_norm(idx_other), 'omitnan');
    
    TopCh_surr_raw  = mu_surr(ch_max);
    TopCh_surr_norm = mu_surr(ch_max) / PLV_raw(ch_max);
    
    Pct_above = 100 * mean(PLV_raw(idx_top) > mu_surr(idx_top));
    
    % ---------- Output ----------
    fprintf('\n===== PLV Summary =====\n');
    fprintf('Highest AVG raw PLV region:      %s\n', topRegion);
    fprintf('Top region avg raw PLV:          %.4f\n', Top_mean_raw);
    fprintf('Top region avg normalized PLV:   %.3f\n', Top_mean_norm);
    fprintf('Other regions avg raw PLV:       %.4f\n', Other_mean_raw);
    fprintf('Other regions avg normalized:    %.3f\n', Other_mean_norm);
    fprintf('Highest PLV channel surrogate:   raw=%.4f | norm=%.3f\n', ...
            TopCh_surr_raw, TopCh_surr_norm);
    fprintf('%% channels in %s > surrogate:    %.1f%%\n', topRegion, Pct_above);
    highestPLV_channel_name = Channel_Names{channels(ch_max)};
    highestPLV_raw  = PLV_raw(ch_max);
    highestPLV_norm = PLV_norm(ch_max);
    fprintf('Highest PLV channel:             %s (raw=%.4f, norm=%.3f)\n', ...
            highestPLV_channel_name, highestPLV_raw, highestPLV_norm);
    
    % ---------- Z-PLV for highest PLV channel ----------
    
    % Surrogate mean and std for that channel
    mu_ch = mean(PLV_surrogate(ch_max, :));
    sd_ch = std( PLV_surrogate(ch_max, :));
    
    % Z-PLV
    Z_PLV = (PLV_raw(ch_max) - mu_ch) / sd_ch;
    
    fprintf('Highest PLV channel Z-PLV:       %.2f (%s)\n', ...
            Z_PLV, Channel_Names{channels(ch_max)});
    
    
    %% Plot PLV wrt Reference Channel with Brainnetome labeling
    
    %electrode_data_name = 'PD24N009_electrode_data'; %%PD25N006 - 'electrode_data(Sheet1).csv' PD25N008 - 'PD25N008_electrode_data' PD24N009 - PD24N009_electrode_data
    electrode_data_name = P.electrode_data_name;
    electrode_data_filepath = fullfile(Patient_Folder_Path, Patient_ID, electrode_data_name);
    disp(electrode_data_filepath)
    meta = readtable(electrode_data_filepath, 'VariableNamingRule','preserve');
    
    badChannel = P.badChannel;
    plot_PLV_Relative_Brainnetome(Phase_Data(1:63), Phase_Data{ref_channel}, PLV_Relative_ref, ref_channel, Channel_Names, meta, badChannel);
    
    
    %%
    
    plot_PLV_Relative(Phase_Data(1:63),Phase_Data{ref_channel},PLV_Relative_ref, ref_channel, Channel_Names);  %Plot a 21x3 grid - color corresponds to PLV wrt ref channel, # corresponds to avg. phase difference wrt ref channel
    
    %% Calculating Ideal MVAR model order using AIC
    maxOrder  = 50;
    fs        = 1000;
    blockSize = 250;
    TT_all    = TT_longest_refined;
    chIdx     = [61 70];
    
    nTotal  = height(TT_all);
    nBlocks = floor(nTotal / blockSize);
    
    stepBlocks = 20;                          % <-- CHANGED
    evalBlocks = 1:stepBlocks:nBlocks;
    nEval      = numel(evalBlocks);
    
    AIC_all = NaN(maxOrder, nEval);           % <-- NEW: store full AIC curves
    bestOrders = NaN(nEval,1);
    
    for k = 1:nEval
        b = evalBlocks(k);
    
        idxStart = (b-1)*blockSize + 1;
        idxEnd   = b*blockSize;
    
        TT = TT_all(idxStart:idxEnd, :);
    
        % --- Extract data matrix ---
        X = TT{:, chIdx}';
        nSamp = size(X,2);
        nCh   = size(X,1);
    
        AIC_vals = NaN(maxOrder,1);
    
        for p = 1:maxOrder
            try
                [~,E] = armorf(X, 1, nSamp, p);
                num_params = nCh^2 * p;
                AIC_vals(p) = log(det(E)) + (2 * num_params) / nSamp;
            catch
                AIC_vals(p) = NaN;
            end
        end
    
        AIC_all(:,k) = AIC_vals;              % <-- STORE CURVE
        [~, bestOrders(k)] = min(AIC_vals);
    
        fprintf('Block %d (%d:%d): best order = %d\n', ...
                b, idxStart, idxEnd, bestOrders(k));
    end
    
    figure; hold on;
    plot(1:maxOrder, AIC_all, 'LineWidth', 1);
    xlabel('Model order p');
    ylabel('AIC');
    title(sprintf('AIC curves (every %dth block)', stepBlocks));
    grid on;
    
    
    for k = 1:nEval
        [~, pmin] = min(AIC_all(:,k));
        plot(pmin, AIC_all(pmin,k), 'ko', 'MarkerFaceColor','k');
    end
    
    
    avgOrder    = mean(bestOrders, 'omitnan');
    medianOrder = median(bestOrders, 'omitnan');
    maxOrderEst = max(bestOrders);
    fprintf('\nModel order stats (every %dth block):\n', stepBlocks);
    fprintf('  Average order: %.2f\n', avgOrder);
    fprintf('  Median order : %.1f\n', medianOrder);
    fprintf('  Max order    : %d\n', maxOrderEst);
    
    
    
    %% Mann-Whitney U Test of GC (Alternative to 2 sample t test)
    
    
    
    %OLD VERSION BEFORE 1/9/26
    segment_length_sec = 4;    % segment length (seconds)
    fs = 1000;
    order = 50; %model order you get from AIC
    freq_band = 'beta';
    chIdx = [33 69];            % [A B] [61 70]
    
    % === Use channel names for labeling ===
    A_name = Channel_Names{chIdx(1)};
    B_name = Channel_Names{chIdx(2)};
    
    TT = TT_longest_refined;
    
    % Pre-filter the full-length signals once (safe length) NEW
    filtOrder = 1024;
    
    A_full = double(TT{:, chIdx(1)});  A_full = A_full(:).';   % force row
    B_full = double(TT{:, chIdx(2)});  B_full = B_full(:).';   % force row
    
    A_beta = Myeegfilt(A_full, fs, 13, 30, 0, filtOrder);
    B_beta = Myeegfilt(B_full, fs, 13, 30, 0, filtOrder);
    % End of new
    
    
    samples_per_segment = segment_length_sec * fs;
    nTotalSamples = height(TT);
    nSegments = floor(nTotalSamples / samples_per_segment);
    
    fprintf('Splitting into %d segments of %.1f s each (%.0f samples)\n', ...
            nSegments, segment_length_sec, samples_per_segment);
    
    GC_all = zeros(2,2,nSegments);
    GC_AB_seg = zeros(nSegments,1);
    GC_BA_seg = zeros(nSegments,1);
    
    for k = 1:nSegments
        idxStart = (k-1)*samples_per_segment + 1;
        idxEnd   = k*samples_per_segment;
        Tseg = TT(idxStart:idxEnd, chIdx);   % only keep the two channels
        Tseg{:,1} = A_beta(idxStart:idxEnd).';
        Tseg{:,2} = B_beta(idxStart:idxEnd).';
        [GC,~,~,~] = compute_static_GC(Tseg, fs, order, [], 'none', false, false);
    
        GC_all(:,:,k) = GC;
        GC_AB_seg(k)  = GC(2,1);
        GC_BA_seg(k)  = GC(1,2);
    
        fprintf('Segment %2d/%2d done | GC(%s→%s)=%.4f | GC(%s→%s)=%.4f\n', ...
            k, nSegments, A_name, B_name, GC_AB_seg(k), B_name, A_name, GC_BA_seg(k));
    end
    
    GC_mean = mean(GC_all, 3, 'omitnan');
    
    % === Statistics ===
    [p_val,~,stats] = signrank(GC_AB_seg, GC_BA_seg);
    zval = stats.zval;
    fprintf('\nWilcoxon signed-rank: z = %.3f, p = %.4g\n', zval, p_val);
    
    [~,p_ttest,~,stats_t] = ttest(GC_AB_seg, GC_BA_seg);
    fprintf('Paired t-test: t(%d) = %.3f, p = %.4g\n', stats_t.df, stats_t.tstat, p_ttest);
    
    diffs = GC_AB_seg - GC_BA_seg;
    qqplot(diffs);
    title(sprintf('QQ Plot of paired differences: %s↔%s', A_name, B_name));
    %% Mean and SD
    
    % === Mean and standard deviation of segment-wise GC ===
    mean_AB = mean(GC_AB_seg, 'omitnan');
    std_AB  = std(GC_AB_seg, 0, 'omitnan');
    sem_AB  = std_AB / sqrt(sum(~isnan(GC_AB_seg)));
    
    mean_BA = mean(GC_BA_seg, 'omitnan');
    std_BA  = std(GC_BA_seg, 0, 'omitnan');
    sem_BA  = std_BA / sqrt(sum(~isnan(GC_BA_seg)));
    
    fprintf('\nGC %s→%s: mean = %.4f | std = %.4f | sem = %.4f\n', ...
        A_name, B_name, mean_AB, std_AB, sem_AB);
    fprintf('GC %s→%s: mean = %.4f | std = %.4f | sem = %.4f\n', ...
        B_name, A_name, mean_BA, std_BA, sem_BA);
    
    
    %% === Visualization ===
    
    A_label = 'ECoG Ch.33';
    B_label = 'STN';
    
    avg_sec = 4;  % <-- desired averaging window in seconds
    
    figure('Name','Segmented GC Summary','Position',[100 100 1100 900]);
    
    t_seg = (0:nSegments-1) * segment_length_sec;
    
    % (1) Segment-wise GC vs time (RAW)
    subplot(4,1,1);
    plot(t_seg, GC_AB_seg, 'b', 'LineWidth',1, ...
        'DisplayName', sprintf('%s→%s', A_label, B_label)); hold on;
    plot(t_seg, GC_BA_seg, 'r', 'LineWidth',1, ...
        'DisplayName', sprintf('%s→%s', B_label, A_label));
    
    xlabel('Segment start (s)');
    ylabel('GC value');
    legend('Location','best'); grid on;
    title(sprintf('Granger Causality over %.3g s windows', segment_length_sec));
    
    % (1b) GC averaged over avg_sec (NON-overlapping blocks)
    subplot(4,1,2);
    
    nPerBlock = round(avg_sec / segment_length_sec);
    nBlocks   = floor(nSegments / nPerBlock);
    
    GC_AB_blk = zeros(nBlocks,1);
    GC_BA_blk = zeros(nBlocks,1);
    t_blk     = zeros(nBlocks,1);
    
    for b = 1:nBlocks
        idx = (b-1)*nPerBlock + (1:nPerBlock);
        GC_AB_blk(b) = mean(GC_AB_seg(idx), 'omitnan');
        GC_BA_blk(b) = mean(GC_BA_seg(idx), 'omitnan');
        t_blk(b)     = mean(t_seg(idx));
    end
    
    plot(t_blk, GC_AB_blk, 'b', 'LineWidth',2, ...
        'DisplayName', sprintf('%s→%s', A_label, B_label)); hold on;
    plot(t_blk, GC_BA_blk, 'r', 'LineWidth',2, ...
        'DisplayName', sprintf('%s→%s', B_label, A_label));
    
    xlabel('Time (s)');
    ylabel('GC value');
    legend('Location','best'); grid on;
    title(sprintf('GC averaged over %.1f s (non-overlapping)', avg_sec));
    
    % (2) Mean 2×2 GC matrix (true average)
    subplot(4,1,3);
    imagesc(GC_mean);
    axis square; colormap(parula); colorbar;
    caxis([0 max(GC_mean(:))*1.1]);
    xlabel('Source'); ylabel('Target');
    set(gca,'XTick',1:2,'XTickLabel',{A_label, B_label}, ...
            'YTick',1:2,'YTickLabel',{A_label, B_label}, 'FontSize',9);
    title('Mean GC Matrix');
    
    % (3a) Bar comparison
    subplot(4,2,7);
    mean_AB = mean(GC_AB_seg,'omitnan');
    mean_BA = mean(GC_BA_seg,'omitnan');
    barVals = [mean_AB, mean_BA];
    b = bar(barVals, 0.5);
    b.FaceColor = 'flat';
    b.CData = [0 0 1; 1 0 0];
    
    set(gca,'XTickLabel',{sprintf('%s→%s',A_label,B_label), ...
                          sprintf('%s→%s',B_label,A_label)});
    ylabel('Mean GC');
    title(sprintf('Wilcoxon signed-rank: p=%.4g (z=%.3f)', p_val, zval));
    grid on;
    ylim([0, max(barVals)*1.25]);
    
    xt = 1:2;
    for i = 1:numel(barVals)
        text(xt(i), barVals(i) + 0.05*max(barVals), sprintf('%.4f', barVals(i)), ...
            'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
    end
    
    % (3b) Boxchart distribution
    subplot(4,2,8);
    xgroups = categorical([repmat({sprintf('%s→%s',A_label,B_label)}, nSegments, 1); ...
                           repmat({sprintf('%s→%s',B_label,A_label)}, nSegments, 1)]);
    yvals = [GC_AB_seg; GC_BA_seg];
    boxchart(xgroups, yvals, 'BoxWidth',0.6, 'MarkerStyle','o');
    ylabel('GC value');
    title('Segment-wise Distribution');
    grid on;
    
    
    
    %% Plotting Everything Together
    
    ECoG_ch = 'LS63';  
    
    signals = {"targetDistance", NO_to_BR_ch, BR_to_NO_ch, ECoG_ch};
    
    col_ECoG_to_STN = [0 0.2 0.6];   % dark blue
    col_STN_to_ECoG = [0.85 0 0];   % red
    
    titles  = { ...
        "Target Distance (thresholds shown)", ...
        sprintf("NO signal: %s", NO_to_BR_ch), ...
        sprintf("BR signal: %s", BR_to_NO_ch), ...
        sprintf("ECoG channel: %s", ECoG_ch) ...
        };
    
    colors  = {'k','b','r','m'};   % target, NO, BR, ECoG
    
    % --- Compute GC segment centers in time (seconds) ---
    t_seg = (0:nSegments-1) * segment_length_sec + segment_length_sec/2;
    
    % --- Set up figure with one extra subplot ---
    figure('Name','Synchronized Signals + GC','Position',[100 100 1200 900]);
    
    % --- (1) Synced signals ---
    for i = 1:numel(signals)
        subplot(numel(signals)+1, 1, i);
        plot(TT_longest_refined.Time, TT_longest_refined.(signals{i}), colors{i}, 'LineWidth', 1);
        title(titles{i}, 'Interpreter','none');
        ylabel(signals{i}, 'Interpreter','none');
        grid on;
        if i == 1
            yline(win(1),'--r'); yline(win(2),'--r');
        end
    end
    xlabel('Time');
    
    % --- (2) Granger Causality time series (bottom plot) ---
    A_name = Channel_Names{chIdx(1)};
    B_name = Channel_Names{chIdx(2)};
    
    subplot(numel(signals)+1, 1, numel(signals)+1);
    plot(t_seg, GC_AB_seg, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s→%s', A_name, B_name)); hold on;
    plot(t_seg, GC_BA_seg, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s→%s', B_name, A_name));
    xlabel('Segment start time (s)');
    ylabel('GC value');
    legend('Location','best');
    title(sprintf('Segment-wise Granger Causality (%s band, %d-s segments)', freq_band, segment_length_sec));
    grid on;
    
    
    
    
    
    
    
    %%
    
    fs = 1000;              % Hz
    order = 15;             % MVAR order
    chIdx = [56 69]; %21 42 61 62 63 65 74
    [GC, A, E, chNames] = compute_static_GC(TT_longest_refined, fs, order, chIdx, 'beta');
    
    %%
    fs = 1000;              % Hz
    order = 15;             % MVAR order
    chIdx = [56 69]; %21 42 61 62 63 65 74
    windowWidth = 100;
    
    [GC_all, GC_mean, GC_AB, GC_BA, seg_t] = ...
        compute_segmented_GC(TT_longest_refined, fs, order, chIdx, 'beta', windowWidth);
    
    
    
    %%
    win_lengths = [4, 5, 10, 15, 20, 25, 50, 100]; % seconds
    fs = 1000;
    chIdx = [56 69];
    order = 15;
    
    GC_AB_vals = zeros(size(win_lengths));  % mean A->B per window length
    GC_BA_vals = zeros(size(win_lengths));  % mean B->A per window length
    
    for i = 1:numel(win_lengths)
        wl = win_lengths(i);
        nsamp = wl * fs;
        nTotal = height(TT_longest_refined);
        nSeg = floor(nTotal / nsamp);
    
        if nSeg < 1
            warning('Skipping %.0f s window (too short)', wl);
            continue;
        end
    
        fprintf('\n===== %.0f-sec window (%d segments) =====\n', wl, nSeg);
        GC_AB_seg = zeros(nSeg,1);
        GC_BA_seg = zeros(nSeg,1);
    
        for k = 1:nSeg
            idx = (k-1)*nsamp + (1:nsamp);
            Tseg = TT_longest_refined(idx,:);
            [GC,~,~,~] = compute_static_GC(Tseg, fs, order, chIdx, 'beta', false);
    
            % GC(target, source)
            GC_AB_seg(k) = GC(2,1);  % A -> B
            GC_BA_seg(k) = GC(1,2);  % B -> A
            fprintf('  Segment %d/%d done\n', k, nSeg);
        end
    
        GC_AB_vals(i) = mean(GC_AB_seg,'omitnan');
        GC_BA_vals(i) = mean(GC_BA_seg,'omitnan');
    
        fprintf('  Avg GC (A->B)=%.3f | (B->A)=%.3f\n', GC_AB_vals(i), GC_BA_vals(i));
    end
    
    % --- Plot results ---
    figure;
    plot(win_lengths, GC_AB_vals, '-s', 'LineWidth', 1.5); hold on;
    plot(win_lengths, GC_BA_vals, '-s', 'LineWidth', 1.5);
    xlabel('Window length (s)');
    ylabel('Mean GC strength');
    title('Directional GC vs. Window Length (averaged across segments)');
    legend('A → B','B → A','Location','best');
    grid on;
    
    
    
    %%
    
    fs = 1000;
    order = 12;
    lags = -30:1:30; % test ±30 ms in 1 ms steps
    sweep_gc_vs_lag(TT_longest_refined, fs, order, "LS63", "CMacro_LFP_AP_T", lags, 'beta');
    
    
    
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
    
    ft_defaults
    
    %% FieldTrip Visualizations
    
    Patient_Folder_Path = 'C:\Users\cours\Downloads\AndersonLabDataAndPackages\Patient Data';
    
    % --- Check for lh.pial under Imaging/freesurfer/surf ---
    %lh_pial_default = fullfile(Patient_Folder_Path, Patient_ID, 'freesurfer', 'surf', 'lh.pial');
    %lh_pial_path = fullfile(Patient_Folder_Path, Patient_ID, 'Imaging', 'freesurfer', 'surf', 'lh.pial.T1'); %%rh vs lh
    
    %lh_pial_path = fullfile(Patient_Folder_Path, Patient_ID, 'Imaging', 'fsaverage', 'surf', 'lh.pial'); %%rh vs lh
    
    
    %rh_pial_default = fullfile(Patient_Folder_Path, Patient_ID, 'freesurfer', 'surf', 'rh.pial');
    %rh_pial_path = fullfile(Patient_Folder_Path, Patient_ID, 'Imaging', 'freesurfer', 'surf', 'rh.pial.T1');
    
    %rh_pial_path = fullfile(Patient_Folder_Path, Patient_ID, 'Imaging', 'fsaverage', 'surf', 'rh.pial');
    
    
    [ftver, ftpath] = ft_version;
    
    lh_pial_mni_path = fullfile(ftpath, 'template', 'anatomy', 'surface_pial_left.mat');
    rh_pial_mni_path = fullfile(ftpath, 'template', 'anatomy', 'surface_pial_right.mat');
    
    
    % if exist(lh_pial_default,'file')
    %     lh_pial_path = lh_pial_default;
    % else
    %     error('No lh.pial found.');
    % end
    
    % if exist(rh_pial_default,'file')
    %     rh_pial_path = rh_pial_default;
    % else
    %     error('No rh.pial found.');
    % end
    
    
    % --- ACPC electrode file ---
    acpc_imaging = fullfile(Patient_Folder_Path, Patient_ID, 'Imaging', [Patient_ID '_elec_acpc_fr.mat']);
    acpc_default = fullfile(Patient_Folder_Path, Patient_ID, ['Subject' Patient_ID '_elec_acpc_fr.mat']);
    
    if exist(acpc_imaging, 'file')
        ACPC_Path = acpc_imaging;
    elseif exist(acpc_default, 'file')
        ACPC_Path = acpc_default;
    else
        error('No ACPC electrode file found for %s.', Patient_ID);
    end
    
    % --- FSAVG electrode file ---
    %FSAVG_Path = fullfile(Patient_Folder_Path, Patient_ID, ['Subject' Patient_ID '_elec_fsavg_frs.mat']);
    % ACPC_Path = fullfile(Patient_Folder_Path, Patient_ID, [Patient_ID '_elec_acpc_fr.mat']);
    FS_PATH = fullfile(Patient_Folder_Path, Patient_ID, 'Imaging', [Patient_ID '_elec_fsavg_frs.mat']);
    MNI_FRV_Path = fullfile(Patient_Folder_Path, Patient_ID, 'Imaging', [Patient_ID '_elec_mni_frv.mat']);
    
    % disp('lh_pial_path:'), disp(lh_pial_path)
    % disp('ACPC_Path:'), disp(ACPC_Path)
    % disp('FSAVG_Path:'), disp(FSAVG_Path)
    
    % --- Load meshes dynamically ---
    disp('--- Parameter sanity check ---')
    disp(['ACPC file exists: ', string(exist(ACPC_Path,'file')==2)]);
    % disp(['Pial file exists: ', string(exist(lh_pial_path,'file')==2)]);
    
    elec = load(ACPC_Path).elec_acpc_fr;
    disp(['Electrode count: ', num2str(size(elec.chanpos,1))]);
    disp(['Data vector length: ', num2str(length(PLV_Relative_ref(1:63)))]);
    
    %% Plot PLV wrt Macro LFP on ACPC FR Mesh
    
    data = PLV_Relative_ref(1:63);
    data(48) = 0;
    
    %load_ACPC_FR_mesh(FS_PATH, rh_pial_path, data, []);
    
    %load_ACPC_FR_mesh(ACPC_Path, lh_pial_path, data);
    %load_FSAVG_FRS_mesh(FS_PATH, rh_pial_path, data);
    
    load_MNI_FRV_mesh(MNI_FRV_Path, rh_pial_mni_path, data, [] );
    
    
    
    %%
    
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
    
    
    
    
    function TT_window = extractNthLongestWindow(TT_master, sigName, fs, win, trimFrontSec, trimBackSec, nSelect, subplotChan)
    %EXTRACTNTHLONGESTWINDOW Keep the n-th longest contiguous window within threshold,
    %   with optional separate trimming of front and back, and optional subplot overlay.
    %
    %   TT_window = extractNthLongestWindow(TT_master, sigName, fs, win, trimFrontSec, trimBackSec, nSelect, subplotChan)
    %
    %   Inputs:
    %     TT_master    - timetable (full data)
    %     sigName      - name of signal variable to threshold (string/char)
    %     fs           - sampling rate (Hz)
    %     win          - [lower upper] threshold window
    %     trimFrontSec - seconds to cut from the beginning of the window (default = 0)
    %     trimBackSec  - seconds to cut from the end of the window (default = 0)
    %     nSelect      - which longest window to select (1 = longest, 2 = 2nd longest, etc.)
    %     subplotChan  - (optional) another channel name to plot in a subplot above
    %
    %   Output:
    %     TT_window - timetable (subset of TT_master)
    
        if nargin < 5, trimFrontSec = 0; end
        if nargin < 6, trimBackSec  = 0; end
        if nargin < 7, nSelect = 1; end
        if nargin < 8, subplotChan = ''; end  % optional extra channel
    
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
    
        % --- Step 6: Plot ---
        t = TT_master.Properties.RowTimes;
    
        if ~isempty(subplotChan) && ismember(subplotChan, TT_master.Properties.VariableNames)
            % Two subplots: upper = subplotChan, lower = sigName
            figure;
            tiledlayout(2,1);
            
            % Top plot
            nexttile;
            plot(t, TT_master.(subplotChan), 'b'); grid on;
            ylabel(subplotChan, 'Interpreter', 'none');
            title(sprintf('%d-th longest window overlay: %s & %s', nSelect, subplotChan, sigName));
            xlim([t(1) t(end)]);
            
            % Bottom plot
            nexttile;
            plot(t, sig, 'k'); hold on; grid on;
            yline(win(1), '--r','Lower'); yline(win(2), '--r','Upper');
            yl = ylim;
            patch([t(S) t(E) t(E) t(S)], [yl(1) yl(1) yl(2) yl(2)], ...
                  [0.85 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');
            uistack(findobj(gca,'Type','line'),'top');
            xlabel('Time'); ylabel(sigName, 'Interpreter', 'none');
            xlim([t(1) t(end)]);
        else
            % Single plot
            figure;
            plot(t, sig, 'k'); hold on; grid on;
            yline(win(1), '--r','Lower'); yline(win(2), '--r','Upper');
            yl = ylim;
            patch([t(S) t(E) t(E) t(S)], [yl(1) yl(1) yl(2) yl(2)], ...
                  [0.85 0.3 0.3],'FaceAlpha',0.2,'EdgeColor','none');
            uistack(findobj(gca,'Type','line'),'top');
            xlabel('Time'); ylabel(sigName, 'Interpreter', 'none');
            title(sprintf('%d-th longest segment = %.2f s (trim = %.2f s front, %.2f s back)', ...
                nSelect, durSorted(nSelect), trimFrontSec, trimBackSec));
        end
    end
    
    
    
    
    
    
    
    function load_ACPC_FR_mesh(acpc_path, pial_lh_path, data, colorRange, pial_rh_path)
        % load_ACPC_FR_mesh(acpc_path, pial_lh_path, data, colorRange, pial_rh_path)
        %
        % pial_rh_path is OPTIONAL.
        % If provided and the file exists, the function will also plot the RH pial
        % surface using the SAME interpolated source data.
    
        % === Handle optional argument ===
        if nargin < 5
            pial_rh_path = [];
        end
    
        if nargin < 4 || isempty(colorRange)
            colorRange = [min(data(:)) max(data(:))];
        end
    
        % === Step 1: Load electrode positions and mesh (LH) ===
        elec_acpc_fr = load(acpc_path).elec_acpc_fr;
        pial_lh = ft_read_headshape(pial_lh_path);
        pial_lh.coordsys = 'acpc';
    
        pial_lh = ft_convert_units(pial_lh, 'mm');
        elec_acpc_fr = ft_convert_units(elec_acpc_fr, 'mm');
    
        % === Step 2: Handle NaNs ===
        ref_idx = isnan(data);
        data_clean = data;
        data_clean(ref_idx) = 0;
    
        % === Step 3: Normalize data ===
        norm_data = (data_clean - colorRange(1)) ./ (colorRange(2) - colorRange(1));
        norm_data = max(min(norm_data, 1), 0);
    
        % === Step 4: Build source structure ===
        source = struct();
        source.pos = elec_acpc_fr.chanpos;
        source.pow = norm_data;
        source.dimord = 'pos';
        source.coordsys = 'acpc';
    
        % === Step 5: Interpolate electrode data onto LH ===
        cfg = [];
        cfg.parameter = 'pow';
        cfg.interpmethod = 'sphere_weighteddistance';
        cfg.sphereradius = 3;
        interp_source = ft_sourceinterpolate(cfg, source, pial_lh);
        cfg_interp = cfg; 
    
        % === Step 6: Plot interpolated data on LH ===
        cfg = [];
        cfg.funparameter = 'pow';
        cfg.funcolorlim = [0 1];
        cfg.method = 'surface';
        cfg.interpmethod = 'nearest';
        cfg.sphereradius = 3;
        cfg.camlight = 'no';
        ft_sourceplot(cfg, interp_source, pial_lh);
        hold on;   % <-- allow plotting RH next if applicable
    
        % =======================================================
        % === OPTIONAL: Load & plot RH pial (NO other changes) ===
        % =======================================================
        % OPTIONAL: Load & plot RH pial
        if ~isempty(pial_rh_path) && exist(pial_rh_path, 'file')
            pial_rh = ft_read_headshape(pial_rh_path);
            pial_rh.coordsys = 'acpc';
            pial_rh = ft_convert_units(pial_rh, 'mm');
        
            % Re-use interpolation cfg
            interp_source_rh = ft_sourceinterpolate(cfg_interp, source, pial_rh);
        
            % Plot RH surface using SAME plotting cfg as LH
            ft_sourceplot(cfg, interp_source_rh, pial_rh);
        end
    
        % =======================================================
    
        view([-55 10]);
        material dull;
        lighting gouraud;
        camlight;
        colormap(parula);
    
        % === Step 7: Plot electrodes ===
        elec_plot_colors = norm_data;
        ft_plot_sens(elec_acpc_fr, 'elecsize', 35, ...
            'facecolor', elec_plot_colors, 'edgecolor', 'black');
    
        % === Step 8: Overlay reference electrodes ===
        if any(ref_idx)
            ref_positions = elec_acpc_fr.chanpos(ref_idx, :);
            scatter3(ref_positions(:,1), ref_positions(:,2), ref_positions(:,3), ...
                100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        end
    
        % === Step 9: Add colorbar ===
        cb = colorbar;
        tick_vals = linspace(0, 1, 5);
        cb.Ticks = tick_vals;
        cb.TickLabels = string(linspace(colorRange(1), colorRange(2), 5));
        ylabel(cb, 'PLV value');
    end
    
    
    
    
    function load_FSAVG_FRS_mesh(fsavg_path, pial_lh_path, data, colorRange, pial_rh_path)
    % load_FSAVG_FRS_mesh(fsavg_path, pial_lh_path, data, colorRange, pial_rh_path)
    %
    % Mirrors load_ACPC_FR_mesh but for FSAVG (fsaverage) space.
    % - fsavg_path: path to *_elec_fsavg_frs.mat
    % - pial_lh_path: path to fsaverage lh.pial
    % - data: vector with one value per electrode (same length as chanpos)
    % - colorRange: [min max] for mapping data -> colormap (optional)
    % - pial_rh_path: OPTIONAL fsaverage rh.pial for RH plot
    
        % === Handle optional args ===
        if nargin < 5
            pial_rh_path = [];
        end
    
        if nargin < 4 || isempty(colorRange)
            colorRange = [min(data(:)) max(data(:))];
        end
    
        % === Step 1: Load FSAVG electrodes and LH pial mesh ===
        S = load(fsavg_path);
        if ~isfield(S, 'elec_fsavg_frs')
            error('File %s does not contain elec_fsavg_frs.', fsavg_path);
        end
        elec_fsavg = S.elec_fsavg_frs;
    
        pial_lh = ft_read_headshape(pial_lh_path);
        pial_lh.coordsys = 'fsaverage';
    
        elec_fsavg = ft_convert_units(elec_fsavg, 'mm');
        pial_lh    = ft_convert_units(pial_lh, 'mm');
    
        % === Step 2: Sanity check on data length ===
        nElec = size(elec_fsavg.chanpos, 1);
        data  = data(:);  % ensure column
    
        if numel(data) ~= nElec
            error('Data length (%d) does not match number of FSAVG electrodes (%d).', ...
                  numel(data), nElec);
        end
    
        % === Step 3: Handle NaNs (reference / bad channels) ===
        ref_idx    = isnan(data);
        data_clean = data;
        data_clean(ref_idx) = 0;
    
        % === Step 4: Normalize to [0,1] ===
        norm_data = (data_clean - colorRange(1)) ./ (colorRange(2) - colorRange(1));
        norm_data = max(min(norm_data, 1), 0);
    
        % === Step 5: Build FieldTrip source structure at electrode locations ===
        source = struct();
        source.pos      = elec_fsavg.chanpos;
        source.pow      = norm_data;
        source.dimord   = 'pos';
        source.coordsys = 'fsaverage';
    
        % === Step 6: Interpolate onto LH pial surface ===
        cfg_interp = [];
        cfg_interp.parameter    = 'pow';
        cfg_interp.interpmethod = 'sphere_weighteddistance';
        cfg_interp.sphereradius = 3;
    
        interp_source_lh = ft_sourceinterpolate(cfg_interp, source, pial_lh);
    
        % === Step 7: Plot config ===
        cfg_plot = [];
        cfg_plot.funparameter = 'pow';
        cfg_plot.funcolorlim  = [0 1];
        cfg_plot.method       = 'surface';
        cfg_plot.interpmethod = 'nearest';
        cfg_plot.sphereradius = 3;
        cfg_plot.camlight     = 'no';
    
        % === Step 8: Plot LH ===
        ft_sourceplot(cfg_plot, interp_source_lh);
        hold on;
    
        % === Step 9: OPTIONAL: RH pial ===
        if ~isempty(pial_rh_path) && exist(pial_rh_path, 'file')
            pial_rh = ft_read_headshape(pial_rh_path);
            pial_rh.coordsys = 'fsaverage';
            pial_rh = ft_convert_units(pial_rh, 'mm');
    
            interp_source_rh = ft_sourceinterpolate(cfg_interp, source, pial_rh);
            ft_sourceplot(cfg_plot, interp_source_rh);
        end
    
        % === Step 10: View, lighting, colormap ===
        view([-55 10]);
        material dull;
        lighting gouraud;
        camlight;
        colormap(parula);
    
        % === Step 11: Plot electrodes on top ===
        elec_plot_colors = norm_data;
        ft_plot_sens(elec_fsavg, 'elecsize', 35, ...
            'facecolor', elec_plot_colors, 'edgecolor', 'black');
    
        % === Step 12: Overlay reference/bad electrodes in red ===
        if any(ref_idx)
            ref_positions = elec_fsavg.chanpos(ref_idx, :);
            scatter3(ref_positions(:,1), ref_positions(:,2), ref_positions(:,3), ...
                100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        end
    
        % === Step 13: Colorbar in original data units ===
        cb = colorbar;
        tick_vals = linspace(0, 1, 5);
        cb.Ticks = tick_vals;
        cb.TickLabels = string(linspace(colorRange(1), colorRange(2), 5));
        ylabel(cb, 'PLV value');
    
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
                    title_str = sprintf('PLV Map with respect to STN (Raw)');
                    label_val = PLV_relative_vector(i);
                else
                    % Normalized PLV scaling
                    color_val = plv_norm(i);
                    title_str = sprintf('PLV Map with respect to STN (Normalized)');
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
    
    
    function plot_PLV_Relative_Brainnetome(phase_data, ref_phase, PLV_relative_vector, ref_channel, channelNames, meta, badCh)
        % Plot PLV relative to a reference channel on a 21x3 grid
        % Face color = PLV (parula colormap)
        % Border color = Brainnetome region (high-contrast categorical palette)
        % Also computes region-level PLV statistics
    
        numChannels = length(phase_data);
        ref_name = channelNames{ref_channel};
    
        if nargin < 7
            badCh = [];   % default: no bad channels
        end
    
        % Set bad channels to 0 in the PLV vector
        PLV_relative_vector(badCh) = 0;
    
        % Normalize PLV values
        PLV_copy = PLV_relative_vector;
        PLV_copy(ref_channel) = -1;
        PLV_copy(64:end) = -1;
        max_plv = max(PLV_copy);
    
        disp(max_plv)
    
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
    
        plv_norm = PLV_relative_vector;
        plv_norm(isnan(plv_norm)) = 0;
        plv_norm = plv_norm / max_plv;
    
        % === Brainnetome mapping (numeric-based) ===
        % === Brainnetome mapping (numeric-based with zero-padding) ===
        cleanLabels = regexprep(string(meta.Brainnetome), '\s*\(\d+%?\)\s*$', '');
        
        % Extract electrode numbers from meta (works for "1", "L01", "LS01")
        metaNums = regexp(string(meta.Electrode), '\d+', 'match', 'once');
        metaNums = string(pad(metaNums, 2, 'left', '0'));   % "1" -> "01"
        
        % Build map: electrode number -> Region
        brain_map = containers.Map(metaNums, cleanLabels);
    
        regions = unique(cleanLabels);
    
        % High-contrast categorical palette (ColorBrewer Set1 style)
        baseColors = [
            228,26,28;    % red
            55,126,184;   % blue
            77,175,74;    % green
            152,78,163;   % purple
            255,127,0;    % orange
            166,86,40;    % brown
            247,129,191;  % pink
            153,153,153;  % grey
        ] / 255;
    
        regionColors = baseColors(mod(0:numel(regions)-1, size(baseColors,1))+1, :);
        regionColor  = containers.Map(regions, num2cell(regionColors,2));
    
        legendHandles = gobjects(numel(regions),1);
    
        % ---- Manual subplot positions ----
        pos1 = [0.05, 0.1, 0.4, 0.8];  % left PLV map
        pos2 = [0.25, 0.1, 0.4, 0.8];
        figure;
    
        for subplot_idx = 1:2
            if subplot_idx == 1
                ax = axes('Position', pos1);
                color_vals = PLV_relative_vector;
                title_str = sprintf('PLV Map with respect to STN (raw)');
                labels = PLV_relative_vector;
            else
                ax = axes('Position', pos2);
                color_vals = plv_norm;
                title_str = sprintf('PLV Map with respect to STN (normalized)');
                labels = plv_norm;
            end
            hold(ax,'on');
    
            % Draw electrodes
            for i = 1:numChannels
                x_i = electrode_positions{i}(1);
                y_i = electrode_positions{i}(2);
    
                % Fill color
                color_idx = round(color_vals(i) * 255) + 1;
                color_idx = min(max(color_idx,1),256);
                faceCol   = cmap(color_idx,:);
    
                % Border color (use numeric channel match)
                chNum = regexp(string(channelNames{i}), '\d+', 'match', 'once');
                if isKey(brain_map, chNum)
                    region  = brain_map(chNum);
                    edgeCol = regionColor(region);
    
                    % Save legend proxy
                    idx = find(strcmp(regions, region));
                    if isempty(legendHandles(idx)) || ~isgraphics(legendHandles(idx))
                        legendHandles(idx) = plot(ax, nan,nan,'s', ...
                            'MarkerSize',10,'MarkerEdgeColor',edgeCol, ...
                            'MarkerFaceColor','w','LineWidth',2);
                    end
                else
                    edgeCol = 'k';
                    fprintf('No match for channel: "%s"\n', channelNames{i});
                    region = "Unlabeled";
                end
    
                % Electrode square
                rectangle('Position',[x_i,y_i,1,1], ...
                          'FaceColor',faceCol,'EdgeColor',edgeCol,'LineWidth',2, ...
                          'Parent',ax);
    
                % Labels
                text(x_i+0.5,y_i+0.65,string(channelNames{i}), ...
                     'HorizontalAlignment','center','FontSize',7,'FontWeight','bold','Parent',ax);
                text(x_i+0.5,y_i+0.35,sprintf('%.3f',labels(i)), ...
                     'HorizontalAlignment','center','FontSize',6,'Parent',ax);
            end
    
            axis(ax,'equal'); xlim(ax,[0,grid_width]); ylim(ax,[0,grid_height]);
            set(ax,'XTick',[],'YTick',[],'YDir','normal');
            title(ax,title_str, 'Interpreter','none');
            colormap(ax,cmap);
            cbar = colorbar(ax);
    
            
            %ylabel(cbar, sprintf('PLV wrt %s', ref_name), 'Interpreter','none');
    
            raw_min = min(PLV_relative_vector(1:63));
            raw_max = max(PLV_relative_vector(1:63));
    
            disp(raw_min)
            disp(raw_max)
            
            % Tell MATLAB that the color range represents the RAW PLV range
            caxis(ax, [raw_min raw_max]);
            
            % Create colorbar and label with raw PLV tick labels
            cbar = colorbar(ax);
            cbar.Ticks = linspace(raw_min, raw_max, 5);
            cbar.TickLabels = arrayfun(@(v) sprintf('%.3f', v), cbar.Ticks, 'UniformOutput', false);
            
            ylabel(cbar, sprintf('PLV wrt %s', ref_name), 'Interpreter','none');
        end
    
        % ---- Legend outside far right ----
        validIdx = isgraphics(legendHandles);
        lgd = legend(legendHandles(validIdx), regions(validIdx), ...
            'Interpreter','none','Orientation','vertical');
        lgd.Position = [0.6, 0.2, 0.15, 0.6];
    
        % ---- Region statistics ----
        PLV_relative_vector = PLV_relative_vector(1:numChannels);
    
        chanRegions = strings(numChannels,1);
        for i = 1:numChannels
            chNum = regexp(string(channelNames{i}), '\d+', 'match', 'once');
            if isKey(brain_map, chNum)
                chanRegions(i) = brain_map(chNum);
            else
                chanRegions(i) = "Unlabeled";
            end
        end
    
        % Build table and compute average PLV per region
        T = table(chanRegions, PLV_relative_vector, 'VariableNames', {'Region','PLV'});
        regionPLVs = groupsummary(T, "Region", "mean");
        regionPLVs.Properties.VariableNames = {'Region','GroupCount','MeanPLV'};
    
        
        regionPLVs = sortrows(regionPLVs, "MeanPLV", "descend");
    
        [sortedVals, sortedIdx] = sort(PLV_relative_vector, "descend");
        top5ch  = channelNames(sortedIdx(1:5));
        top5val = sortedVals(1:5);
        top5reg = chanRegions(sortedIdx(1:5));
    
        fprintf('\n=== Region Ranking (by Avg PLV) ===\n');
        disp(regionPLVs);
    
        fprintf('\n=== Top 5 Channels by PLV ===\n');
        for k = 1:5
            fprintf('Ch %s: PLV = %.3f, Region = %s\n', top5ch{k}, top5val(k), top5reg(k));
        end
    end
    
    
    function [GC_matrix, A, E, chNames] = compute_static_GC(TT, fs, order, chIdx, freq_range, plotFlag, doFilter)
    % compute_static_GC  Estimate pairwise (static) Granger causality with optional band filtering and plotting.
    %
    %   [GC_matrix, A, E, chNames] = compute_static_GC(TT, fs, order, chIdx, freq_range, plotFlag, doFilter)
    %
    % Inputs:
    %   TT         - timetable, table, or numeric matrix (channels × samples) or (samples × channels) if TT is timetable/table
    %   fs         - sampling frequency (Hz)
    %   order      - AR model order (e.g., 5–20)
    %   chIdx      - optional vector of channel indices (default = all)
    %   freq_range - 'beta','theta','alpha','none' (default='none')
    %   plotFlag   - true/false (default=true)
    %   doFilter   - true/false: if true, filter INSIDE this function (default=false)
    %
    % Outputs:
    %   GC_matrix  - nCh × nCh matrix (GC(i,j): j→i)
    %   A, E       - AR coefficients and residual covariance from full model
    %   chNames    - channel names (if available)
    
        if nargin < 5 || isempty(freq_range)
            freq_range = 'none';
        end
        if nargin < 6 || isempty(plotFlag)
            plotFlag = true;
        end
        if nargin < 7 || isempty(doFilter)
            doFilter = false;  % <-- NEW DEFAULT: do NOT filter unless explicitly requested
        end
    
        % --- Extract numeric data and names ---
        if istimetable(TT)
            X = TT.Variables';                      % channels × samples
            allNames = string(TT.Properties.VariableNames);
        elseif istable(TT)
            X = table2array(TT)';                   % channels × samples
            allNames = string(TT.Properties.VariableNames);
        else
            X = TT;
            if size(X,1) > size(X,2)
                % If user passed samples×channels accidentally, flip to channels×samples
                % (optional safeguard; remove if you don't want it)
                % X = X';
            end
            allNames = "Ch" + string(1:size(X,1));
        end
    
        % --- Channel selection ---
        nAll = size(X,1);
        if nargin < 4 || isempty(chIdx)
            chIdx = 1:nAll;
        end
        X = X(chIdx,:);
        chNames = allNames(chIdx);
        nCh = numel(chIdx);
        nSamp = size(X,2);
    
        % --- Optional band-filter (explicitly controlled) ---
        if doFilter
            if strcmpi(freq_range, 'none')
                % user asked to filter but gave 'none' -> treat as no filtering
            else
                filtOrder = 1024;
    
                % Myeegfilt uses filtfilt => requires epochframes >= 3*filtorder
                minLen = 3 * filtOrder;
                if nSamp < minLen
                    error(['compute_static_GC:SegmentTooShortForFilter\n' ...
                           'Requested %s filtering with filtOrder=%d, but nSamp=%d < 3*filtOrder=%d.\n' ...
                           'Fix: (1) prefilter outside compute_static_GC, then call with doFilter=false,\n' ...
                           'or (2) increase window length, or (3) reduce filtOrder.\n'], ...
                          freq_range, filtOrder, nSamp, minLen);
                    % If you prefer warn+skip instead of error, replace the error(...) with:
                    % warning('Segment too short for FIR(%d). Skipping filter.', filtOrder);
                end
    
                fprintf('Applying %s-band filter inside compute_static_GC...\n', freq_range);
                for ch = 1:nCh
                    ch_raw = double(X(ch,:));
                    switch lower(freq_range)
                        case 'theta'
                            X(ch,:) = Myeegfilt(ch_raw, fs, 4, 9, 0, filtOrder);
                        case 'beta'
                            X(ch,:) = Myeegfilt(ch_raw, fs, 13, 30, 0, filtOrder);
                        case 'alpha'
                            X(ch,:) = Myeegfilt(ch_raw, fs, 8, 12, 0, filtOrder);
                        otherwise
                            error('Unknown frequency range: %s', freq_range);
                    end
                end
                fprintf('Filtering complete.\n');
            end
        end
    
        % --- Preprocess ---
        X = X - mean(X,2);
        X = X ./ std(X,0,2);
    
        keep = var(X,0,2) > 1e-6 & all(isfinite(X),2);
        X = X(keep,:);
        chNames = chNames(keep);
        nCh = sum(keep);
    
        if nCh < 2
            error('compute_static_GC:NotEnoughValidChannels', ...
                  'After preprocessing, fewer than 2 valid channels remain.');
        end
    
        % --- Fit full model ---
        fprintf('Fitting full MVAR model (%d channels, order %d)...\n', nCh, order);
        t0 = tic;
        try
            [A,E] = armorf(X, 1, nSamp, order);
        catch
            warning('armorf failed; retrying with small jitter...');
            X = X + 1e-6*randn(size(X));
            [A,E] = armorf(X, 1, nSamp, order);
        end
        fprintf('Full model fit done in %.2f s.\n', toc(t0));
    
        % --- GC computation ---
        fprintf('Computing pairwise GC (%d sources)...\n', nCh);
        GC_matrix = zeros(nCh,nCh);
        fprintf('Progress: ');
        for j = 1:nCh
            idx = setdiff(1:nCh, j);
            Xr = X(idx,:);
            [~,E_r] = armorf(Xr, 1, nSamp, order);
    
            for iIdx = 1:numel(idx)
                i = idx(iIdx);
                GC_matrix(i,j) = log(E_r(iIdx,iIdx) / E(i,i));
            end
    
            if mod(j, max(1, round(nCh/20))) == 0 || j==nCh
                fprintf('|');
            end
        end
        fprintf(' done.\n');
    
        % --- Visualization (optional) ---
        if plotFlag
            figure;
            imagesc(GC_matrix);
            axis square;
            colormap(parula);
            colorbar;
            caxis([0 max(GC_matrix(:))]);
            xlabel('Source (causal) channel');
            ylabel('Target (influenced) channel');
            title(sprintf('Static Granger Causality (%s band, order=%d)', freq_range, order));
            set(gca, 'XTick', 1:nCh, 'XTickLabel', chNames, ...
                     'YTick', 1:nCh, 'YTickLabel', chNames, ...
                     'XTickLabelRotation', 90);
        end
    
        fprintf('Total runtime: %.2f s\n', toc(t0));
    end
    
    
    
    function sweep_gc_vs_lag(TT, fs, order, chA, chB, lag_ms_vec, freq_range)
    % sweep_gc_vs_lag  Slide channel B relative to A and compute bivariate GC vs lag.
    %
    % Inputs
    %   TT            : timetable with signals (rows=Time, vars=channels)
    %   fs            : sampling rate (Hz)
    %   order         : AR order for GC
    %   chA, chB      : channel variable names (strings) or indices
    %   lag_ms_vec    : vector of lags in ms (e.g., -20:1:20)
    %   freq_range    : 'none'|'beta'|'theta' (uses your Myeegfilt if set)
    %
    % Notes
    %   * Positive lag means **B is delayed** by that amount (A leads).
    %   * Uses linear interpolation (sub-sample accurate), trims edges equally.
    %   * Computes bivariate GC only on the (A,B) pair for speed/clarity.
    %   * Only one summary plot is shown (no intermediate GC matrix plots).
    
    % --- Extract signals and names ---
    if isstring(chA) || ischar(chA)
        aSig = TT.(chA); aName = string(chA);
    else
        aSig = TT.Variables(:,chA);
        aName = string(TT.Properties.VariableNames{chA});
    end
    
    if isstring(chB) || ischar(chB)
        bSig = TT.(chB); bName = string(chB);
    else
        bSig = TT.Variables(:,chB);
        bName = string(TT.Properties.VariableNames{chB});
    end
    
    aSig = double(aSig(:)); 
    bSig = double(bSig(:));
    t    = (0:numel(aSig)-1)'/fs;
    
    % --- Optional bandpass filter ---
    if nargin >= 7 && ~strcmpi(freq_range,'none')
        switch lower(freq_range)
            case 'beta'
                aSig = Myeegfilt(aSig', fs, 13, 30, 0, 1024)'; 
                bSig = Myeegfilt(bSig', fs, 13, 30, 0, 1024)';
            case 'theta'
                aSig = Myeegfilt(aSig', fs, 4, 9, 0, 1024)';  
                bSig = Myeegfilt(bSig', fs, 4, 9, 0, 1024)';
            otherwise
                error('Unknown freq_range');
        end
    end
    
    % --- Preallocate ---
    GC_AB = zeros(size(lag_ms_vec));
    GC_BA = zeros(size(lag_ms_vec));
    GC_sum = zeros(size(lag_ms_vec));
    ASYM   = zeros(size(lag_ms_vec));
    
    fprintf('Sweeping lags (%d total)...\n', numel(lag_ms_vec));
    
    for k = 1:numel(lag_ms_vec)
        lag_s = lag_ms_vec(k) / 1000;   % convert ms to seconds
        % Shift B by +lag (delay B)
        tB = t + lag_s;
        bShift = interp1(tB, bSig, t, 'linear', NaN);
    
        % Trim NaNs so signals are same length
        valid = ~isnan(bShift);
        a = aSig(valid); 
        b = bShift(valid);
    
        % Z-score normalization
        a = (a - mean(a)) / std(a);
        b = (b - mean(b)) / std(b);
    
        % Build temporary timetable
        Tpair = array2timetable([a b], ...
            'RowTimes', seconds((0:numel(a)-1)'/fs), ...
            'VariableNames', [aName, bName]);
    
        % Run GC silently (no matrix plotting)
        [GC, ~, ~, ~] = compute_static_GC(Tpair, fs, order, [], 'none', false);
    
        % Extract relevant pairwise values
        GC_AB(k) = GC(2,1);  % A → B
        GC_BA(k) = GC(1,2);  % B → A
        GC_sum(k)= GC_AB(k) + GC_BA(k);
        if GC_sum(k) > 0
            ASYM(k) = (GC_AB(k) - GC_BA(k)) / GC_sum(k);
        else
            ASYM(k) = 0;
        end
    
        if mod(k, round(numel(lag_ms_vec)/10)) == 0 || k == numel(lag_ms_vec)
            fprintf('  completed %d / %d lags...\n', k, numel(lag_ms_vec));
        end
    end
    
    fprintf('Sweep complete. Generating summary plot...\n');
    
    % --- Final summary plot ---
    figure('Name','GC vs lag');
    subplot(3,1,1);
    plot(lag_ms_vec, GC_AB, '-o', lag_ms_vec, GC_BA, '-s'); grid on;
    xlabel('Lag applied to B (ms)'); ylabel('GC');
    legend(sprintf('%s→%s', aName, bName), sprintf('%s→%s', bName, aName), 'Location','best');
    title('Directional GC vs lag');
    
    subplot(3,1,2);
    plot(lag_ms_vec, GC_sum, '-o'); grid on;
    xlabel('Lag applied to B (ms)'); ylabel('GC(A→B)+GC(B→A)');
    title('Total coupling strength (direction-agnostic)');
    
    subplot(3,1,3);
    plot(lag_ms_vec, ASYM, '-o'); yline(0,'k-'); grid on;
    xlabel('Lag applied to B (ms)'); ylabel('Asymmetry index');
    title('Asymmetry = (A→B − B→A) / (sum)');
    
    % --- Summary printout ---
    [~,kMaxSum] = max(GC_sum);
    [~,kMaxAB]  = max(GC_AB);
    [~,kMaxBA]  = max(GC_BA);
    fprintf('Max total strength at lag %+d ms\n', lag_ms_vec(kMaxSum));
    fprintf('Max %s→%s at lag %+d ms\n', aName, bName, lag_ms_vec(kMaxAB));
    fprintf('Max %s→%s at lag %+d ms\n', bName, aName, lag_ms_vec(kMaxBA));
    fprintf('Done.\n');
    end
    
    
    
    
    function [GC_all, GC_mean, GC_AB, GC_BA, seg_times] = compute_segmented_GC(TT, fs, order, chIdx, freq_range, seg_len_sec)
    % compute_segmented_GC
    %   Computes Granger causality over multiple short segments and aggregates results.
    %
    % Inputs:
    %   TT           - timetable or table of data (Time x Channels)
    %   fs           - sampling frequency (Hz)
    %   order        - AR model order
    %   chIdx        - vector of channel indices (e.g. [63 74])
    %   freq_range   - 'beta','theta','none'
    %   seg_len_sec  - length of each segment in seconds
    %
    % Outputs:
    %   GC_all   - 3D array [nCh x nCh x nSegments]
    %   GC_mean  - average GC matrix across segments
    %   GC_AB    - vector of GC(A→B) per segment (if 2 channels)
    %   GC_BA    - vector of GC(B→A) per segment
    %   seg_times - start time of each segment (in seconds)
    %
    % Example:
    %   [GC_all,GC_mean,GC_AB,GC_BA,t] = compute_segmented_GC(TT_longest_refined,1000,10,[63 74],'beta',5);
    
    % --- setup
    X = TT.Variables';
    nCh = numel(chIdx);
    seg_len_samples = seg_len_sec * fs;
    nSamples = size(X, 2);
    nSeg = floor(nSamples / seg_len_samples);
    
    fprintf('Computing GC on %d segments of %.1f s each...\n', nSeg, seg_len_sec);
    
    GC_all = zeros(nCh, nCh, nSeg);
    GC_AB = zeros(1, nSeg);
    GC_BA = zeros(1, nSeg);
    seg_times = (0:nSeg-1) * seg_len_sec;
    
    % --- loop over segments
    for k = 1:nSeg
        idx = (k-1) * seg_len_samples + (1:seg_len_samples);
        Tseg = TT(idx, :);
        [GC, ~, ~, ~] = compute_static_GC(Tseg, fs, order, chIdx, freq_range, false);
        GC_all(:, :, k) = GC;
    
        if nCh == 2
            GC_AB(k) = GC(2,1); % A->B
            GC_BA(k) = GC(1,2); % B->A
        end
    
        fprintf('Segment %d/%d done\n', k, nSeg);
    end
    
    GC_mean = mean(GC_all, 3, 'omitnan');
    
    % --- plotting
    figure('Name', 'GC Segments');
    
    % (1) Segment-by-segment GC time plot if only 2 channels
    if nCh == 2
        subplot(2, 1, 1);
        plot(seg_times, GC_AB, '-o', seg_times, GC_BA, '-s', 'LineWidth', 1.5);
        xlabel('Segment start (s)');
        ylabel('GC value');
        legend('A→B', 'B→A', 'Location', 'best');
        title(sprintf('GC per %.1f s segment', seg_len_sec));
        grid on;
    end
    
    % (2) Mean GC heatmap
    subplot(2, 1, 2);
    imagesc(GC_mean);
    axis square;
    colormap(parula);
    colorbar;
    
    % --- try to extract channel names
    if istimetable(TT) || istable(TT)
        allNames = TT.Properties.VariableNames;
        chanLabels = allNames(chIdx);
    else
        chanLabels = compose('Ch%d', chIdx);
    end
    
    % --- apply channel name ticks
    set(gca, 'XTick', 1:nCh, 'XTickLabel', chanLabels, ...
             'YTick', 1:nCh, 'YTickLabel', chanLabels, ...
             'TickLabelInterpreter', 'none', 'FontSize', 10);
    
    xtickangle(45);
    xlabel('Source channel');
    ylabel('Target channel');
    title(sprintf('Mean GC matrix across %d segments (%.1f s each)', nSeg, seg_len_sec));
    
    fprintf('All segments done. Mean GC matrix plotted.\n');
    
    end
    
    
    
    
    function load_MNI_FRV_mesh(mni_path, pial_lh_path, data, colorRange, pial_rh_path)
    % load_MNI_FRV_mesh(mni_path, pial_lh_path, data, colorRange, pial_rh_path)
    %
    % Mirrors load_ACPC_FR_mesh and load_FSAVG_FRS_mesh, but for MNI-space FRV
    % electrode coordinates.
    %
    % Inputs:
    %   mni_path      - path to *_elec_mni_frv.mat
    %   pial_lh_path  - path to MNI LH pial surface (surface_pial_left.mat)
    %   data          - vector of values (one per electrode)
    %   colorRange    - [min max] range for normalization (optional)
    %   pial_rh_path  - OPTIONAL: path to RH pial surface
    %
    
        if nargin < 5
            pial_rh_path = [];
        end
    
        if nargin < 4 || isempty(colorRange)
            colorRange = [min(data(:)), max(data(:))];
        end
    
    
     
        S = load(mni_path);
        if ~isfield(S, 'elec_mni_frv')
            error('File %s does not contain elec_mni_frv', mni_path);
        end
    
        elec_mni = S.elec_mni_frv;
    
        % Ensure units
        elec_mni = ft_convert_units(elec_mni, 'mm');
        elec_mni.coordsys = 'mni';
    
    
    
        pial_lh = load(pial_lh_path);
        if isfield(pial_lh, 'mesh')
            pial_lh = pial_lh.mesh;
        end
        pial_lh = ft_convert_units(pial_lh, 'mm');
        pial_lh.coordsys = 'mni';
    
    
     
        nElec = size(elec_mni.chanpos, 1);
        data = data(:);
    
        if numel(data) ~= nElec
            error('Data length (%d) does not match MNI electrodes (%d).',...
                  numel(data), nElec);
        end
    
    
        ref_idx = isnan(data);
        data_clean = data;
        data_clean(ref_idx) = 0;
    
    
        
        norm_data = (data_clean - colorRange(1)) ./ (colorRange(2) - colorRange(1));
        norm_data = max(min(norm_data, 1), 0);
    
    
    
        source = struct();
        source.pos = elec_mni.chanpos;
        source.pow = norm_data;
        source.dimord = 'pos';
        source.coordsys = 'mni';
    
    
        cfg_interp = [];
        cfg_interp.parameter    = 'pow';
        cfg_interp.interpmethod = 'sphere_weighteddistance';
        cfg_interp.sphereradius = 3;
    
        interp_lh = ft_sourceinterpolate(cfg_interp, source, pial_lh);
    
    
    
        cfg_plot = [];
        cfg_plot.method       = 'surface';
        cfg_plot.funparameter = 'pow';
        cfg_plot.funcolorlim  = [0 1];
        cfg_plot.interpmethod = 'nearest';
        cfg_plot.sphereradius = 3;
        cfg_plot.camlight     = 'no';
    
    
    
        ft_sourceplot(cfg_plot, interp_lh);
        hold on;
    
    
    
        if ~isempty(pial_rh_path) && exist(pial_rh_path, 'file')
            pial_rh = load(pial_rh_path);
            if isfield(pial_rh, 'mesh')
                pial_rh = pial_rh.mesh;
            end
            pial_rh = ft_convert_units(pial_rh, 'mm');
            pial_rh.coordsys = 'mni';
    
            interp_rh = ft_sourceinterpolate(cfg_interp, source, pial_rh);
            ft_sourceplot(cfg_plot, interp_rh);
        end
    
    
        view([-55 10]);
        material dull;
        lighting gouraud;
        camlight;
        colormap(parula);
    
    
    
        ft_plot_sens(elec_mni, 'elecsize', 35,...
            'facecolor', norm_data, 'edgecolor', 'black');
    
    
    
        if any(ref_idx)
            scatter3(elec_mni.chanpos(ref_idx,1),...
                     elec_mni.chanpos(ref_idx,2),...
                     elec_mni.chanpos(ref_idx,3),...
                     100, 'r', 'filled', 'MarkerEdgeColor','k','LineWidth',1.5);
        end
    
    
    
        cb = colorbar;
        ticks = linspace(0, 1, 5);
        cb.Ticks = ticks;
        cb.TickLabels = string(linspace(colorRange(1), colorRange(2), 5));
        ylabel(cb, 'PLV value');
    
    end
