%% CORTICAL-STN CONNECTIVITY | faithful, organized reproduction
% MATLAB R2025a. Requires original project helpers/toolboxes on MATLAB path:
% myRetime, openNSx, ns2timetable, Myeegfilt, instPhaseFreq,
% MVGC tsdata_to_var, optional FieldTrip.
% Preserves ORIGINAL alignment, epoch selection, PLV filtering and GC fit.
% Only intended scientific departures: remove unconditional wrong-patient
% overwrite; correct uninitialized meta/ACPC_Path and plotting errors.
% This has not been executed against the original patient recordings.
clearvars; close all;
set(groot,'DefaultFigureVisible','on');

%% USER CONFIGURATION (change here only)
Patient_Folder_Path = 'C:\Users\cours\Downloads\AndersonLabDataAndPackages\Patient Data';
whichPatient = 'PD25N006';
FieldTrip_Path = 'C:\Users\cours\Downloads\AndersonLabDataAndPackages\Packages\fieldtrip-20251023';
showFieldTrip = true;
showBrainnetomeGrid = true;
runGC = true;
runAIC = true;
runADF = true;
showNOEvents = true;
runPLVSurrogates = true;
nSur = 100; seed = 12345;
channels = [19 20 21 40 41 42 61 62 63]; % original reporting subset
chIdx = [21 70];  % ORIGINAL GC positional indices (not electrode numbers)
segment_length_sec = 0.25; order = 25; filtOrder = 1024;
ECoG_ch = 'LS33';

%% PATIENT PRESETS (original exact values)


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
        "backCut",     3800, ... %3820 or 3800
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


P = presets.(whichPatient);
Patient_ID = char(P.Patient_ID);
NO_timetable_name = P.NO_timetable_name;
NO_SamplingFreq = P.NO_SamplingFreq;
ns_files = P.ns_files;
use_ns2_only = P.use_ns2_only;
NO_to_BR_ch = P.NO_to_BR_ch;
BR_to_NO_ch = P.BR_to_NO_ch;
BR_SamplingFreq = 1000;
fs = BR_SamplingFreq;

%% 1. LOAD RECORDINGS (original algorithm)


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
    if strcmp(whichPatient,'PD24N009')
        ns_timetables = {}; % ns2/ns5 have different channel sets and are joined below
    else
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
    end
    


% Unlike the original, load the ns5/ns2 combination ONLY for PD24N009.
if strcmp(whichPatient,'PD24N009')


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


end


%% 2. BUILD RESAMPLED DATA (original algorithm)


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


%% 3. PRE-ALIGNMENT DIAGNOSTICS


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


%% 4. FIRST CROSS-CORRELATION (unchanged helper)


    fs        = 1000;
    br_win    = P.br_win;
    no_win    = P.no_win;
    BR_to_NO_ch       = P.BR_to_NO_ch;
    %   [600 850], ...
       % "no_win",       [2900 3150], ....
    
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


%% 5. BUILD INITIAL SYNCHRONIZED TIMETABLE


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


%% 6. SELECT EPOCH, THEN REFINE ALIGNMENT


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


%% 7. PLOT FINAL SYNCHRONIZED SIGNALS


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


% Remove file label only if it exists: otherwise phase extraction
% attempts to bandpass-filter the nonnumeric SourceFile strings.
if ismember('SourceFile', TT_longest_refined.Properties.VariableNames)
    TT_longest_refined = removevars(TT_longest_refined,'SourceFile');
end
assert(~isempty(TT_longest_refined),'Selected epoch is empty.');
fprintf('Final epoch: %s to %s | samples: %d\n', ...
    string(TT_longest_refined.Time(1)),string(TT_longest_refined.Time(end)),height(TT_longest_refined));
fprintf('First lag %.6f s | refined lag %.6f s\n',lag_s,refined_lag_s);

%% 7B. OPTIONAL ORIGINAL ADF DIAGNOSTIC
if runADF
    if exist('adftest','file') == 2 && all(ismember({ECoG_ch,NO_to_BR_ch},TT_longest_refined.Properties.VariableNames))
        N=min(1000,height(TT_longest_refined));
        xe=TT_longest_refined.(ECoG_ch); xl=TT_longest_refined.(NO_to_BR_ch);
        x_ecog_w=xe(1:N);x_lfp_w=xl(1:N);
        [he,pe,se]=adftest(x_ecog_w);
        [hl,pl,sl]=adftest(x_lfp_w);
        fprintf('ADF first %.3f s: ECoG h=%d p=%.5f stat=%.4f | STN h=%d p=%.5f stat=%.4f\n', ...
            N/fs,he,pe,se,hl,pl,sl);
        figure('Name','ADF diagnostic signals');
        subplot(2,1,1);plot(TT_longest_refined.Time(1:N),x_ecog_w);title(ECoG_ch);
        subplot(2,1,2);plot(TT_longest_refined.Time(1:N),x_lfp_w);title(NO_to_BR_ch);
    else
        warning('ADF diagnostic skipped (missing adftest or signal column).');
    end
end

%% 7C. OPTIONAL ORIGINAL NO EVENT/WINDOW VIEW
if showNOEvents
    NO_to_BR_offset=NO_on_BR_refined.Time(1)-NO_timetable.Time(1);
    tNO_start=TT_longest_refined.Time(1)-NO_to_BR_offset;
    tNO_end=TT_longest_refined.Time(end)-NO_to_BR_offset;
    figure('Name','Final selected epoch in native NO time');
    plot(NO_timetable.Time,NO_timetable.(NO_to_BR_ch),'k');hold on;
    xline(tNO_start,'r','LineWidth',2);xline(tNO_end,'r','LineWidth',2);
    if ~isempty(NO_timetable.Properties.Events)
        E=NO_timetable.Properties.Events;
        if ~isempty(E),xline(E.Time,'--b');end
    end
    title('Final analysis epoch on Neural Omega timebase');grid on;
    fprintf('NO native epoch from %s to %s\n',string(tNO_start),string(tNO_end));
end

%% 8. PHASE AND PLV (original exact helpers)


    assert(all(varfun(@isnumeric,TT_longest_refined,'OutputFormat','uniform')), ...
    'Synchronized timetable contains nonnumeric variables; inspect its columns.');
[Phase_Data, Channel_Names] = obtainPhase_tt(TT_longest_refined, 1000, 'beta');
    disp(Channel_Names)
    disp('Obtained Phase Data')
    PLV_VectorMatrix = compute_PLV_VectorMatrix(Phase_Data); %Get the PLV vector matrix for every channel (64x64)
    PLV_Average = getAveragePLV(PLV_VectorMatrix); %Get the Average PLV for every channel (64x1)
    disp('Calculated PLV Vector Matrix')
    
    
    %% Visualize PLV Vector Matrix
    visualize_PLV_Matrix(PLV_VectorMatrix, Channel_Names)
    


ref_channel = find(strcmp(Channel_Names, NO_to_BR_ch),1);
assert(~isempty(ref_channel),'Reference %s missing from Phase_Data.',NO_to_BR_ch);
PLV_Relative_ref = getPLV_wrt_Ref(PLV_VectorMatrix, ref_channel);
assert(numel(PLV_Relative_ref)>=63,'Expected 63 ECoG channels.');
plot_PLV_Relative(Phase_Data(1:63),Phase_Data{ref_channel}, ...
    PLV_Relative_ref,ref_channel,Channel_Names);

%% 9. ELECTRODE METADATA + BRAINNETOME 21x3 PLV GRID
meta = table();
if showBrainnetomeGrid
    electrode_data_filepath = findElectrodeMetadata(fullfile(Patient_Folder_Path,Patient_ID),P.electrode_data_name);
    meta = readtable(electrode_data_filepath,'VariableNamingRule','preserve');
    plot_PLV_Relative_Brainnetome(Phase_Data(1:63),Phase_Data{ref_channel}, ...
        PLV_Relative_ref,ref_channel,Channel_Names,meta,P.badChannel);
end

%% 10. ORIGINAL PLV SURROGATE TEST AND SUMMARY
if runPLVSurrogates
    rng(seed,'twister');
    PLV = getPLV_wrt_Ref(PLV_VectorMatrix,ref_channel);
    PLV_raw = PLV(channels);
    PLV_norm = PLV_raw/max(PLV_raw);
    PLV_surrogate = nan(numel(channels),nSur);
    ref_phase = Phase_Data{ref_channel};
    T = numel(ref_phase);
    for s=1:nSur
        shift=randi(T-1);
        phase_surr=Phase_Data;
        phase_surr{ref_channel}=circshift(ref_phase,shift);
        vm=compute_PLV_VectorMatrix(phase_surr);
        v=getPLV_wrt_Ref(vm,ref_channel);
        PLV_surrogate(:,s)=v(channels);
        if mod(s,10)==0 || s==1 || s==nSur
            fprintf('Surrogate %d/%d\n',s,nSur);
        end
    end
    mu_surr=mean(PLV_surrogate,2);
    p_emp_all=(sum(PLV_surrogate>=PLV_raw,2)+1)/(nSur+1);
    [~,ch_max]=max(PLV_raw);
    fprintf('Highest PLV in reporting subset: %s, raw %.4f, normalized %.3f, surrogate p %.5f\n', ...
        Channel_Names{channels(ch_max)},PLV_raw(ch_max),PLV_norm(ch_max),p_emp_all(ch_max));
    % Region summaries correspond to the original reporting subset, not all 63.
    if ~isempty(meta)
        region_labels=brainnetomeRegions(meta,Channel_Names(channels));
        regions=unique(region_labels);
        means=nan(numel(regions),1);
        for r=1:numel(regions)
            means(r)=mean(PLV_raw(region_labels==regions(r)),'omitnan');
        end
        [~,ridx]=max(means); topRegion=regions(ridx);
        idx_top=region_labels==topRegion;
        fprintf('Highest mean PLV region in selected subset: %s, mean %.4f\n',topRegion,means(ridx));
        fprintf('Channels in top region above their surrogate means: %.1f%%\n', ...
            100*mean(PLV_raw(idx_top)>mu_surr(idx_top)));
    end
end

%% 11. FIELDTRIP MNI PLV SURFACE (original mesh function)
if showFieldTrip
    if ~isfolder(FieldTrip_Path)
        warning('FieldTrip not found: %s. Skipping surface.',FieldTrip_Path);
    else
        addpath(FieldTrip_Path); ft_defaults;
        [~,ftpath]=ft_version;
        MNI_FRV_Path=fullfile(Patient_Folder_Path,Patient_ID,'Imaging',[Patient_ID '_elec_mni_frv.mat']);
        rh_pial_mni_path=fullfile(ftpath,'template','anatomy','surface_pial_right.mat');
        if exist(MNI_FRV_Path,'file') && exist(rh_pial_mni_path,'file')
            try
                load_MNI_FRV_mesh(MNI_FRV_Path,rh_pial_mni_path,PLV_Relative_ref(1:63),[]);
            catch ME
                warning('FieldTrip surface failed: %s',ME.message);
            end
        else
            warning('MNI electrodes or RH pial template missing; skipping FieldTrip figure.');
        end
    end
end

%% 12. MVAR ORDER DIAGNOSTIC (original non-filtered AIC)
if runGC && runAIC


    maxOrder  = 50;
    fs        = 1000;
    blockSize = 250;
    TT_all    = TT_longest_refined;
    chIdx     = [21 70]; %61 70 %63 74 %33 69 %21 70
    
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
    


end


%% 13. BIVARIATE BETA GC (original MVGC estimator)
if runGC


    % GC settings are defined in USER CONFIGURATION above.

    fs = 1000;
    
    A_name = Channel_Names{chIdx(1)};
    B_name = Channel_Names{chIdx(2)};
    TT = TT_longest_refined;
    fprintf('GC columns [%d %d] = {%s, %s}\n',chIdx(1),chIdx(2), ...
        TT.Properties.VariableNames{chIdx(1)},TT.Properties.VariableNames{chIdx(2)});
    
    % --- Prefilter the full-length signals once ---
    filtOrder = 1024;
    
    A_full = double(TT{:, chIdx(1)});  A_full = A_full(:).';
    B_full = double(TT{:, chIdx(2)});  B_full = B_full(:).';
    
    A_beta = Myeegfilt(A_full, fs, 13, 30, 0, filtOrder);
    B_beta = Myeegfilt(B_full, fs, 13, 30, 0, filtOrder);
    
    % --- Segmenting ---
    samples_per_segment = round(segment_length_sec * fs);
    nTotalSamples = height(TT);
    nSegments = floor(nTotalSamples / samples_per_segment);
    
    fprintf('Splitting into %d segments of %.2f s each (%d samples)\n', ...
            nSegments, segment_length_sec, samples_per_segment);
    
    GC_all    = nan(2,2,nSegments);
    GC_AB_seg = nan(nSegments,1);
    GC_BA_seg = nan(nSegments,1);
    
    ok_seg    = false(nSegments,1);
    fail_msgs = strings(nSegments,1);
    
    for k = 1:nSegments
        idxStart = (k-1)*samples_per_segment + 1;
        idxEnd   = k*samples_per_segment;
    
        % MVGC expects variables x time
        Xseg = [A_beta(idxStart:idxEnd);
                B_beta(idxStart:idxEnd)];
    
        [GC, ok, msg] = mvgc_uncond_time_gc_2ch(Xseg, order);
    
        if ~ok
            fail_msgs(k) = string(msg);
            % uncomment if you want to see each failure:
            % warning('Segment %d/%d failed: %s', k, nSegments, msg);
            continue;
        end
    
        ok_seg(k) = true;
    
        GC_all(:,:,k) = GC;
        GC_AB_seg(k)  = GC(2,1);    % A -> B
        GC_BA_seg(k)  = GC(1,2);    % B -> A
    
        fprintf('Segment %2d/%2d done | GC(%s→%s)=%.4f | GC(%s→%s)=%.4f\n', ...
            k, nSegments, A_name, B_name, GC_AB_seg(k), B_name, A_name, GC_BA_seg(k));
    end
    
    fprintf('\nMVGC: %d/%d segments succeeded.\n', nnz(ok_seg), nSegments);
    
    % show top failure reasons if lots failed
    if nnz(ok_seg) < nSegments
        msgs = fail_msgs(fail_msgs ~= "");
        if ~isempty(msgs)
            [u,~,ic] = unique(msgs);
            counts = accumarray(ic,1);
            [counts,idx] = sort(counts,'descend');
            u = u(idx);
            fprintf('Top failure reasons:\n');
            for i = 1:min(5,numel(u))
                fprintf('  %dx  %s\n', counts(i), u(i));
            end
        end
    end
    
    %Calculate mean GC
    GC_mean = mean(GC_all, 3, 'omitnan');
    
    % === Statistics  ===
    valid = isfinite(GC_AB_seg) & isfinite(GC_BA_seg);
    fprintf('Valid paired segments for stats: %d/%d\n', nnz(valid), nSegments);
    
    if nnz(valid) >= 3
        [p_val,~,stats] = signrank(GC_AB_seg(valid), GC_BA_seg(valid));
        if isfield(stats,'zval'), fprintf('\nWilcoxon signed-rank: z = %.3f, p = %.4g\n', stats.zval, p_val);
        else, fprintf('\nWilcoxon signed-rank: p = %.4g\n',p_val); end
    
        [~,p_ttest,~,stats_t] = ttest(GC_AB_seg(valid), GC_BA_seg(valid));
        fprintf('Paired t-test: t(%d) = %.3f, p = %.4g\n', stats_t.df, stats_t.tstat, p_ttest);
    
    else
        warning('Not enough valid segments for signrank/ttest. Skipping stats.');
    end


% Use paired valid segments for means, uncertainty, and the test.
    nValid=nnz(valid);
    p_val=NaN; zval=NaN; r_effect=NaN;
    if nValid>=3
        [p_val,~,stats]=signrank(GC_AB_seg(valid),GC_BA_seg(valid));
        if isfield(stats,'zval'), zval=stats.zval; end
        if isfinite(zval),r_effect=abs(zval)/sqrt(nValid);end
    end
    fprintf('Directional effect size r = %.4f (n=%d paired segments)\n',r_effect,nValid);


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


%% NEW VISUALIZE (MGVC Version)
    % ===== Labels =====
    A_label = A_name;
    B_label = B_name;

    % ===== Color specification =====
    col_AB = [0.15 0.55 0.87];   % blue (A → B)
    col_BA = [0.87 0.33 0];   % red  (B → A)

    
    avg_sec = 4;  % desired averaging window (seconds)
    
    % ===== Valid mask (prevents NaNs from breaking plots/stats) =====
    valid = isfinite(GC_AB_seg) & isfinite(GC_BA_seg);
    
    % Segment start times
    t_seg = (0:nSegments-1) * segment_length_sec;
    
    % If you haven't computed stats yet, do it here (robust)
    if nnz(valid) >= 3
        [p_val,~,stats] = signrank(GC_AB_seg(valid), GC_BA_seg(valid));
        if isfield(stats,'zval'), zval = stats.zval; else, zval=NaN; end
    else
        p_val = NaN; zval = NaN;
    end
    
    figure('Name','Segmented GC Summary','Position',[100 100 1100 900]);
    
    % (1) Segment-wise GC vs time (RAW)
    subplot(3,1,1);
    plot(t_seg, GC_AB_seg, 'Color', col_AB, 'LineWidth',1, ...
        'DisplayName', sprintf('%s→%s', A_label, B_label)); hold on;
    plot(t_seg, GC_BA_seg, 'Color', col_BA, 'LineWidth',1, ...
        'DisplayName', sprintf('%s→%s', B_label, A_label));
    xlabel('Segment start (s)');
    ylabel('GC value');
    legend('Location','best'); grid on;
    title(sprintf('Granger Causality over %.3g s windows', ...
        segment_length_sec));
    
    % (1b) GC averaged over avg_sec (NON-overlapping blocks)
    subplot(3,1,2);
    
    nPerBlock = max(1, round(avg_sec / segment_length_sec));
    nBlocks   = floor(nSegments / nPerBlock);
    
    GC_AB_blk = nan(nBlocks,1);
    GC_BA_blk = nan(nBlocks,1);
    t_blk     = nan(nBlocks,1);
    
    for b = 1:nBlocks
        idx = (b-1)*nPerBlock + (1:nPerBlock);
        GC_AB_blk(b) = mean(GC_AB_seg(idx), 'omitnan');
        GC_BA_blk(b) = mean(GC_BA_seg(idx), 'omitnan');
        t_blk(b)     = mean(t_seg(idx));
    end
    
    plot(t_blk, GC_AB_blk, 'Color', col_AB, 'LineWidth',2, ...
        'DisplayName', sprintf('%s→%s', A_label, B_label)); hold on;
    plot(t_blk, GC_BA_blk, 'Color', col_BA, 'LineWidth',2, ...
        'DisplayName', sprintf('%s→%s', B_label, A_label));
    xlabel('Time (s)');
    ylabel('GC value');
    legend('Location','best'); grid on;
    title(sprintf('GC Averaged Over %.1f s', avg_sec));
    

    
    % (3a) Bar comparison (means)
    subplot(3,2,5);
    
    mean_AB = mean(GC_AB_seg, 'omitnan');
    mean_BA = mean(GC_BA_seg, 'omitnan');
    barVals = [mean_AB, mean_BA];
    
    bb = bar(barVals, 0.5);
    
    hold on;
    errorbar(1, barVals(1), sem_AB, 'Color',[0 0 0.5], ...
    'LineStyle','none','LineWidth',1,'CapSize',12);   % bright cyan-blue
    errorbar(2, barVals(2), sem_BA, 'Color',[0.5 0 0], ...
    'LineStyle','none','LineWidth',1,'CapSize',12);   % bright salmon-red


    bb.FaceColor = 'flat';
    bb.CData = [col_AB; col_BA];
    
    set(gca,'XTickLabel',{sprintf('%s→%s',A_label,B_label), ...
                          sprintf('%s→%s',B_label,A_label)});
    ylabel('Mean GC');
    
    if isfinite(p_val)
        title(sprintf('Wilcoxon signed-rank: p=%.4g (z=%.3f)', p_val, zval));
    else
        title('Wilcoxon signed-rank: not enough valid segments');
    end
    
    grid on;
    
    mbar = max(barVals);
    if isfinite(mbar) && mbar > 0
        ylim([mbar*0.75, mbar*1.1]);
    end
    
    xt = 1:2;
    for i = 1:numel(barVals)
        if isfinite(barVals(i))
            text(xt(i), barVals(i) + 0.04*mbar, sprintf('%.4f', barVals(i)), ...
                'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
        end
    end
    
    % Mean 2×2 GC matrix (true average)
    subplot(3,2,6);
        imagesc(GC_mean);
    axis square; colormap(parula); colorbar;
    
    mx = max(GC_mean(:), [], 'omitnan');
    if isfinite(mx) && mx > 0
        caxis([0 mx*1.1]);
    end
    
    xlabel('Source'); ylabel('Target');
    set(gca,'XTick',1:2,'XTickLabel',{A_label, B_label}, ...
            'YTick',1:2,'YTickLabel',{A_label, B_label}, 'FontSize',8);
    title('Mean GC Matrix');


end


%% 14. SAVE RESULTS (without altering calculated quantities)
resultsFile=fullfile(Patient_Folder_Path,Patient_ID,[Patient_ID '_Cortical_STN_reproduced.mat']);
if runGC
    save(resultsFile,'whichPatient','P','lag_s','refined_lag_s', ...
        'TT_longest_refined','Channel_Names','PLV_Relative_ref','PLV_VectorMatrix', ...
        'GC_mean','GC_AB_seg','GC_BA_seg','chIdx','order','segment_length_sec','-v7.3');
else
    save(resultsFile,'whichPatient','P','lag_s','refined_lag_s', ...
        'TT_longest_refined','Channel_Names','PLV_Relative_ref','PLV_VectorMatrix','-v7.3');
end
fprintf('Saved %s\n',resultsFile);

%% LOCAL FUNCTIONS (preserved from original where scientific output matters)


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
        PLV_copy(64:end) = -1;
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
    
            for i = 1:min(numChannels,63)
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
            ylabel(cbar, sprintf('PLV with respect to %s',ref_name),'Interpreter','none');
        end
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
    
        % --- choose electrode positions robustly ---
        if isfield(elec_mni,'elecpos') && ~isempty(elec_mni.elecpos)
            pos = elec_mni.elecpos;
        else
            pos = elec_mni.chanpos;
        end

    
        pial_lh = load(pial_lh_path);
        if isfield(pial_lh, 'mesh')
            pial_lh = pial_lh.mesh;
        end
        pial_lh = ft_convert_units(pial_lh, 'mm');
        pial_lh.coordsys = 'mni';
    
    
     
        nElec = size(pos, 1);
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
        source.pos = pos;
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
            scatter3(pos(ref_idx,1),...
                     pos(ref_idx,2),...
                     pos(ref_idx,3),...
                     100, 'r', 'filled', 'MarkerEdgeColor','k','LineWidth',1.5);
        end
    
    
    
        cb = colorbar;
        ticks = linspace(0, 1, 5);
        cb.Ticks = ticks;
        cb.TickLabels = string(linspace(colorRange(1), colorRange(2), 5));
        ylabel(cb, 'PLV value');
    
    end







function [F, ok, msg] = mvgc_uncond_time_gc_2ch(X, p)
% Unconditional (pairwise) time-domain GC via MVGC VAR fitting.
% X: 2 x T, p: VAR order. Returns F(i,j)=j->i

    % Error Checks
    ok = false; msg = '';
    F = nan(2);
    
    if size(X,1) ~= 2
        msg = 'X must be 2xT.';
        return;
    end
    if any(~isfinite(X(:)))
        msg = 'Non-finite values in segment.';
        return;
    end
    T = size(X,2);
    if T <= p + 2
        msg = sprintf('Too few samples (T=%d) for p=%d.', T, p);
        return;
    end

    % Main Calculation %

    %Preprocessing - Demean and z-score
    X = X - mean(X,2);
    s = std(X,0,2);
    if any(s < 1e-12)
        msg = 'Near-zero std in a channel.';
        return;
    end
    X = X ./ s;

    % Fit full MVAR model (using MVGC function)
    try
        [~, SIGfull] = tsdata_to_var(X, p, 'OLS');      % MVGC time-series data to VAR model
    catch ME
        msg = ['tsdata_to_var(full) failed: ' ME.message];
        return;
    end
    vfull = diag(SIGfull);

    % Fit reduced MVAR model (using MVGC function)
    try
        [~, SIG1] = tsdata_to_var(X(1,:), p, 'OLS');    % channel 1 alone
        [~, SIG2] = tsdata_to_var(X(2,:), p, 'OLS');    % channel 2 alone
    catch ME
        msg = ['tsdata_to_var(reduced) failed: ' ME.message];
        return;
    end
    
    % Calculate variance of residual errors
    v1_red = SIG1(1);  % variance of residual for ch1-only model
    v2_red = SIG2(1);  % variance of residual for ch2-only model

    % GC = log(var_reduced / var_full)
    F = zeros(2);
    F(1,2) = log(v1_red / vfull(1));  % 2 -> 1  (B -> A)
    F(2,1) = log(v2_red / vfull(2));  % 1 -> 2  (A -> B)

    if any(~isfinite(F(:)))
        msg = 'Non-finite GC values.';
        F(:) = nan;
        return;
    end

    ok = true;
end




function path = findElectrodeMetadata(folder,name)
    candidates={fullfile(folder,char(name)),fullfile(folder,[char(name) '.csv']),fullfile(folder,[char(name) '.xlsx'])};
    path='';
    for k=1:numel(candidates)
        if isfile(candidates{k}),path=candidates{k};return;end
    end
    error('Electrode metadata file not found for %s in %s',char(name),folder);
end

function regions = brainnetomeRegions(meta, names)
    assert(all(ismember({'Electrode','Brainnetome'},meta.Properties.VariableNames)), ...
        'Metadata needs Electrode and Brainnetome columns.');
    keys=regexp(string(meta.Electrode),'\d+','match','once');
    cleanLabels=regexprep(string(meta.Brainnetome),'\s*\(\d+%?\)\s*$','');
    m=containers.Map('KeyType','char','ValueType','char');
    for k=1:numel(keys)
        if ~ismissing(keys(k)) && strlength(keys(k))>0
            num=str2double(keys(k));
            if isfinite(num),m(sprintf('%d',num))=char(cleanLabels(k));end
        end
    end
    regions=strings(numel(names),1);
    for k=1:numel(names)
        token=regexp(string(names{k}),'\d+','match','once');
        num=str2double(token);
        key=sprintf('%d',num);
        if isfinite(num)&&isKey(m,key),regions(k)=string(m(key));
        else,regions(k)="Unlabeled";end
    end
end

function plot_PLV_Relative_Brainnetome(~,~,PLV_relative_vector,~,channelNames,meta,badCh)
    % Same 21x3 column-major electrode arrangement as the old visualizer.
    % Cell fill=PLV (raw or normalized); cell edge=Brainnetome region.
    n=63; values=PLV_relative_vector(1:n);
    values(badCh(badCh>=1 & badCh<=n))=NaN;
    regions=brainnetomeRegions(meta,channelNames(1:n));
    valid=isfinite(values);
    if ~any(valid),warning('No finite cortical PLVs to plot.');return;end
    vmax=max(values(valid));
    normVals=values/vmax;
    uniqueRegions=unique(regions);
    palette=lines(numel(uniqueRegions));
    f=figure('Name','PLV 21x3 Brainnetome grids','Position',[100 100 1300 850]);
    for panel=1:2
        ax=subplot(1,2,panel,'Parent',f);hold(ax,'on');
        for i=1:n
            x=floor((i-1)/21);y=mod(i-1,21);
            r=find(uniqueRegions==regions(i),1);
            v=values(i);if panel==2,v=normVals(i);end
            if isfinite(v)
                j=max(1,min(256,round(v*255)+1));
                cm=parula(256);face=cm(j,:);
            else
                face=[.75 .75 .75];
            end
            rectangle(ax,'Position',[x y 1 1],'FaceColor',face, ...
                'EdgeColor',palette(r,:),'LineWidth',2);
            if isfinite(v),txt=sprintf('%.3f',v);else,txt='NaN';end
            text(ax,x+.5,y+.5,txt,'HorizontalAlignment','center','FontSize',7);
        end
        axis(ax,'equal');xlim(ax,[0 3]);ylim(ax,[0 21]);
        set(ax,'XTick',[],'YTick',[],'YDir','normal');
        colormap(ax,parula(256));caxis(ax,[0 1]);colorbar(ax);
        if panel==1,title(ax,'PLV wrt STN (raw)');else,title(ax,'PLV wrt STN (normalized)');end
    end
    % Proxy legend uses colored region-border samples (no data modification).
    ax=subplot(1,2,2);hold(ax,'on');h=gobjects(numel(uniqueRegions),1);
    for r=1:numel(uniqueRegions)
        h(r)=plot(ax,nan,nan,'s','MarkerEdgeColor',palette(r,:),'LineWidth',2);
    end
    legend(ax,h,cellstr(uniqueRegions),'Location','eastoutside','Interpreter','none');
    T=table(regions,values,'VariableNames',{'Region','PLV'});
    stats=groupsummary(T,'Region','mean','PLV');
    disp('Brainnetome PLV summary:');disp(sortrows(stats,'mean_PLV','descend'));
end

