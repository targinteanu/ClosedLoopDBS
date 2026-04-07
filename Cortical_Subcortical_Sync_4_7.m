%% ========================================================================
%  Patient Presets
%  ========================================================================

presets = struct();

presets.PD25N006 = struct( ...
    "Patient_ID", "PD25N006", ...
    "NO_timetable_name", "SavedTable_Distances_1375Hz1.mat", ...
    "NO_SamplingFreq", 1375, ...
    "NO_to_BR_ch", "CMacro_LFP_01", ...
    "BR_to_NO_ch", "AINP3", ...
    "ns_files", {{'datafile001.ns2'}}, ...
    "use_ns2_only", true, ...
    "target_window", [-5 3], ...
    "target_channel", "CMacro_LFP_01", ...
    "minDurSec", 20, ...
    "br_win", [1800 2500], ...
    "no_win", [1700 2400], ...
    "frontCut", 60, ...
    "backCut", 3800, ...
    "electrode_data_name", "electrode_data(Sheet1).csv", ...
    "badChannel", [] ...
);

presets.PD25N008 = struct( ...
    "Patient_ID", "PD25N008", ...
    "NO_timetable_name", "SavedTable_Distances_1375Hz1.mat", ...
    "NO_SamplingFreq", 1375, ...
    "NO_to_BR_ch", "CMacro_LFP_AP_T", ...
    "BR_to_NO_ch", "AINP2", ...
    "ns_files", {{'baseline001.ns2', 'DepthStim002.ns2', 'nback004.ns2'}}, ...
    "use_ns2_only", true, ...
    "target_window", [-5 3], ...
    "target_channel", "CMacro_LFP_AP_T", ...
    "minDurSec", 20, ...
    "br_win", [800 1500], ...
    "no_win", [400 1100], ...
    "frontCut", 20, ...
    "backCut", 20, ...
    "electrode_data_name", "PD25N008_electrode_data", ...
    "badChannel", [48] ...
);

presets.PD24N009 = struct( ...
    "Patient_ID", "PD24N009", ...
    "NO_timetable_name", "SavedTable_Distances_1375Hz1.mat", ...
    "NO_SamplingFreq", 1375, ...
    "NO_to_BR_ch", "CMacro_LFP_01", ...
    "BR_to_NO_ch", "AINP4_ns5", ...
    "ns_files", {{'datafile001.ns2', 'datafile001.ns5'}}, ...
    "use_ns2_only", false, ...
    "target_window", [-5 3], ...
    "target_channel", "CMacro_LFP_01", ...
    "minDurSec", 20, ...
    "br_win", [600 850], ...
    "no_win", [2900 3150], ...
    "frontCut", 600, ...
    "backCut", 440, ...
    "electrode_data_name", "PD24N009_electrode_data", ...
    "badChannel", [] ...
);

%% ========================================================================
%  Patient Selection
%  ========================================================================

patientFolder = 'C:\Users\cours\Downloads\AndersonLabDataAndPackages\Patient Data';
whichPatient  = 'PD25N006';

P = presets.(whichPatient);

patientID         = char(P.Patient_ID);
noTimetableName   = P.NO_timetable_name;
noFs              = P.NO_SamplingFreq;
brFs              = 1000;

nsFiles           = P.ns_files;
useNs2Only        = P.use_ns2_only;

noToBrChannel     = P.NO_to_BR_ch;
brToNoChannel     = P.BR_to_NO_ch;

targetChannel     = P.target_channel;
targetWindow      = P.target_window;
minDurSec         = P.minDurSec;
frontCut          = P.frontCut;
backCut           = P.backCut;

%% ========================================================================
%  Load NeuroOmega (NO) Data

noTimetablePath = fullfile(patientFolder, patientID, noTimetableName);
noDataRaw       = load(noTimetablePath);

noTimetableRaw  = sortrows(noDataRaw.Tbl, 'Time');
noTimetable     = myRetime(noTimetableRaw, noFs, 0, 2);

noTargetDistRaw = noTimetable.targetDistance;
noSignalRaw     = noTimetable.(noToBrChannel);

%% ========================================================================
%  Load Blackrock (BR) Data

nsTimetables = cell(size(nsFiles));

for k = 1:numel(nsFiles)
    nsPath = fullfile(patientFolder, patientID, nsFiles{k});
    disp("Loading: " + nsPath)

    ns = openNSx(nsPath);
    fsCurrent = ns.MetaTags.SamplingFreq;

    tt = ns2timetable(nsPath);

    if fsCurrent > brFs
        disp("Downsampling " + nsFiles{k} + " from " + num2str(fsCurrent) + " Hz to 1 kHz...")
        dataDownsampled = resample(double(tt{:, :}), brFs, fsCurrent);
        tDownsampled    = tt.Time(1) + seconds((0:size(dataDownsampled,1)-1)' / brFs);

        tt = array2timetable( ...
            dataDownsampled, ...
            'RowTimes', tDownsampled, ...
            'VariableNames', tt.Properties.VariableNames ...
        );
    else
        tt.SourceFile = repmat(string(nsFiles{k}), height(tt), 1);
    end

    nsTimetables{k} = tt;
end

combinedBrTimetable = vertcat(nsTimetables{:});
combinedBrTimetable = sortrows(combinedBrTimetable, 'Time');

%% ========================================================================
%  Special Handling for PD24N009 (run this only for PD24N009)

if strcmp(patientID, 'PD24N009')
    ns2File = fullfile(patientFolder, patientID, 'datafile001.ns2');
    ns5File = fullfile(patientFolder, patientID, 'datafile001.ns5');

    T2 = sortrows(ns2timetable(ns2File), 'Time');
    T5 = sortrows(ns2timetable(ns5File), 'Time');

    overlapStart = max(T2.Time(1), T5.Time(1));
    overlapEnd   = min(T2.Time(end), T5.Time(end));

    T2 = T2(T2.Time >= overlapStart & T2.Time <= overlapEnd, :);
    T5 = T5(T5.Time >= overlapStart & T5.Time <= overlapEnd, :);

    T5resampled = retime(T5, T2.Time, 'linear');
    T5resampled = renamevars( ...
        T5resampled, ...
        T5resampled.Properties.VariableNames, ...
        append(T5resampled.Properties.VariableNames, "_ns5") ...
    );

    combinedBrTimetable = [T2 T5resampled];
end

%% ========================================================================
%  Resample BR + NO Data to 1 kHz 

brTimetable = myRetime(combinedBrTimetable, brFs, 0, 2);

tNoRaw      = (0:height(noTargetDistRaw)-1)' / noFs;
signalDur   = (height(noTargetDistRaw)-1) / noFs;
tBrCommon   = (0:round(signalDur * brFs))' / brFs;

validLen = min(numel(noTargetDistRaw), numel(tNoRaw));
noTargetDistRaw = noTargetDistRaw(1:validLen);
tNoRaw          = tNoRaw(1:validLen);

noTargetDist = interp1(tNoRaw, noTargetDistRaw, tBrCommon, 'linear');
noSignal     = resample(noSignalRaw, brFs, noFs);
brSignal     = brTimetable.(brToNoChannel);

allVars          = noTimetable.Properties.VariableNames;
signalVars       = setdiff(allVars, "targetDistance");
resampledNoCells = cellfun( ...
    @(v) resample(double(noTimetable.(v)), brFs, noFs), ...
    signalVars, ...
    'UniformOutput', false ...
);
noResampledMatrix = [resampledNoCells{:}];

if ~useNs2Only
    commonLen = min(numel(noTargetDist), size(noResampledMatrix, 1));
    noTargetDist      = noTargetDist(1:commonLen);
    noResampledMatrix = noResampledMatrix(1:commonLen, :);
end

%% ========================================================================
%  Pre-Alignment Visualization

brSignal = clipAndDemean(brSignal, 2000);
noSignal = clipAndDemean(noSignal, 2000);

tBrPlot     = (0:length(brSignal)-1)' / brFs;
tNoPlot     = (0:length(noSignal)-1)' / brFs;
tNoDistPlot = (0:length(noTargetDist)-1)' / brFs;

%% ========================================================================
%  Visualize Raw BR + NO

figure;

subplot(3,1,1);
plot(tBrPlot, brSignal, 'b'); hold on;
xlabel('Time (s)');
ylabel('Amplitude');
title('BR\_signal (clipped + demeaned, pre-alignment)');
ylim([-2000 2000]);

segmentLengths = cellfun(@height, nsTimetables);
segmentEdges   = [0, cumsum(segmentLengths(:).')] / brFs;

for k = 1:numel(nsFiles)
    if k > 1
        xline(segmentEdges(k), '--r');
    end

    midTime = (segmentEdges(k) + segmentEdges(k+1)) / 2;
    text( ...
        midTime, 1800, nsFiles{k}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 8, ...
        'Color', 'k', ...
        'BackgroundColor', 'w' ...
    );
end

subplot(3,1,2);
plot(tNoPlot, noSignal, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('NO\_signal (clipped + demeaned, pre-alignment)');
ylim([-2000 2000]);

subplot(3,1,3);
plot(tNoDistPlot, noTargetDist, 'r');
xlabel('Time (s)');
ylabel('Amplitude');
title('NO\_targetDist (pre-alignment)');
ylim([-5 25]);

%% ========================================================================
%  First Cross-Correlation Alignment

fs    = brFs;
brWin = P.br_win;
noWin = P.no_win;

[lagSeconds, noOnBr, overlapInfo] = align_NO_to_BR( ...
    brSignal, noSignal, fs, brWin, noWin, ...
    brTimetable, noTimetable, noResampledMatrix, ...
    signalVars, noTargetDist, brToNoChannel, noToBrChannel, true ...
);

fprintf('Lag = %.3f s\n', lagSeconds);
disp(overlapInfo)

%% ========================================================================
%  Synchronized Master Timetable (NO + BR unified time table)

tCommon = ( ...
    max(brTimetable.Time(1), noOnBr.Time(1)) : ...
    seconds(1 / fs) : ...
    min(brTimetable.Time(end), noOnBr.Time(end)) ...
)';

TT_master = synchronizeWithOptionalSourceFile(brTimetable, noOnBr, tCommon);

figure;
plot(TT_master.Time, TT_master.(brToNoChannel), 'b'); hold on;
plot(TT_master.Time, TT_master.(noToBrChannel), 'r');
legend('BR', 'NO');
grid on;
xlabel('Time');
ylabel('Signal');
title('NO full span with BR zero-padded');
ylim([-1000 1000]);

%% ========================================================================
%  Extract largest time window at the defined target distance window

TT_longest = extractNthLongestWindow( ...
    TT_master, 'targetDistance', fs, targetWindow, frontCut, backCut, 1, brToNoChannel ...
);

if any(strcmp('SourceFile', TT_longest.Properties.VariableNames))
    TT_longest_clean = removevars(TT_longest, {'SourceFile'});
else
    TT_longest_clean = TT_longest;
end

%% ========================================================================
%  2nd Cross Correlation on extracted window

brSegment = TT_longest.(brToNoChannel);
noSegment = TT_longest.(noToBrChannel);

[c, lags] = xcorr(noSegment - mean(noSegment), brSegment - mean(brSegment), 'coeff');
[~, iMax] = max(c);

refinedLagSamples = lags(iMax);
refinedLagSeconds = refinedLagSamples / fs;

fprintf('Refined lag in TT_longest: %d samples (%.6f s)\n', ...
    refinedLagSamples, refinedLagSeconds);

noOnBrRefined = noOnBr;
noOnBrRefined.Time = noOnBr.Time - seconds(refinedLagSeconds);

tCommonRefined = ( ...
    max(brTimetable.Time(1), noOnBrRefined.Time(1)) : ...
    seconds(1 / fs) : ...
    min(brTimetable.Time(end), noOnBrRefined.Time(end)) ...
)';

TT_master_refined = synchronizeWithOptionalSourceFile(brTimetable, noOnBrRefined, tCommonRefined);

TT_longest_refined = extractNthLongestWindow( ...
    TT_master_refined, 'targetDistance', fs, targetWindow, frontCut, backCut, 1 ...
);

figure;
subplot(2,1,1);
plot(lags / fs, c);
grid on;
xlabel('Lag (s)');
ylabel('r');
title('Refined XCorr in TT\_longest');

subplot(2,1,2);
plot(TT_longest.Time, TT_longest.(brToNoChannel), 'b'); hold on;
plot(TT_longest.Time, TT_longest.(noToBrChannel), 'r');
title('Before refined shift');
legend('BR', 'NO');
grid on;
ylim([-1000 1000]);

figure;
plot(TT_master_refined.Time, TT_master_refined.(brToNoChannel), 'b'); hold on;
plot(TT_master_refined.Time, TT_master_refined.(noToBrChannel), 'r');
title('After refined shift (full span)');
legend('BR', 'NO (refined)');
grid on;
ylim([-1000 1000]);

