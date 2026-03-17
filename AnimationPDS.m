%% Phase-Dependent Stimulation Animation 
% This program animates an idealized phase-dependent stimulation (PDS)
% compared with traditional therapeutic deep brain stimulation (DBS). 
% Input signal will be specified by the user, and the program will load a
% specified spreadsheet, matlab .mat, or blackrock .nsX . 
% In a future version, this program may be made into a function that
% accepts the desired signal as input. 

% set these parameters: 

    % signal processing 
    DBSfreq = 130; % Hz
    fbnd = [13, 30]; % frequency bounds [low, high] cutoff (Hz)
    phtarget = 0; % target phase (Hz) 

    % display
    playbackspeed = 1; % relative to real time
    displaywin = 0.5; % seconds 
    packetsize = 10; % sample 
    nbins = 18; % polar histogram (rose) bins

%% Load data 

[fnfe, fp] = uigetfile('*.*'); 
[fp, fn, fe] = fileparts(fullfile(fp,fnfe));
fe = lower(fe);

Fs = nan; % sample rate 
unitname = '';

if contains(fe, '.ns')
    % load blackrock file 
    NS = openNSx(fullfile(fp,fnfe), 'uV');
    T = ns2timetable(NS);
    unitname = ' (\muV)';

elseif strcmp(fe, '.mat')
    % load mat file
    M = load(fullfile(fp,fnfe));
    Mn = fieldnames(M); Mn = lower(Mn);
    MFs = contains(Mn, 'samplerate') | contains(Mn, 'samplingfreq');
    MFs = find(MFs);
    if ~isempty(MFs)
        Fs = M.(Mn{MFs(1)});
    end

    Mtbls = contains(Mn, 'tbl');
    M = struct2cell(M);
    if any(Mtbls)
        % there is at least one variable labeled as table(s)
        M = M(Mtbls); Mn = Mn(Mtbls);
        Msel = cellfun(@(v) iscell(v)|istable(v)|istimetable(v), M);
        M = M(Msel); Mn = Mn(Msel);
        [Msel, selmade] = listdlg("SelectionMode","single", ...
            "PromptString","Select Variable of Interest", ...
            "ListString",Mn);
        if ~selmade
            error('Selection must be made.')
        end
        M = M{Msel}; Mn = Mn{Msel};
        if iscell(M)
            M = M{1}; 
        end
        T = M;

    else
        % first: check for timetables 
        Msel = cellfun(@istimetable, M);
        if ~any(Msel)
            % second: check for any tables 
            Msel = cellfun(@istable, M);
        end
        if ~any(Msel)
            % third: check for any numeric arrays 
            Msel = cellfun(@isnumeric, M);
        end
        if any(Msel)
            M = M(Msel); Mn = Mn(Msel);
            [Msel, selmade] = listdlg("SelectionMode","single", ...
                "PromptString","Select Variable of Interest", ...
                "ListString",Mn);
            if ~selmade
                error('Selection must be made.')
            end
            T = M{Msel};
        else
            error('No suitable variables found in file.')
        end
    end

else
    % load spreadsheet, txt, etc. 
    try
       T = readtimetable(fullfile(fp,fnfe));
    catch ME1
        warning(['Could not load as timetable due to error: ',ME1.message, ...
            '. Loading as table instead.'])
        try
            T = readtable(fullfile(fp,fnfe));
        catch ME2
            warning(['Could not load as table due to error: ',ME2.message, ...
                '. Loading as matrix instead.'])
            try
                X = readmatrix(fullfile(fp,fnfe));
                T = X;
            catch ME3
                rethrow(ME3);
            end
        end
    end
end

%% identify timing of data 

% try to interpret loaded data as timetable: 
if istable(T)
    varnames = T.Properties.VariableNames;
    varnames = string(varnames);
    varnames = lower(varnames);
    timevar = contains(varnames, 'time');
    timevar = find(timevar);
    if ~isempty(timevar)
        timevar = timevar(1);
        timevar = varnames(timevar);
        try
            t = T.(timevar); 
            t = seconds(t);
            T = removevars(T, timevar);
            T = table2timetable(T, 'RowTimes', t);
        catch ME4
            warning("Could not interpret column "+varnames(timevar)+ ...
                " as time due to error: "+ME4.message)
        end
    end
end

% determine original sample rate 
if istimetable(T)
    if isnan(Fs)
        Fs = T.Properties.SampleRate; % Attempt to extract sample rate from timetable
        if isnan(Fs)
            Fs = median(seconds(diff(T.Time)));
        end
    end
else
    Fs = inputdlg('Specify Data Original Sample Rate (Hz):', ...
        'Original Sample Rate', 1, {Fs});
    if isempty(Fs)
        error('Sampling rate must be specified.')
    end
    Fs = str2double(Fs{1});
    if isnan(Fs)
        error('Numeric sampling rate must be specified.')
    end
    t = (1:height(T))/Fs; t = seconds(t)';
    if istable(T)
        T = table2timetable(T, 'RowTimes', t);
    else
        T = array2timetable(T, 'RowTimes', t);
    end
end

% downsample? 

%% Choose channel of interest 

if width(T) > 1
% Prompt user to select the channel of interest
channelNames = T.Properties.VariableNames;
[channelIdx, selMade] = listdlg("SelectionMode","single", ...
    "PromptString","Select Channel of Interest", ...
    "ListString", channelNames);
if ~selMade
    error('Channel must be selected.')
end
else
    channelIdx = 1;
end

x = T{:,channelIdx};
t = T.Time; 
t = seconds(t - t(1));

if isempty(unitname)
    unitnames = string(T.Properties.VariableUnits);
    try 
        unitname = char(unitnames(channelIdx));
        unitname = [' (',unitname,')'];
    catch ME5
        warning(['Could not determine signal units due to error: ',ME5.message])
    end
end

%% identify outliers 
y = log10(movmean((x-mean(x)).^2, 201) + eps); % power estimate 
[clusIdx, clusC] = kmeans(y, 3); 
[~,noiseIdx] = max(clusC);
isnoise = clusIdx == noiseIdx;

%% target phase identification 

% filter 
BPF = buildFIRBPF(Fs, fbnd(1), fbnd(2), 2, 201);
x = filtfilt(BPF,1,x);

% hilbert transform 
H = hilbert(x);
ph = angle(H); A = abs(H);

% set amplitude threshold 
%[~,Athresh] = midcross(A); 
%Athresh = median(A(~isoutlier(x))); 
Athresh = 0.5*median(A(~isnoise)); 

% find where phase crosses target 
phtarget = radfix(phtarget);
ph = ph - phtarget; ph = radfix(ph);
sph = sign(ph);
sph = sph(2:end) .* sph(1:(end-1)); 
dph = diff(ph);
iPDS = (sph <= 0) & (dph >= 0); % rising zero-cross 
iPDS = [false; iPDS];

% finalize PDS stim timing 
iPDS = iPDS & (A >= Athresh);

% DBS stim timing 
DBSper = round(Fs/DBSfreq); % period (samples)
iDBS = false(size(x));
iDBS(1:DBSper:end) = true;
iDBS = iDBS & (A >= Athresh); % make like medtronic closed loop DBS

%% setup display 
% For time plots, plot the entire signal, but set the x axis limits to the
% window of interest. 

displaywin = ceil(displaywin * Fs); % samples 
cursample = 1;
curwin = [0, displaywin-1] + cursample;
%{
winidx = curwin(1):curwin(2);
xnow = x(winidx);
phnow = ph(winidx);
tnow = t(winidx);
iPDSnow = iPDS(winidx);
iDBSnow = iDBS(winidx);
xPDSnow = nan(size(xnow)); xPDSnow(iPDSnow) = xnow(iPDSnow);
xDBSnow = nan(size(xnow)); xDBSnow(iDBSnow) = xnow(iDBSnow);
%}
xPDS = nan(size(x)); xPDS(iPDS) = x(iPDS);
xDBS = nan(size(x)); xDBS(iDBS) = x(iDBS);
nPDS = sum(iPDS(1:curwin(2)));
nDBS = sum(iDBS(1:curwin(2)));

myfig = figure('Units','normalized', 'Position',[.05,.05,.9,.9]);
%tiledlayout(2,3);

% PDS time plot 
if (phtarget > -pi/2) && (phtarget < pi/2)
    mkr = '^';
else
    mkr = 'v';
end
%nexttile([1,2]); 
ax(1,1) = subplot(2,1,1);
plot(t, x, 'LineWidth',1.5);
grid on; hold on; 
stem(t, xPDS, mkr);
xlabel('time (s)'); ylabel(['EPhys',unitname]);
xlim(t(curwin));
title_PDS = title({'Phase Dependent Stimulation'; ...
    ['Total Stim Count = ',num2str(nPDS)]});

%{
% PDS rose plot
bedge = linspace(-pi, pi, nbins);
ax(1,2) = nexttile; 
polarhistogram(phnow(iPDSnow), 'BinEdges',bedge);
title('Stimulation Phase')
%}

% DBS time plot 
%nexttile([1,2]); 
ax(2,1) = subplot(2,1,2);
plot(t, x, 'LineWidth',1.5);
grid on; hold on; 
stem(t, xDBS, 's');
xlabel('time (s)'); ylabel(['EPhys',unitname]);
xlim(t(curwin));
title_DBS = title({'Existing Closed Loop DBS'; ...
    ['Total Stim Count = ',num2str(nDBS)]});

%{
% PDS rose plot
ax(2,2) = nexttile; 
polarhistogram(phnow(iDBSnow), 'BinEdges',bedge);
title('Stimulation Phase')
%}

%% loop through remaining samples 

pausetime = packetsize * playbackspeed / Fs; % seconds 

while cursample + packetsize + displaywin <= length(x)

ticStart = tic;

cursample = cursample + packetsize;
curwin = [0, displaywin-1] + cursample;
%{
winidx = curwin(1):curwin(2);
xnow = x(winidx);
phnow = ph(winidx);
tnow = t(winidx);
iPDSnow = iPDS(winidx);
iDBSnow = iDBS(winidx);
xPDSnow = nan(size(xnow)); xPDSnow(iPDSnow) = xnow(iPDSnow);
xDBSnow = nan(size(xnow)); xDBSnow(iDBSnow) = xnow(iDBSnow);
%}
nPDS = sum(iPDS(1:curwin(2)));
nDBS = sum(iDBS(1:curwin(2)));

%{
plt_X_PDS.XData = tnow; plt_X_PDS.YData = xnow; 
plt_X_DBS.XData = tnow; plt_X_DBS.YData = xnow; 
plt_PDS.XData = tnow; plt_PDS.YData = xPDSnow;
plt_DBS.XData = tnow; plt_DBS.YData = xDBSnow;
%}

% update rose plots to reflect current display window
%rose_PDS = polarhistogram(phnow(iPDSnow), 'BinEdges',bedge);
%rose_DBS = polarhistogram(phnow(iDBSnow), 'BinEdges',bedge);

% advance x axis window 
ax(1,1).XLim = t(curwin);
ax(2,1).XLim = t(curwin);

% update stim count display 
title_PDS.String{2} = ['Total Stim Count = ',num2str(nPDS)];
title_DBS.String{2} = ['Total Stim Count = ',num2str(nDBS)];

ticDur = toc(ticStart);

drawnow; 
if ticDur < pausetime
    pause(pausetime - ticDur);
    drawnow;
end

end