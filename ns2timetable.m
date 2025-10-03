function [TT, packetLoss] = ns2timetable(nsOrFilename)
% Extract data from blackrock ns_ files using the NPMK and output data as a
% timetable TT. Input nsOrFilename can be the data structure obtained by
% running the NPMK (i.e. openNSx) or a filepath to a .ns_ or .mat file; If
% omitted, the user will be prompted to choose a .ns_ or .mat file
% graphically. Also output a flag packetLoss for whether packets needed to
% be stitched together; if so, TT may not be uniformly sampled. 

if nargin < 1
    nsOrFilename = [];
end
cl = class(nsOrFilename);

%% load file

if isstruct(nsOrFilename)
    NS = nsOrFilename;
else
    if isempty(nsOrFilename)
        [fn,fp] = uigetfile({'*.ns*'; '*.mat'});
        [~,fn,fe] = fileparts(fn);
    elseif ischar(nsOrFilename) || isstring(nsOrFilename)
        [fp,fn,fe] = fileparts(nsOrFilename);
    end
    if strcmpi(fe,'.mat')
        load(fullfile(fp,[fn,fe]), 'NS2', 'ns2', 'NS5', 'ns5', 'NS', 'NS');
        for vtry = {'NS2', 'ns2', 'NS5', 'ns5', 'NS'}
            if exist(vtry{:}, 'var')
                NS = eval(vtry{:});
            end
        end
    else
        openNSx(fullfile(fp,[fn,fe]), 'uV');
        NS = eval(['NS',fe(end)]);
    end
    clear vtry
end

%% interpret data from loaded file 
% obtain usable data variables and other information from the file. 

% Get timing data:
SamplingFreq = NS.MetaTags.SamplingFreq;
datalen = NS.MetaTags.DataPoints;
try
    t1Rel = NS.MetaTags.Timestamp / NS.MetaTags.TimeRes;
catch
    t1Rel = 0;
end
if isfield(NS.MetaTags, 'DataDurationSec')
    datadur = NS.MetaTags.DataDurationSec; 
elseif isfield(NS.MetaTags, 'DataPointsSec')
    datadur = NS.MetaTags.DataPointsSec;
else
    datadur = datalen/SamplingFreq;
end

% handle packet loss 
if (length(t1Rel) > 1) || (length(datadur) > 1)
    packetLoss = true;
    warning('Packet loss. Splicing packets together.') 
    [t1Rel, t1sortind] = sort(t1Rel); 
    datadur = datadur(t1sortind); datalen = datalen(t1sortind);
    tStartEnd = [t1Rel; t1Rel + datadur];
    packetdur = diff(tStartEnd);
    packetgap = tStartEnd(1,2:end) - tStartEnd(2,1:(end-1));
    Dta = NS.Data; Dta = Dta(t1sortind);
    TRel = arrayfun(@(p) linspace(tStartEnd(1,p), tStartEnd(2,p), datalen(p)), ...
        1:max(length(t1Rel), length(datadur)), 'UniformOutput',false);

    while sum(packetgap < 0)
        % some packets begin before the previous packet has ended. 
        % Solution: trust data from the longer packet
        for p = find(packetgap < 0)
            % packet p and p+1 overlap
            warning(['There is overlap between packets ',num2str(p),' and ',num2str(p+1)]);
            overlapSize = round(-packetgap(p)*SamplingFreq); % samples 
            if packetdur(p) > packetdur(p+1)
                % shave off start of packet p+1
                Dta{p+1} = Dta{p+1}(:,(overlapSize+1):end);
                TRel{p+1} = TRel{p+1}(:,(overlapSize+1):end);
                warning([num2str(overlapSize),' samples have been removed from packet ',num2str(p+1)]);
            else
                % shave off end of packet p
                Dta{p} = Dta{p}(:,1:(end-overlapSize));
                TRel{p} = TRel{p}(:,1:(end-overlapSize));
                warning([num2str(overlapSize),' samples have been removed from packet ',num2str(p)]);
            end
        end

        % discard packets that are now empty
        todiscard = cellfun(@isempty, TRel);
        TRel = TRel(~todiscard); Dta = Dta(~todiscard);
        tStartEnd = tStartEnd(:,~todiscard);
        packetdur = diff(tStartEnd);
        packetgap = tStartEnd(1,2:end) - tStartEnd(2,1:(end-1));
    end

    dta = cell2mat(Dta); tRel = cell2mat(TRel);
else
    packetLoss = false;
    dta = NS.Data;
    tRel = linspace(0, datadur, datalen) + t1Rel;
end

% tRel = time relative to start of recording, in seconds 
% t = absolute date/time 
% t0 = date/time at start of recording 
t = seconds(tRel);
t0 = datetime(NS.MetaTags.DateTime, 'TimeZone','UTC'); 
t = t+t0; 

lbl = {NS.ElectrodesInfo.Label};
lbl = upper(lbl);
unitname = {NS.ElectrodesInfo.AnalogUnits};

%% correct labels

% limit names to characters that will behave well as plot labels 
for l = 1:length(lbl)
    L = lbl{l}; 
    L = L((L >= 33)&(L <= 126));
    L(L=='_') = ' ';
    lbl{l} = L;
end
for l = 1:length(unitname)
    L = unitname{l}; 
    L = L((L >= 33)&(L <= 126));
    L(L=='_') = ' ';
    unitname{l} = L;
end

% handle duplicate labels 
for l = 1:length(lbl)
    L = lbl{l};
    dups = strcmp(lbl, L);
    dups = find(dups);
    if length(dups) > 1
        warning(['Channel name ',L,' is duplicated.'])
        n = 1;
        for d = dups
            lbl{d} = [L,'(',num2str(n),')'];
            n = n+1;
        end
    end
end

%% construct table

X = dta'; X = double(X);

TT = array2timetable(X,"RowTimes",t,"VariableNames",lbl); 
TT.Properties.VariableUnits = unitname;
%TT.Properties.SampleRate = NS.MetaTags.SamplingFreq;

TT.Properties.Description = [...
    'file ',NS.MetaTags.FilePath,filesep,NS.MetaTags.Filename,NS.MetaTags.FileExt,...
    ', comment ',NS.MetaTags.Comment];

% sample rate check 
fs_Signal = NS.MetaTags.SamplingFreq;
dT = seconds(diff(t));
dTmax = max(dT); dTmin = min(dT);
fsmax = 1/dTmin; fsmin = 1/dTmax;
errthresh = .01;
if abs((fsmax-fs_Signal)/fs_Signal) > errthresh
    warning('Sample Rate is incorrectly reported or inconsistent.')
end
if abs((fsmin-fs_Signal)/fs_Signal) > errthresh
    warning('Sample Rate is incorrectly reported or inconsistent.')
end

end