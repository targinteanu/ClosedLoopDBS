function [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS(nsOrFilename, channelIndex, channelIndexStim)

if nargin < 1
    nsOrFilename = [];
end
if nargin < 2
    channelIndex = [];
end
if nargin < 3
    channelIndexStim = [];
end
cl = class(nsOrFilename);

%% Select file and channel with user input

% Access file:  
% Find and load a blackrock ns2 or ns5 file or a mat file with recorded
% brain data. 

if isstruct(nsOrFilename)
    ns = nsOrFilename;
else
    if isempty(nsOrFilename)
        [fn,fp] = uigetfile({'*.ns*'; '*.mat'});
        [~,fn,fe] = fileparts(fn);
    elseif ischar(nsOrFilename) || isstring(nsOrFilename)
        [fp,fn,fe] = fileparts(nsOrFilename);
    end
    if strcmpi(fe,'.mat')
        load(fullfile(fp,[fn,fe]), 'NS2', 'ns2', 'NS5', 'ns5', 'NS', 'ns');
        for vtry = {'NS2', 'ns2', 'NS5', 'ns5', 'NS'}
            if exist(vtry{:}, 'var')
                ns = eval(vtry{:});
            end
        end
    else
        openNSx(fullfile(fp,[fn,fe]), 'uV');
        ns = eval(['NS',fe(end)]);
    end
    clear vtry
end

% User selects channel: 
% If no channels specified, user selects both channels. 
% If only recording channel specified, stim channel defaults to ainp1. 
channelNames = {ns.ElectrodesInfo.Label}; 
if isempty(channelIndex)
    channelIndex = listdlg('ListString', channelNames, ...
        'PromptString', 'Select Recording Channel');
    if isempty(channelIndexStim)
        channelIndexStim = listdlg('ListString', channelNames, ...
            'PromptString', 'Select Stim Trigger Channel (usually ainp1)');
    end
end
if isempty(channelIndexStim)
    channelIndexStim = find(contains(channelNames, 'ainp1'));
end
if isempty(channelIndex)
    channelName = '';
else
    channelName = channelNames{channelIndex};
end

%% interpret data from loaded file 
% obtain usable data variables and other information from the file. 

% Get timing data:
SamplingFreq = ns.MetaTags.SamplingFreq;
datalen = ns.MetaTags.DataPoints;
try
    t1Rel = ns.MetaTags.Timestamp / ns.MetaTags.TimeRes;
catch
    t1Rel = 0;
end
if isfield(ns.MetaTags, 'DataDurationSec')
    datadur = ns.MetaTags.DataDurationSec; 
elseif isfield(ns.MetaTags, 'DataPointsSec')
    datadur = ns.MetaTags.DataPointsSec;
else
    datadur = datalen/SamplingFreq;
end

% handle packet loss 
if (length(t1Rel) > 1) || (length(datadur) > 1)
    warning('Packet loss. Splicing packets together.') 
    tStartEnd = [t1Rel; t1Rel + datadur];
    packetdur = diff(tStartEnd);
    packetgap = tStartEnd(1,2:end) - tStartEnd(2,1:(end-1));
    Dta = ns.Data;
    TRel = arrayfun(@(p) linspace(tStartEnd(1,p), tStartEnd(2,p), datalen(p)), ...
        1:max(length(t1Rel), length(datadur)), 'UniformOutput',false);
    if sum(packetgap < 0)
        % some packets begin before the previous packet has ended. 
        % Solution: trust data from the longer packet
        for p = find(packetgap < 0)
            % packet p and p+1 overlap
            warning(['There is overlap between packets ',num2str(p),' and ',num2str(p+1)]);
            overlapSize = packetgap(p)/SamplingFreq; % samples 
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
    end
    dta = cell2mat(Dta); tRel = cell2mat(TRel);
else
    dta = ns.Data;
    tRel = linspace(0, datadur, datalen) + t1Rel;
end

% tRel = time relative to start of recording, in seconds 
% t = absolute date/time 
% t0 = date/time at start of recording 
t = seconds(tRel);
t0 = datetime(ns.MetaTags.DateTime, 'TimeZone','UTC'); 
t = t+t0; 

% Interpret data from ns structure: 
dataAllChannels = double(dta); 
for ch = 1:height(dataAllChannels) % scale to reported unit - can this be replaced with openNSx('uV') ?
    x = dataAllChannels(ch,:); 
    x = x - double( ns.ElectrodesInfo(ch).MinDigiValue ); 
    x = x .* double( ns.ElectrodesInfo(ch).Resolution );
    x = x + double( ns.ElectrodesInfo(ch).MinAnalogValue );
    dataAllChannels(ch,:) = x;
end
dataOneChannel = dataAllChannels(channelIndex,:);
%StimTrainRec = dataAllChannels(channelIndexStim,:) > 3000;
StimTrainRec = [false, diff(dataAllChannels(channelIndexStim,:)) > 2000]; % rising edge
if ~numel(channelIndexStim)
    StimTrainRec = false(size(dataOneChannel));
end

end