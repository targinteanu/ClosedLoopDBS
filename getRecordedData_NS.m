function [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS(nsOrFilename, channelIndex)

if nargin < 1
    nsOrFilename = [];
end
if nargin < 2
    channelIndex = [];
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
% The user selects the recording channel. The stimulus trigger channel is
% assumed to be 'ainp1'
channelNames = {ns.ElectrodesInfo.Label}; 
if isempty(channelIndex)
    channelIndex = listdlg('ListString', channelNames);
end
channelName = channelNames{channelIndex};
channelIndexStim = find(contains(channelNames, 'ainp1'));

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
    warning('Packet loss. Largest packet will be selected.') % replace with something better?
    tStartEnd = [t1Rel; t1Rel + datadur];
    packetdur = diff(tStartEnd);
    [~,p] = max(packetdur); 
    dta = ns.Data{p}; datalen = datalen(p); datadur = datadur(p);
    if length(t1Rel) > 1
        t1Rel = t1Rel(p);
    end
else
    dta = ns.Data;
end

% tRel = time relative to start of recording, in seconds 
% t = absolute date/time 
% t0 = date/time at start of recording 
tRel = linspace(0, datadur, datalen) + t1Rel;
t = seconds(tRel);
t0 = datetime(ns.MetaTags.DateTime); 
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
    StimTrainRec = zeros(size(dataOneChannel));
end

end