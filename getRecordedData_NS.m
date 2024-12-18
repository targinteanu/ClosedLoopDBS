function [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
    channelName, channelIndex, channelIndexStim, channelNames]...
    = getRecordedData_NS(nsOrFilename)

if nargin < 1
    nsOrFilename = [];
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
        openNSx(fullfile(fp,[fn,fe]));
        ns = eval(['NS',fe(end)]);
    end
    clear vtry
end

% User selects channel: 
% The user selects the recording channel. The stimulus trigger channel is
% assumed to be 'ainp1'
channelNames = {ns.ElectrodesInfo.Label}; 
channelIndex = listdlg('ListString', channelNames);
channelName = channelNames{channelIndex};
channelIndexStim = find(contains(channelNames, 'ainp1'));

%% interpret data from loaded file 
% obtain usable data variables and other information from the file. 

% Get timing data:
% tRel = time relative to start of recording, in seconds 
% t = absolute date/time 
% t0 = date/time at start of recording 
SamplingFreq = ns.MetaTags.SamplingFreq;
try
    tRel = linspace(0,ns.MetaTags.DataPointsSec,ns.MetaTags.DataPoints);
catch
    try
        tRel = linspace(0,ns.MetaTags.DataDurationSec,ns.MetaTags.DataPoints);
    catch
        tRel = linspace(0,ns.MetaTags.DataPoints/SamplingFreq,ns.MetaTags.DataPoints);
    end
end
t = seconds(tRel);
t0 = datetime(ns.MetaTags.DateTime); 
t = t+t0; 

% Interpret data from ns structure: 
dataAllChannels = double(ns.Data); 
for ch = 1:height(dataAllChannels) % scale to reported unit 
    x = dataAllChannels(ch,:); 
    x = x - double( ns.ElectrodesInfo(ch).MinDigiValue ); 
    x = x .* double( ns.ElectrodesInfo(ch).Resolution );
    x = x + double( ns.ElectrodesInfo(ch).MinAnalogValue );
    dataAllChannels(ch,:) = x;
end
dataOneChannel = dataAllChannels(channelIndex,:);
StimTrainRec = dataAllChannels(channelIndexStim,:) > 1e4;
if ~numel(channelIndexStim)
    StimTrainRec = zeros(size(dataOneChannel));
end

end