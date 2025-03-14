function [emptyData, contData, buffData, chanInfo, startTic] = ...
    initRawData_AO(chsel, bufferSize)
% Initialize the multichannel raw data structure using Alpha Omega
% interface. 
%
% Inputs: chsel is selected channel ID numbers output from BlackRock.
% bufferSize is the size(s) of each corresponding channel buffer (samples). 
% 
% Output data structure is a cell array with columns to each channel. Each
% channel is represented by the matrix [time column, data column]. Time
% column contains time stamps in machine time seconds and nan elsewhere. 
% 
% chanInfo has fields: SampleRate, Name, Unit, IDnumber

%% setup

[Results, channelsData] = AO_GetAllChannels();
if Results
    error(['Failed to obtain channel info with error code ',num2str(Results)])
end

chnum = [channelsData.channelID];
chname = {channelsData.channelName};

% limit to only desired continuous data channels
chincl = ...
    contains(chname, 'LFP') | ...
    contains(chname, 'SEG') | ... SEEG (?)
    contains(chname, 'RAW') | ... ?
    contains(chname, 'AI') ; % Analog In (?)
% currently excluding SPK (spike?), FE/RM Audio, RM AO, UD, InPort,
% StimMarker, Internal Detection, TTL, DOUT, TS Sync, ACC_X/Y/Z
% (acceleration?)
% Channel names on Neuro Omega need to be included here!!! (unsure if they 
% are the same)
chnum = chnum(chincl); chname = chname(chincl);

if isempty(chsel)
    % select all channels
    chsel = chnum;
end

emptyData = cell(1,length(chsel)); 
contData = emptyData; 
buffData = emptyData;
chanInfo = emptyData;

%% handle/check inputs 
if length(bufferSize) < length(chsel)
    if length(bufferSize) == 1
        % assume the one input applies to all channels.
        bufferSize = repmat(bufferSize, size(chsel));
    else
        error('Incompatible input dimensions.')
    end
end

%% set (desired) AO params here
fs = 22000; % sampling rate, Hz 
bufferSizeAO = ceil(1.1*1000*max(bufferSize)/fs); % ms

%% assign data to the structure 

for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch));

    if ~isempty(chInd)
        if length(chInd) > 1
            error('Non-unique channel ID(s).')
        end
        chnum_ch = chsel(ch);
        chname_ch = chname{chInd};

        % initiate the buffering channel 
        Results = AO_AddBufferingChannel(chnum_ch, bufferSizeAO);
        if Results
            msg = [...
                'Failed to initiate channel ',num2str(chnum_ch),': ',chname_ch,' ',...
                'with error code ',num2str(Results)];
            error(msg); % consider changing to a warning and excluding this channel.
        end

        % Create raw data buffer of zeros of the correct length
        emptyData{ch} = [nan(bufferSize(ch),1), zeros(bufferSize(ch),1)];

    end
end

end