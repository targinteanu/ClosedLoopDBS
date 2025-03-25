function [emptyData, contData, buffData, chanInfo, startTic] = ...
    initRawData_AO(chsel, bufferSize)
% Initialize the multichannel raw data structure using Alpha Omega
% interface. 
%
% Inputs: chsel is horizontal selected channel ID numbers output from AO.
% bufferSize is the size(s) of each corresponding channel buffer (samples). 
% 
% Output data structure is a cell array with columns to each channel. Each
% channel is represented by the matrix [time column, data column]. Time
% column contains time stamps in machine time seconds and nan elsewhere. 
% 
% chanInfo has fields: SampleRate, Name, Unit, IDnumber

%% connect to hardware 
connect_AO(); startTic = tic; 
pause(1);

%% setup

[Results, channelsData] = AO_GetAllChannels();
%%{
if Results == 4
    % give it some time and try again 
    disconnect_AO(); 
    pause(1); 
    connect_AO();
    pause(1);
    [Results, channelsData] = AO_GetAllChannels();
    if Results
        error('tried again and still failed to connect')
    end
end
%}
if Results
    error(['Failed to obtain channel info with error code ',num2str(Results)])
end

chnum = [channelsData.channelID]; chnum = double(chnum);
chname = {channelsData.channelName};

% limit to only desired continuous data channels
chincl = false(size(chnum));
for ch = 1:length(chnum)
    chname_ch = chname{ch};
    % 
    % Currently including only LFP channels because they have the same
    % sample rate. TO DO: include other channels with their sample rates 
    %
    % currently excluding SPK (spike?), FE/RM Audio, RM AO, UD, InPort,
    % StimMarker, Internal Detection, TTL, DOUT, TS Sync, ACC_X/Y/Z
    % (acceleration?)
    %
    % Channel names on Neuro Omega need to be included here!!! (unsure if they
    % are the same)
    % 
    if contains(chname_ch, 'LFP')
        n = sscanf(chname_ch, 'LFP %f'); 
        chincl(ch) = n < 97;
        %{
        %{
    elseif contains(chname_ch, 'SEG')
        %SEEG (?)
        n = sscanf(chname_ch, 'SEG %f'); 
        chincl(ch) = n < 97;
        %}
    elseif contains(chname_ch, 'RAW')
        % ?
        n = sscanf(chname_ch, 'RAW %f'); 
        chincl(ch) = n < 11;
    elseif contains(chname_ch, 'AI')
        % Analog In (?)
        n = sscanf(chname_ch, 'AI %f'); 
        chincl(ch) = n < 9;
        %}
    end
    % 
    % TO DO: can above chan inclusion be automated, i.e. try 
    % [Result, contData, DataCapture] = AO_GetChannelData(chsel(ch)) and
    % exclude channels with nonzero error / nan DataCapture? 
    % 
end
%{
chincl = ...
    contains(chname, 'LFP') | ...
    contains(chname, 'SEG') | ... SEEG (?)
    contains(chname, 'RAW') | ... ?
    contains(chname, 'AI') ; % Analog In (?)
%}
chnum = chnum(chincl); chname = chname(chincl);

if isempty(chsel)
    % select all channels
    chsel = chnum;
end
chsel = double(chsel);

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
bufferSizeAO = max(bufferSizeAO, 5000);  % min allowed 
bufferSizeAO = min(bufferSizeAO, 20000); % max allowed
ChannelGain = 20;
BitResolution = 2500000/(2^16*ChannelGain);

%% initialize all channels

% add buffering channels 
W = length(chsel);
for ch = 1:W
    chnum_ch = chsel(ch); chname_ch = chname{ch};
    Results = AO_AddBufferingChannel(chnum_ch, bufferSizeAO);
    if Results
        msg = [...
            'Failed to initiate channel ',num2str(chnum_ch),': ',chname_ch,' ',...
            'with error code ',num2str(Results)];
        error(msg); % consider changing to a warning and excluding this channel.
    end
end

%we clear the old data so we can have the new data
AO_ClearChannelData(); 
pause(1);

% get initial data 
[Results,continuousData,DataCapture,time] = AO_GetAlignedData(chsel); 
if Results == -3
    % no samples; give it more time 
    pause(1); 
    [Results,continuousData,DataCapture,time] = AO_GetAlignedData(chsel);
end
if Results
    msg = ['Failed to acquire data with error code ',num2str(Results)];
    error(msg); % consider trying again until success 
end
time = time/1510; % TO DO: make sure this is correct time in seconds!
continuousData = continuousData(1:DataCapture); 
L = DataCapture/W;
continuousData = reshape(continuousData, ...
    L, W); % columns = channels 
continuousData = double(continuousData);

%% assign data to the structure 

for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch));

    if ~isempty(chInd)
        if length(chInd) > 1
            error('Non-unique channel ID(s).')
        end
        chnum_ch = chsel(ch);
        chname_ch = chname{chInd};

        % Create raw data buffer of zeros of the correct length
        emptyData{ch} = [nan(bufferSize(ch),1), zeros(bufferSize(ch),1)];
        %emptyData{ch}(1,1) = time - (bufferSize)/fs(chInd);
        emptyData{ch}(end,1) = time - 1/fs; % ??
        contData{ch} = [nan(L,1), continuousData(:,chInd)];
        contData{ch}(1,1) = time;

        % channel info 
        ud.SampleRate = fs; 
        ud.Name = chname_ch;
        ud.Unit = '';
        ud.MinDigital = 0; % to negate offset calculation
        ud.MinAnalog = 0;  % to negate offset calculation
        ud.Resolution = BitResolution;
        ud.IDnumber = chnum_ch;
        chanInfo{ch} = ud; 

        buffData{ch} = bufferData(emptyData{ch}, contData{ch});

    end

end

end