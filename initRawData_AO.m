function [emptyData, contData, buffData, chanInfo, startTic] = ...
    initRawData_AO(chsel, bufferSize, chantype)
% Initialize the multichannel raw data structure using Alpha Omega
% interface. 
%
% Inputs: chsel is horizontal selected channel ID numbers output from AO.
% bufferSize is the size(s) of each corresponding channel buffer (samples). 
% chantype is e.g. 'LFP', 'AI' (analog in), 'SPK', 'RAW', 'SEG'. Currently,
% only one type at a time is supported. 
% 
% Output data structure is a cell array with columns to each channel. Each
% channel is represented by the matrix [time column, data column]. Time
% column contains time stamps in machine time seconds and nan elsewhere. 
% 
% chanInfo has fields: SampleRate, Name, Unit, IDnumber

if nargin < 3
    chantype = 'LFP';
end

%% hardcode known channel info according to type 
% There does not seem to be a way to extract this info automatically.
%
% currently excluding FE/RM Audio, RM AO, UD, InPort,
% StimMarker, Internal Detection, TTL, DOUT, TS Sync, ACC_X/Y/Z
% (acceleration?)
%
% Channel names on Neuro Omega need to be included here!!! (unsure if they
% are the same)
% 
if strcmpi(chantype, 'LFP')
    % filtered/downsampled LFP
    fs = 1875;
    channelselector = @(n) (n < 97) && (n >= 33);
elseif strcmpi(chantype, 'AI')
    % analog in
    fs = 3750;
    channelselector = @(n) n < 9;
elseif strcmpi(chantype, 'SEG')
    % timestamp of threshold crossing (spike sorting)
    fs = 30000;
    channelselector = @(n) n < 97;
elseif strcmpi(chantype, 'SPK')
    % filtered spike
    fs = 30000;
    channelselector = @(n) false;
elseif strcmpi(chantype, 'RAW')
    % unprocessed (no filering/downsampling) ver of other signals
    fs = 30000;
    channelselector = @(n) n < 11;
else
    fs = nan; 
    channelselector = @(n) false;
end

%% connect to hardware 
connect_AO(); startTic = tic; 
pause(1);

%% setup

chnamefspec = [chantype,' %f'];

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
chincl = false(size(chnum)); Fs = nan(size(chnum));
for chInd = 1:length(chnum)
    chname_ch = chname{chInd};
    if contains(chname_ch, chantype)
        n = sscanf(chname_ch, chnamefspec);
        chincl(chInd) = true; % alternatively, use chanselector(n)
        Fs(chInd) = fs;
    end
    % 
    % TO DO: can above chan inclusion be automated, i.e. try 
    % [Result, contData, DataCapture] = AO_GetChannelData(chsel(ch)) and
    % exclude channels with nonzero error / nan DataCapture? 
    % 
end

chnum = chnum(chincl); chname = chname(chincl); Fs = Fs(chincl);

if isempty(chsel)
    % select all channels
    chsel = chnum;
end
chsel = double(chsel);

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
bufferSizeAO = ceil(1.1*1000*max(bufferSize)/min(Fs)); % ms
bufferSizeAO = max(bufferSizeAO, 5000);  % min allowed 
bufferSizeAO = min(bufferSizeAO, 20000); % max allowed
ChannelGain = 20;
BitResolution = 2500000/(2^16*ChannelGain);

%% initialize all channels

% add buffering channels 
W = length(chsel);
for chInd = 1:W
    chnum_ch = chsel(chInd); chname_ch = chname{chInd};
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
doretry = true;
while doretry
[Results,continuousData,DataCapture,time] = AO_GetAlignedData(chsel); 
if Results == -3
    % no samples; give it more time 
    pause(.1); 
    [Results,continuousData,DataCapture,time] = AO_GetAlignedData(chsel);
end
if Results
    [~,~,LastError] = AO_GetError();
    [errchan, nerr] = sscanf(LastError, 'ERROR --> AO_GetAlignedData :: Channel %f has no samples');
    if nerr > 0
        warning(['Removing channel ',num2str(errchan),' due to error code ',num2str(Results)])
        remchan = (chsel == errchan);
        chsel = chsel(~remchan);
        bufferSize = bufferSize(~remchan);
        doretry = true;
    else
    msg = ['Failed to acquire data with error code ',num2str(Results)];
    error(msg); % consider trying again until success 
    end
else
    doretry = false;
end
end

W = length(chsel);

time = time/1510; % TO DO: make sure this is correct time in seconds!
continuousData = continuousData(1:DataCapture); 
L = DataCapture/W;
continuousData = reshape(continuousData, ...
    L, W); % columns = channels 
continuousData = double(continuousData);

%% assign data to the structure 

emptyData = cell(1,length(chsel)); 
contData = emptyData; 
buffData = emptyData;
chanInfo = emptyData;

for ch = 1:length(chnum)
    chnum_ch = chnum(ch);
    chInd = find(chsel == chnum_ch);

    if ~isempty(chInd)
        if length(chInd) > 1
            error('Non-unique channel ID(s).')
        end
        chname_ch = chname{ch};
        fs = Fs(ch);

        % Create raw data buffer of zeros of the correct length
        emptyData{chInd} = [nan(bufferSize(chInd),1), zeros(bufferSize(chInd),1)];
        if ~isempty(emptyData{chInd})
            %emptyData{ch}(1,1) = time - (bufferSize)/fs(chInd);
            emptyData{chInd}(end,1) = time - 1/fs; % ??
        end
        contData{chInd} = [nan(L,1), continuousData(:,chInd)];
        contData{chInd}(1,1) = time;

        % channel info 
        ud.SampleRate = fs; 
        ud.Name = chname_ch;
        ud.Unit = '';
        ud.MinDigital = 0; % to negate offset calculation
        ud.MinAnalog = 0;  % to negate offset calculation
        ud.Resolution = BitResolution;
        ud.IDnumber = chnum_ch;
        chanInfo{chInd} = ud; 

        buffData{chInd} = bufferData(emptyData{chInd}, contData{chInd});

    end

end

end